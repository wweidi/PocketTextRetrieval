#!/usr/bin/env python3
"""Build pocket-level PDB files and manifests from BioLiP chain structures.

The input chain files are downloaded receptor structures.  For each BioLiP
row this script selects the annotated binding-site residues using PDB residue
numbering, writes a pocket-only PDB, and adds the resulting path and quality
metadata to a new manifest.  Original manifests and raw chain files are not
modified.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


RESIDUE_TOKEN = re.compile(r"^([A-Za-z]{1,3})(-?\d+)([A-Za-z]?)$")


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def residue_keys(text: str) -> set[tuple[int, str]]:
    result: set[tuple[int, str]] = set()
    for token in re.split(r"[\s,;]+", text or ""):
        match = RESIDUE_TOKEN.match(token.strip())
        if match:
            result.add((int(match.group(2)), match.group(3).upper()))
    return result


def pdb_residue_key(line: str) -> tuple[int, str] | None:
    if len(line) < 27:
        return None
    number = line[22:26].strip()
    if not number:
        return None
    try:
        return int(number), line[26:27].strip().upper()
    except ValueError:
        return None


def load_chain(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n\r") for line in handle if line.startswith("ATOM  ")]


def select_atoms(lines: list[str], wanted: set[tuple[int, str]]) -> list[str]:
    return [line for line in lines if (key := pdb_residue_key(line)) in wanted]


def write_pocket(path: Path, atoms: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for line in atoms:
            handle.write(line + "\n")
        handle.write("TER\n")


def fallback_text(row: dict[str, str]) -> str:
    ligand = (row.get("ligand_name") or row.get("ligand_id") or "unknown ligand").strip()
    count = len(residue_keys(row.get("binding_residues_pdb", "")))
    if not count:
        count = int(row.get("residue_count") or 0)
    catalytic = " The annotated site also contains catalytic residues." if row.get("catalytic_residues", "").strip() else ""
    return (
        f"A protein binding pocket for {ligand} with {count} observed contact residues."
        f" The pocket is defined by the experimentally annotated local binding site."
        f"{catalytic}"
    )


def process_manifest(
    input_path: Path,
    output_path: Path,
    raw_dir: Path,
    pocket_dir: Path,
    log_handle,
    counters: Counter[str],
) -> None:
    with input_path.open("r", encoding="utf-8", errors="replace", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        extra_fields = [
            "pocket_structure_path",
            "pocket_atom_count",
            "pocket_residue_count",
            "pocket_numbering_source",
            "text_semantic_status",
        ]
        for field in extra_fields:
            if field not in fields:
                fields.append(field)
        rows = list(reader)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        total = len(rows)
        for index, row in enumerate(rows, start=1):
            counters["rows"] += 1
            pdb_id = (row.get("pdb_id") or "").strip().lower()
            chain = (row.get("receptor_chain") or "").strip()
            pocket_id = (row.get("pocket_id") or f"{pdb_id}_{chain}_{index}").strip()
            chain_path = raw_dir / f"{pdb_id}{chain}.pdb"
            pocket_path = pocket_dir / f"{pocket_id}.pdb"
            status = "ok"
            numbering_source = "pdb"
            atom_lines: list[str] = []

            if not chain_path.exists():
                status = "missing_chain_file"
            else:
                lines = load_chain(chain_path)
                pdb_wanted = residue_keys(row.get("binding_residues_pdb", ""))
                atom_lines = select_atoms(lines, pdb_wanted)
                if not atom_lines:
                    renumbered_wanted = residue_keys(row.get("binding_residues_renumbered", ""))
                    atom_lines = select_atoms(lines, renumbered_wanted)
                    if atom_lines:
                        numbering_source = "renumbered_fallback"
                if not atom_lines:
                    status = "no_matching_residues"

            if atom_lines:
                write_pocket(pocket_path, atom_lines)
                counters["pockets"] += 1
                counters["atoms"] += len(atom_lines)
                counters[f"numbering_{numbering_source}"] += 1
            else:
                counters["failed"] += 1

            text_value = (row.get("text_semantic") or "").strip()
            if text_value:
                text_status = "existing"
                counters["text_existing"] += 1
            else:
                text_value = fallback_text(row)
                text_status = "generated_fallback"
                counters["text_generated"] += 1
            row["text_semantic"] = text_value
            row["pocket_structure_path"] = str(pocket_path)
            row["pocket_atom_count"] = str(len(atom_lines))
            row["pocket_residue_count"] = str(len({pdb_residue_key(line) for line in atom_lines if pdb_residue_key(line)}))
            row["pocket_numbering_source"] = numbering_source if atom_lines else "none"
            row["text_semantic_status"] = text_status
            writer.writerow(row)

            log_handle.write(
                f"{now()} split={input_path.stem} index={index}/{total} status={status} "
                f"pocket_id={pocket_id} atoms={len(atom_lines)} "
                f"residues={row['pocket_residue_count']} numbering={row['pocket_numbering_source']} "
                f"text={text_status}\n"
            )
            if index % 100 == 0 or index == total:
                log_handle.flush()
                print(
                    f"progress split={input_path.stem} {index}/{total} "
                    f"pockets={counters['pockets']} failed={counters['failed']}",
                    flush=True,
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--val", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--pocket-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()

    args.pocket_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    counters: Counter[str] = Counter()
    with args.log.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"{now()} event=start\n")
        process_manifest(
            args.train,
            args.output_dir / "train_pocket_manifest.tsv",
            args.raw_dir,
            args.pocket_dir,
            log_handle,
            counters,
        )
        process_manifest(
            args.val,
            args.output_dir / "val_pocket_manifest.tsv",
            args.raw_dir,
            args.pocket_dir,
            log_handle,
            counters,
        )
        process_manifest(
            args.test,
            args.output_dir / "test_biolip_pocket_manifest.tsv",
            args.raw_dir,
            args.pocket_dir,
            log_handle,
            counters,
        )
        log_handle.write(f"{now()} event=summary counts={dict(counters)}\n")
        log_handle.flush()
    print(f"summary={dict(counters)}")
    if counters["failed"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
