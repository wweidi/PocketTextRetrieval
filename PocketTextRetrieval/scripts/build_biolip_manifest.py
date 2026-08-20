#!/usr/bin/env python3
"""Build an initial pocket-text manifest from BioLiP annotations.

The generated semantic text is intentionally local to the binding pocket. It
does not use the protein-level UniProt function text. Exact residue labels are
kept only in text_full for an ablation, not in the main semantic text.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path


HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQCY")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
AROMATIC = set("FWY")


def load_ligands(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ccd = (row.get("#CCD") or row.get("CCD") or "").strip()
            if ccd:
                result[ccd] = row
    return result


def load_excluded_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    if not path.exists():
        return ids
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            value = line.strip().split(",", 1)[0].strip()
            if value and value.lower() not in {"id", "protein_id"}:
                ids.add(value)
    return ids


def residue_summary(residue_text: str) -> tuple[int, Counter[str], str]:
    tokens = re.findall(r"([A-Za-z])(-?\d+)", residue_text or "")
    residues = [letter.upper() for letter, _ in tokens if letter.upper() in set("ACDEFGHIKLMNPQRSTVWY")]
    counts = Counter(residues)
    n = len(residues)
    if not n:
        return 0, counts, ""
    groups = []
    if sum(counts[a] for a in HYDROPHOBIC) / n >= 0.35:
        groups.append("hydrophobic")
    if sum(counts[a] for a in POLAR) / n >= 0.20:
        groups.append("polar")
    if sum(counts[a] for a in POSITIVE) / n >= 0.12:
        groups.append("positively charged")
    if sum(counts[a] for a in NEGATIVE) / n >= 0.12:
        groups.append("negatively charged")
    if sum(counts[a] for a in AROMATIC) / n >= 0.12:
        groups.append("aromatic")
    return n, counts, ", ".join(groups) if groups else "mixed-composition"


def ligand_class(ccd: str, name: str) -> str:
    text = f"{ccd} {name}".lower()
    if any(x in text for x in ("atp", "adp", "amp", "gdp", "gtp", "nad", "nucleotide")):
        return "nucleotide-like ligand"
    if any(x in text for x in ("heme", "porphyrin", "chlorophyll")):
        return "porphyrin-like ligand"
    if any(x in text for x in ("cofactor", "flavin", "fmn", "fad")):
        return "cofactor-like ligand"
    if any(x in text for x in ("sugar", "glucose", "fructose", "saccharide")):
        return "carbohydrate-like ligand"
    if any(x in text for x in ("peptide", "leu", "gly", "ala")) and len(name) > 20:
        return "peptide-like ligand"
    return "small-molecule ligand"


def build_text(ccd: str, name: str, residue_text: str, catalytic_text: str) -> tuple[str, str, int]:
    count, counts, composition = residue_summary(residue_text)
    display_name = name.strip() if name.strip() else ccd
    display_name = re.sub(r"\s+", " ", display_name).split(";", 1)[0][:160]
    category = ligand_class(ccd, display_name)
    interaction_hint = ""
    if catalytic_text.strip():
        interaction_hint = " The annotated site also contains catalytic residues."
    semantic = (
        f"A protein binding pocket for {display_name} ({category}) with {count} observed "
        f"contact residues. The pocket has a {composition} local chemical environment."
        f"{interaction_hint}"
    )
    full = semantic
    if residue_text.strip():
        full += f" Contact residues include {residue_text.strip()}."
    return semantic, full, count


def split_for_group(group: str) -> str:
    digest = hashlib.sha1(group.encode("utf-8")).hexdigest()
    value = int(digest[:8], 16) / 0xFFFFFFFF
    if value < 0.80:
        return "train"
    if value < 0.90:
        return "val"
    return "test"


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, value: str) -> None:
        self.parent.setdefault(value, value)

    def find(self, value: str) -> str:
        self.add(value)
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--biolip", type=Path, required=True)
    parser.add_argument("--ligands", type=Path, required=True)
    parser.add_argument("--oneprot-test-pocket", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ligands = load_ligands(args.ligands)
    excluded = load_excluded_ids(args.oneprot_test_pocket)
    fields = [
        "pocket_id", "pdb_id", "receptor_chain", "resolution", "site_code",
        "ligand_id", "ligand_chain", "ligand_serial", "uniprot_id",
        "binding_residues_pdb", "binding_residues_renumbered", "catalytic_residues",
        "ec_number", "go_terms", "affinity_literature", "pubmed_id",
        "sequence", "ligand_name", "ligand_formula", "ligand_smiles",
        "residue_count", "text_semantic", "text_full", "protein_group",
        "split", "structure_url", "structure_path", "source",
    ]
    counts = Counter()
    records = []
    union_find = UnionFind()
    with gzip.open(args.biolip, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 21:
                counts["malformed"] += 1
                continue
            (
                pdb_id, receptor_chain, resolution, site_code, ligand_id,
                ligand_chain, ligand_serial, binding_pdb, binding_seq,
                catalytic_pdb, catalytic_seq, ec_number, go_terms,
                affinity_lit, affinity_moad, affinity_pdbbind, affinity_db,
                uniprot_id, pubmed_id, ligand_seq, sequence,
            ) = columns[:21]
            ligand_lower = ligand_id.lower()
            if not uniprot_id or not binding_pdb or not sequence:
                counts["missing_required"] += 1
                continue
            if ligand_lower in {"rna", "dna", "pep", "peptide", "hoh", "wat"}:
                counts["non_small_molecule"] += 1
                continue
            try:
                resolution_value = float(resolution)
            except ValueError:
                counts["bad_resolution"] += 1
                continue
            if resolution_value <= 0 or resolution_value > 3.5:
                counts["resolution_filter"] += 1
                continue
            record = {
                "pdb_id": pdb_id,
                "receptor_chain": receptor_chain,
                "resolution": resolution,
                "site_code": site_code,
                "ligand_id": ligand_id,
                "ligand_chain": ligand_chain,
                "ligand_serial": ligand_serial,
                "uniprot_id": uniprot_id,
                "binding_pdb": binding_pdb,
                "binding_seq": binding_seq,
                "catalytic_pdb": catalytic_pdb,
                "ec_number": ec_number,
                "go_terms": go_terms,
                "affinity_lit": affinity_lit,
                "pubmed_id": pubmed_id,
                "sequence": sequence,
            }
            records.append(record)
            union_find.union(f"u:{uniprot_id}", f"p:{pdb_id.lower()}")

    component_members: dict[str, set[str]] = {}
    for node in union_find.parent:
        component_members.setdefault(union_find.find(node), set()).add(node)
    excluded_components = {
        root for root, members in component_members.items()
        if any(node.startswith("u:") and node[2:] in excluded for node in members)
    }

    output_handles = {}
    writers = {}
    try:
        for name in ("all", "train", "val", "test"):
            handle = (args.output_dir / f"biolip_nr_pockettext_{name}.tsv").open(
                "w", encoding="utf-8", newline=""
            )
            output_handles[name] = handle
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writers[name] = writer

        for record in records:
            component_root = union_find.find(f"u:{record['uniprot_id']}")
            if component_root in excluded_components:
                counts["oneprot_test_excluded"] += 1
                continue
            ligand_id = record["ligand_id"]
            ligand = ligands.get(ligand_id, {})
            ligand_name = ligand.get("name", "")
            semantic, full, residue_count = build_text(
                ligand_id, ligand_name, record["binding_pdb"], record["catalytic_pdb"]
            )
            pocket_id = "{}_{}_{}_{}_{}_{}".format(
                record["pdb_id"].lower(), record["receptor_chain"], record["site_code"],
                ligand_id, record["ligand_chain"], record["ligand_serial"]
            )
            protein_group = f"component:{component_root}"
            split = split_for_group(protein_group)
            row = {
                "pocket_id": pocket_id,
                "pdb_id": record["pdb_id"].lower(),
                "receptor_chain": record["receptor_chain"],
                "resolution": record["resolution"],
                "site_code": record["site_code"],
                "ligand_id": ligand_id,
                "ligand_chain": record["ligand_chain"],
                "ligand_serial": record["ligand_serial"],
                "uniprot_id": record["uniprot_id"],
                "binding_residues_pdb": record["binding_pdb"],
                "binding_residues_renumbered": record["binding_seq"],
                "catalytic_residues": record["catalytic_pdb"],
                "ec_number": record["ec_number"],
                "go_terms": record["go_terms"],
                "affinity_literature": record["affinity_lit"],
                "pubmed_id": record["pubmed_id"],
                "sequence": record["sequence"],
                "ligand_name": ligand_name,
                "ligand_formula": ligand.get("formula", ""),
                "ligand_smiles": ligand.get("SMILES", ""),
                "residue_count": str(residue_count),
                "text_semantic": semantic,
                "text_full": full,
                "protein_group": protein_group,
                "split": split,
                "structure_url": f"https://files.rcsb.org/download/{record['pdb_id'].upper()}.pdb",
                "structure_path": f"data/raw/biolip/receptor_nr/{record['pdb_id'].lower()}{record['receptor_chain']}.pdb",
                "source": "BioLiP_nr",
            }
            writers["all"].writerow(row)
            writers[split].writerow(row)
            counts["kept"] += 1
    finally:
        for handle in output_handles.values():
            handle.close()
    print("\n".join(f"{key}={value}" for key, value in sorted(counts.items())))


if __name__ == "__main__":
    main()
