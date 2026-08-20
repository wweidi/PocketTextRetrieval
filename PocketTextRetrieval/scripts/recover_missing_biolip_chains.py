#!/usr/bin/env python3
"""Recover BioLiP receptor chains from RCSB mmCIF files.

BioLiP's receptor batch archives use chain identifiers such as ``AAA`` or
``SPE1`` that cannot always be represented losslessly in legacy PDB columns.
The recovery step therefore matches the exact mmCIF ``auth_asym_id`` and
writes a single-chain PDB-like file in the same directory as the BioLiP
chain files.  The output is intentionally limited to ATOM records, matching
the existing BioLiP receptor files.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import shlex
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_missing(path: Path) -> dict[str, set[str]]:
    targets: dict[str, set[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key = line.strip().lower()
        if len(key) < 5:
            continue
        # PDB ids are normalized to lowercase for filenames; mmCIF chain
        # identifiers are matched case-insensitively but written in the
        # conventional uppercase form used by the BioLiP files.
        targets.setdefault(key[:4], set()).add(key[4:].upper())
    return targets


def fetch_cif(pdb_id: str, cache_dir: Path, retries: int) -> tuple[str, Path | None, str]:
    target = cache_dir / f"{pdb_id}.cif"
    if target.exists() and target.stat().st_size > 100:
        return pdb_id, target, "cached"
    temporary = cache_dir / f".{pdb_id}.cif.part"
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    request = urllib.request.Request(url, headers={"User-Agent": "PocketTextRetrieval/0.1"})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            if len(payload) <= 100:
                raise OSError("RCSB response is unexpectedly small")
            temporary.write_bytes(payload)
            temporary.replace(target)
            return pdb_id, target, "downloaded"
        except (OSError, urllib.error.URLError, TimeoutError) as exc:
            if attempt >= retries:
                try:
                    temporary.unlink()
                except FileNotFoundError:
                    pass
                return pdb_id, None, f"failed:{type(exc).__name__}"
            time.sleep(2**attempt)
    return pdb_id, None, "failed:unknown"


def atom_rows(cif_path: Path) -> list[dict[str, str]]:
    lines = cif_path.read_text(encoding="utf-8", errors="replace").splitlines()
    rows: list[dict[str, str]] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        header_start = index + 1
        headers: list[str] = []
        while header_start < len(lines) and lines[header_start].lstrip().startswith("_"):
            headers.append(lines[header_start].strip().split()[0])
            header_start += 1
        if "_atom_site.group_PDB" not in headers:
            index = header_start
            continue
        index = header_start
        while index < len(lines):
            line = lines[index].strip()
            if not line or line.startswith("#"):
                index += 1
                if line.startswith("#"):
                    break
                continue
            if line == "loop_" or line.startswith("_"):
                break
            fields = shlex.split(line, comments=False, posix=True)
            if len(fields) >= len(headers):
                rows.append(dict(zip(headers, fields)))
            index += 1
        if rows:
            return rows
    return rows


def first_model(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    models = [row.get("_atom_site.pdbx_PDB_model_num", "1") for row in rows]
    model = "1" if "1" in models else models[0]
    return [row for row in rows if row.get("_atom_site.pdbx_PDB_model_num", "1") == model]


def sequence_number(value: str) -> tuple[str, str]:
    match = re.match(r"^(-?\d+)([A-Za-z]?)$", value)
    if not match:
        return "0", " "
    return match.group(1), match.group(2) or " "


def to_pdb_line(serial: int, row: dict[str, str], chain: str) -> str:
    atom_name = row.get("_atom_site.auth_atom_id") or row.get("_atom_site.label_atom_id") or "X"
    residue = row.get("_atom_site.auth_comp_id") or row.get("_atom_site.label_comp_id") or "UNK"
    residue = residue[:3]
    seq, insertion = sequence_number(row.get("_atom_site.auth_seq_id", "0"))
    chain_field = chain[:2].rjust(2)
    x = float(row.get("_atom_site.Cartn_x", "0"))
    y = float(row.get("_atom_site.Cartn_y", "0"))
    z = float(row.get("_atom_site.Cartn_z", "0"))
    occupancy = float(row.get("_atom_site.occupancy", "1.0"))
    bfactor = float(row.get("_atom_site.B_iso_or_equiv", "0.0"))
    element = (row.get("_atom_site.type_symbol") or "").upper()[:2].rjust(2)
    return (
        f"ATOM  {serial:5d} {atom_name[:4].rjust(4)} {residue:>3}{chain_field}"
        f"{seq:>4}{insertion:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{bfactor:6.2f}          {element}\n"
    )


def extract_chain(cif_path: Path, pdb_id: str, chain: str, output_dir: Path) -> tuple[str, str]:
    rows = first_model(atom_rows(cif_path))
    atom_rows_for_chain = [
        row for row in rows
        if row.get("_atom_site.group_PDB") == "ATOM"
        and row.get("_atom_site.auth_asym_id", "").upper() == chain.upper()
    ]
    if not atom_rows_for_chain:
        atom_rows_for_chain = [
            row for row in rows
            if row.get("_atom_site.group_PDB") == "ATOM"
            and row.get("_atom_site.label_asym_id", "").upper() == chain.upper()
        ]
    if not atom_rows_for_chain:
        return f"{pdb_id}{chain}", "missing_chain"
    output = output_dir / f"{pdb_id}{chain}.pdb"
    if output.exists() and output.stat().st_size > 100:
        return output.stem, "exists"
    content = "".join(to_pdb_line(index, row, chain) for index, row in enumerate(atom_rows_for_chain, 1))
    output.write_text(content + "TER\n", encoding="ascii", errors="ignore")
    return output.stem, f"recovered:{len(atom_rows_for_chain)}_atoms"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--missing-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()

    targets = read_missing(args.missing_list)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)

    pdb_ids = sorted(targets)
    messages: list[str] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [executor.submit(fetch_cif, pdb_id, args.cache_dir, args.retries) for pdb_id in pdb_ids]
        for future in as_completed(futures):
            pdb_id, cif_path, status = future.result()
            messages.append(f"fetch pdb={pdb_id} status={status}")
            if cif_path is None:
                for chain in sorted(targets[pdb_id]):
                    messages.append(f"unrecovered key={pdb_id}{chain} reason={status}")
                continue
            for chain in sorted(targets[pdb_id]):
                key, chain_status = extract_chain(cif_path, pdb_id, chain, args.output_dir)
                messages.append(f"chain key={key} status={chain_status}")

    args.log.write_text("\n".join(sorted(messages)) + "\n", encoding="utf-8")
    recovered = sum("status=recovered:" in message for message in messages)
    exists = sum("status=exists" in message for message in messages)
    unrecovered = sum(message.startswith("unrecovered ") or "status=failed:" in message for message in messages)
    print(f"pdb_ids={len(pdb_ids)} recovered={recovered} exists={exists} unrecovered={unrecovered}")
    for message in sorted(messages):
        if message.startswith("unrecovered ") or "status=failed:" in message or "status=missing_chain" in message:
            print(message)


if __name__ == "__main__":
    main()
