#!/usr/bin/env python3
"""Convert the Receptor.AI literature pocket JSON files to a TSV index."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.repo_dir / "benchmark-dataset"
    structure_root = root / "structures"
    rows = []
    for json_path in sorted((root / "pockets").rglob("*.json")):
        record = json.loads(json_path.read_text(encoding="utf-8"))
        target = record.get("target", "")
        target_dir = structure_root / json_path.parent.name
        structures = sorted(str(path.relative_to(args.repo_dir)) for path in target_dir.glob("*.pdb"))
        for pocket in record.get("pockets", []):
            pocket_id = pocket.get("pocket_id", "")
            description = " ".join((pocket.get("description") or "").split())
            amino_acids = ";".join(pocket.get("amino_acids") or [])
            ligands = ";".join(pocket.get("ligands") or [])
            rows.append(
                {
                    "pocket_id": f"{record.get('paper_id', json_path.stem)}::{pocket_id}",
                    "paper_id": record.get("paper_id", ""),
                    "target": target,
                    "paper_name": record.get("paper_name", ""),
                    "doi": record.get("DOI", ""),
                    "pocket_description": description,
                    "amino_acids": amino_acids,
                    "ligands": ligands,
                    "structure_paths": ";".join(structures),
                    "source": "Receptor.AI_literature",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "pocket_id", "paper_id", "target", "paper_name", "doi",
        "pocket_description", "amino_acids", "ligands", "structure_paths", "source",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"rows={len(rows)}")
    print(f"papers={len({row['paper_id'] for row in rows})}")
    print(f"targets={len({row['target'] for row in rows})}")
    print(f"rows_with_structure={sum(bool(row['structure_paths']) for row in rows)}")


if __name__ == "__main__":
    main()
