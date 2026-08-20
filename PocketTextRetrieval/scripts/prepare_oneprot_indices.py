#!/usr/bin/env python3
"""Normalize the downloaded OneProt modality files into explicit indices."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_one_column(path: Path) -> list[str]:
    values = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            value = line.strip().split(",", 1)[0]
            if value and value.lower() not in {"id", "protein_id"}:
                values.append(value)
    return values


def read_text(path: Path) -> dict[str, str]:
    result = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) >= 2 and row[0].strip():
                result[row[0].strip()] = row[1]
    return result


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    test_text = read_text(args.raw_dir / "test_text.csv")
    val_text = read_text(args.raw_dir / "val_text.csv")
    test_pocket = read_one_column(args.raw_dir / "test_pocket.csv")
    val_pocket = read_one_column(args.raw_dir / "val_pocket.csv")

    write_rows(
        args.output_dir / "oneprot_test_text.tsv",
        ["protein_id", "protein_text", "source"],
        [{"protein_id": key, "protein_text": value, "source": "OneProt"} for key, value in sorted(test_text.items())],
    )
    write_rows(
        args.output_dir / "oneprot_val_text.tsv",
        ["protein_id", "protein_text", "source"],
        [{"protein_id": key, "protein_text": value, "source": "OneProt"} for key, value in sorted(val_text.items())],
    )
    write_rows(
        args.output_dir / "oneprot_test_pocket.tsv",
        ["protein_id", "pocket_available", "source"],
        [{"protein_id": key, "pocket_available": "1", "source": "OneProt"} for key in sorted(set(test_pocket))],
    )
    write_rows(
        args.output_dir / "oneprot_val_pocket.tsv",
        ["protein_id", "pocket_available", "source"],
        [{"protein_id": key, "pocket_available": "1", "source": "OneProt"} for key in sorted(set(val_pocket))],
    )

    paired_fields = ["protein_id", "protein_text", "pocket_available", "source"]
    paired_rows = [
        {
            "protein_id": key,
            "protein_text": test_text[key],
            "pocket_available": "1",
            "source": "OneProt_test_all",
        }
        for key in sorted(set(test_text) & set(test_pocket))
    ]
    write_rows(args.output_dir / "oneprot_test_paired.tsv", paired_fields, paired_rows)

    # Preserve the official combined alignment rows as a separate file. This
    # file is intentionally not renamed as a pocket-text dataset: its text is
    # protein-level UniProt text.
    with (args.raw_dir / "test_all.csv").open("r", encoding="utf-8", errors="replace", newline="") as source, (
        args.output_dir / "oneprot_test_all.tsv"
    ).open("w", encoding="utf-8", newline="") as target:
        reader = csv.DictReader(source)
        if reader.fieldnames:
            fields = [field for field in ("id", "func_text", "structure", "sequence", "pocket") if field in reader.fieldnames]
            writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in reader:
                writer.writerow({field: row.get(field, "") for field in fields})

    print(f"test_text={len(test_text)}")
    print(f"test_pocket={len(set(test_pocket))}")
    print(f"val_text={len(val_text)}")
    print(f"val_pocket={len(set(val_pocket))}")
    print(f"paired_test_intersection={len(paired_rows)}")


if __name__ == "__main__":
    main()
