#!/usr/bin/env python3
"""Deduplicate pocket manifests by pocket_id without modifying the originals."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from datetime import datetime
from pathlib import Path


def timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def deduplicate(input_path: Path, output_path: Path, log_handle) -> Counter[str]:
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    with input_path.open("r", encoding="utf-8", errors="replace", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows: list[dict[str, str]] = []
        for row in reader:
            counts["rows_in"] += 1
            pocket_id = (row.get("pocket_id") or "").strip()
            if not pocket_id:
                counts["empty_id"] += 1
                continue
            if pocket_id in seen:
                counts["duplicates_removed"] += 1
                log_handle.write(f"{timestamp()} duplicate split={input_path.stem} pocket_id={pocket_id}\n")
                continue
            seen.add(pocket_id)
            rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts["rows_out"] = len(rows)
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=Path, required=True)
    parser.add_argument("--val", type=Path, required=True)
    parser.add_argument("--test", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    specs = [
        ("train", args.train, args.output_dir / "train_pocket_manifest_dedup.tsv"),
        ("val", args.val, args.output_dir / "val_pocket_manifest_dedup.tsv"),
        ("test", args.test, args.output_dir / "test_biolip_pocket_manifest_dedup.tsv"),
    ]
    total = Counter()
    with args.log.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"{timestamp()} event=start\n")
        for name, input_path, output_path in specs:
            counts = deduplicate(input_path, output_path, log_handle)
            for key, value in counts.items():
                total[key] += value
            log_handle.write(
                f"{timestamp()} split={name} input={input_path} output={output_path} "
                f"counts={dict(counts)}\n"
            )
        log_handle.write(f"{timestamp()} event=summary counts={dict(total)}\n")
    print(f"summary={dict(total)}")


if __name__ == "__main__":
    main()
