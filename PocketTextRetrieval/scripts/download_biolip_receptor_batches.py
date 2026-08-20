#!/usr/bin/env python3
"""Collect BioLiP chain-level receptor structures from official batch archives."""

from __future__ import annotations

import argparse
import csv
import io
import tarfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_targets(manifests: list[Path]) -> dict[str, set[str]]:
    targets: dict[str, set[str]] = {}
    for manifest in manifests:
        with manifest.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pdb_id = (row.get("pdb_id") or "").strip().lower()
                chain = (row.get("receptor_chain") or "").strip()
                if pdb_id and chain:
                    targets.setdefault(pdb_id[1:3], set()).add(pdb_id + chain)
    return targets


def fetch(url: str, retries: int = 3) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "PocketTextRetrieval/0.1"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except (OSError, urllib.error.URLError):
            if attempt + 1 == retries:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")


def process_prefix(prefix: str, wanted: set[str], output_dir: Path) -> tuple[str, dict[str, int], str]:
    local = {"archives": 0, "extracted": 0, "missing": 0, "errors": 0}
    url = f"https://seq2fun.dcmb.med.umich.edu/BioLiP/weekly/receptor_{prefix}_nr.tar.bz2"
    try:
        payload = fetch(url)
        local["archives"] = 1
        extracted = set()
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:bz2") as archive:
            for member in archive.getmembers():
                name = Path(member.name).name
                if not name.endswith(".pdb"):
                    continue
                key = name[:-4]
                if key not in wanted:
                    continue
                output = output_dir / name
                if output.exists() and output.stat().st_size > 100:
                    extracted.add(key)
                    continue
                source = archive.extractfile(member)
                if source is None:
                    continue
                output.write_bytes(source.read())
                extracted.add(key)
                local["extracted"] += 1
        local["missing"] = len(wanted - extracted)
        message = f"prefix={prefix} wanted={len(wanted)} extracted={len(extracted)} missing={len(wanted - extracted)}"
    except Exception as exc:
        local["errors"] = 1
        message = f"error prefix={prefix} type={type(exc).__name__} detail={exc}"
    return prefix, local, message


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--keep-archives", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    targets = read_targets(args.manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts = {"archives": 0, "extracted": 0, "missing": 0, "errors": 0}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [executor.submit(process_prefix, prefix, targets[prefix], args.output_dir) for prefix in sorted(targets)]
        for future in as_completed(futures):
            _, local, message = future.result()
            for key, value in local.items():
                counts[key] += value
            print(message, flush=True)
    print("summary=" + repr(counts))


if __name__ == "__main__":
    main()
