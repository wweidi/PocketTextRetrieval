#!/usr/bin/env python3
"""Download the unique RCSB PDB structures referenced by pocket manifests."""

from __future__ import annotations

import argparse
import csv
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def read_pdb_ids(manifests: list[Path]) -> list[str]:
    ids: set[str] = set()
    for manifest in manifests:
        with manifest.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pdb_id = (row.get("pdb_id") or "").strip().lower()
                if pdb_id:
                    ids.add(pdb_id)
    return sorted(ids)


def download_one(pdb_id: str, output_dir: Path, retries: int) -> tuple[str, str]:
    target = output_dir / f"{pdb_id}.pdb"
    if target.exists() and target.stat().st_size > 100:
        return pdb_id, "exists"
    temporary = output_dir / f".{pdb_id}.pdb.part"
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    request = urllib.request.Request(url, headers={"User-Agent": "PocketTextRetrieval/0.1"})
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=90) as response, temporary.open("wb") as handle:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
            if temporary.stat().st_size <= 100:
                raise OSError("downloaded file is unexpectedly small")
            temporary.replace(target)
            return pdb_id, "downloaded"
        except (OSError, urllib.error.URLError, TimeoutError) as exc:
            if attempt >= retries:
                try:
                    temporary.unlink()
                except FileNotFoundError:
                    pass
                return pdb_id, f"failed:{type(exc).__name__}"
            time.sleep(2 ** attempt)
    return pdb_id, "failed:unknown"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--retries", type=int, default=2)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pdb_ids = read_pdb_ids(args.manifest)
    print(f"unique_pdb={len(pdb_ids)}", flush=True)
    counts = {"exists": 0, "downloaded": 0, "failed": 0}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [executor.submit(download_one, pdb_id, args.output_dir, args.retries) for pdb_id in pdb_ids]
        for index, future in enumerate(as_completed(futures), start=1):
            pdb_id, status = future.result()
            category = status.split(":", 1)[0]
            counts[category] = counts.get(category, 0) + 1
            if category == "failed":
                print(f"failed {pdb_id} {status}", flush=True)
            if index % 100 == 0 or index == len(futures):
                print(f"progress={index}/{len(futures)} counts={counts}", flush=True)


if __name__ == "__main__":
    main()
