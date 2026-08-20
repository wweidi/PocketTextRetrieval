#!/usr/bin/env python3
"""Download and prepare the three pretrained encoders used by this project.

The script is intentionally resumable.  Every model has its own ``.part`` file
and JSON status file, so an interrupted process can be started again safely.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import tarfile
import time
from pathlib import Path
from typing import BinaryIO, Iterable

import requests


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT_ROOT = PROJECT_ROOT / "checkpoints"
LOG_ROOT = PROJECT_ROOT / "logs"
STATE_ROOT = LOG_ROOT / "model_weights_state"

PROFSA_FILE_ID = "1lFBe4ak7QXS4LS-qAemvWJatT9AL8huf"
UNIMOL_URL = (
    "https://github.com/deepmodeling/Uni-Mol/releases/download/"
    "v0.1/pocket_pre_220816.pt"
)
BIOMEDBERT_REPO = "microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext"


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def write_state(model: str, **values: object) -> None:
    STATE_ROOT.mkdir(parents=True, exist_ok=True)
    path = STATE_ROOT / f"{model}.json"
    state: dict[str, object] = {"model": model, "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    if path.exists():
        try:
            state.update(json.loads(path.read_text()))
        except (OSError, json.JSONDecodeError):
            pass
    state.update(values)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(state, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def response_total(response: requests.Response, initial_size: int = 0) -> int | None:
    content_range = response.headers.get("Content-Range", "")
    match = re.search(r"/([0-9]+)$", content_range)
    if match:
        return int(match.group(1))
    length = response.headers.get("Content-Length")
    if length and length.isdigit():
        return initial_size + int(length)
    return None


def copy_response(
    response: requests.Response,
    destination_part: Path,
    model: str,
    initial_chunk: bytes = b"",
    append: bool = False,
    total: int | None = None,
) -> None:
    mode = "ab" if append else "wb"
    written = destination_part.stat().st_size if append and destination_part.exists() else 0
    last_report = time.monotonic()
    last_bytes = written
    with destination_part.open(mode) as output:
        if initial_chunk:
            output.write(initial_chunk)
            written += len(initial_chunk)
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if not chunk:
                continue
            output.write(chunk)
            written += len(chunk)
            now = time.monotonic()
            if now - last_report >= 5:
                speed = (written - last_bytes) / max(now - last_report, 1e-6) / 1024 / 1024
                if total:
                    percent = written * 100 / total
                    log(f"{model}: {written / 1024**3:.2f}/{total / 1024**3:.2f} GiB ({percent:.1f}%), {speed:.1f} MiB/s")
                else:
                    log(f"{model}: {written / 1024**3:.2f} GiB downloaded, {speed:.1f} MiB/s")
                last_report, last_bytes = now, written
                write_state(model, status="downloading", bytes_downloaded=written, total_bytes=total)
    log(f"{model}: download stream finished at {written / 1024**3:.2f} GiB")
    write_state(model, status="downloaded", bytes_downloaded=written, total_bytes=total)


def download_http(url: str, destination: Path, model: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix(destination.suffix + ".part")
    existing = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    log(f"{model}: requesting {url}")
    with requests.get(url, headers=headers, stream=True, timeout=(30, 120), allow_redirects=True) as response:
        response.raise_for_status()
        if existing and response.status_code == 206:
            total = response_total(response, existing)
            log(f"{model}: resuming at {existing / 1024**3:.2f} GiB")
            copy_response(response, part, model, append=True, total=total)
        else:
            if existing:
                log(f"{model}: server did not honor Range; restarting partial file")
            total = response_total(response)
            copy_response(response, part, model, total=total)
    part.replace(destination)


def download_google_drive(destination: Path, model: str) -> None:
    """Download a large Google Drive file, handling the confirmation page."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix(destination.suffix + ".part")
    existing = part.stat().st_size if part.exists() else 0
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0 PocketTextRetrieval model downloader"})
    url = "https://drive.usercontent.google.com/download"
    params = {"id": PROFSA_FILE_ID, "export": "download", "confirm": "t"}
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    log(f"{model}: requesting Google Drive archive")
    response = session.get(url, params=params, headers=headers, stream=True, timeout=(30, 120), allow_redirects=True)
    response.raise_for_status()

    if existing and response.status_code == 206:
        total = response_total(response, existing)
        copy_response(response, part, model, append=True, total=total)
    else:
        if existing:
            log(f"{model}: Google Drive did not honor Range; restarting partial archive")
        first = response.raw.read(128 * 1024)
        if not first.startswith(b"\x1f\x8b"):
            response.close()
            page = first.decode("utf-8", errors="ignore")
            token_match = re.search(r"confirm=([0-9A-Za-z_-]+)", page)
            if not token_match:
                token_match = re.search(r"name=\"confirm\" value=\"([^\"]+)", page)
            if not token_match:
                raise RuntimeError("Google Drive returned an HTML confirmation/error page")
            params["confirm"] = token_match.group(1)
            response = session.get(url, params=params, stream=True, timeout=(30, 120), allow_redirects=True)
            response.raise_for_status()
            first = response.raw.read(128 * 1024)
            if not first.startswith(b"\x1f\x8b"):
                raise RuntimeError("Google Drive confirmation did not return a gzip archive")
        total = response_total(response)
        copy_response(response, part, model, initial_chunk=first, total=total)
    response.close()
    part.replace(destination)


def find_member(tar: tarfile.TarFile, suffixes: Iterable[str]) -> tarfile.TarInfo:
    members = tar.getmembers()
    for suffix in suffixes:
        for member in members:
            if member.isfile() and member.name.endswith(suffix):
                return member
    raise FileNotFoundError(f"No archive member matched: {list(suffixes)}")


def extract_profsa_weights(archive: Path, output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = output_dir / "profsa_last.ckpt"
    mol_pretrain = output_dir / "mol_pre_no_h_220816.pt"
    with tarfile.open(archive, "r:gz") as tar:
        checkpoint_member = find_member(
            tar,
            ("data/log/train/profsa/profsa_release/checkpoints/last.ckpt", "profsa_release/checkpoints/last.ckpt"),
        )
        mol_member = find_member(tar, ("data/pretrain/mol_pre_no_h_220816.pt", "mol_pre_no_h_220816.pt"))
        for member, destination in ((checkpoint_member, checkpoint), (mol_member, mol_pretrain)):
            source = tar.extractfile(member)
            if source is None:
                raise RuntimeError(f"Cannot read archive member {member.name}")
            with source, destination.open("wb") as target:
                shutil.copyfileobj(source, target, length=1024 * 1024)
            log(f"profsa: extracted {member.name} -> {destination} ({destination.stat().st_size / 1024**2:.1f} MiB)")
    return checkpoint, mol_pretrain


def run_profsa() -> None:
    model = "profsa"
    output_dir = CHECKPOINT_ROOT / "profsa"
    archive = output_dir / "profsa.tar.gz"
    checkpoint = output_dir / "profsa_last.ckpt"
    mol_pretrain = output_dir / "mol_pre_no_h_220816.pt"
    if checkpoint.exists() and mol_pretrain.exists():
        log(f"{model}: already prepared; skipping")
        write_state(model, status="complete", files=[str(checkpoint), str(mol_pretrain)])
        return
    write_state(model, status="downloading", destination=str(archive))
    download_google_drive(archive, model)
    write_state(model, status="extracting", archive=str(archive))
    extract_profsa_weights(archive, output_dir)
    write_state(model, status="complete", files=[str(checkpoint), str(mol_pretrain)])
    log(f"{model}: complete")


def run_unimol() -> None:
    model = "unimol"
    destination = CHECKPOINT_ROOT / "unimol" / "pocket_pre_220816.pt"
    if destination.exists() and destination.stat().st_size > 1_000_000:
        log(f"{model}: already prepared; skipping")
        write_state(model, status="complete", files=[str(destination)])
        return
    write_state(model, status="downloading", destination=str(destination), url=UNIMOL_URL)
    download_http(UNIMOL_URL, destination, model)
    if destination.stat().st_size <= 1_000_000:
        raise RuntimeError("Uni-Mol file is unexpectedly small; download may be an error page")
    write_state(model, status="complete", files=[str(destination)])
    log(f"{model}: complete")


def run_biomedbert() -> None:
    model = "biomedbert"
    destination = CHECKPOINT_ROOT / "text" / "BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext"
    weight_candidates = [destination / "pytorch_model.bin", destination / "model.safetensors"]
    if any(path.exists() and path.stat().st_size > 1_000_000 for path in weight_candidates):
        log(f"{model}: already prepared; skipping")
        write_state(model, status="complete", files=[str(path) for path in destination.iterdir()])
        return
    destination.mkdir(parents=True, exist_ok=True)
    write_state(model, status="downloading", destination=str(destination), repo_id=BIOMEDBERT_REPO)
    log(f"{model}: downloading required Hugging Face files from {BIOMEDBERT_REPO}")
    from huggingface_hub import snapshot_download

    snapshot_download(
        repo_id=BIOMEDBERT_REPO,
        local_dir=str(destination),
        allow_patterns=[
            "config.json",
            "tokenizer.json",
            "tokenizer_config.json",
            "special_tokens_map.json",
            "vocab.txt",
            "pytorch_model.bin",
            "model.safetensors",
        ],
    )
    if not any(path.exists() and path.stat().st_size > 1_000_000 for path in weight_candidates):
        raise RuntimeError("BiomedBERT weight file was not found after Hugging Face download")
    write_state(model, status="complete", files=[str(path) for path in destination.iterdir()])
    log(f"{model}: complete")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=("profsa", "unimol", "biomedbert"), required=True)
    args = parser.parse_args()
    try:
        if args.model == "profsa":
            run_profsa()
        elif args.model == "unimol":
            run_unimol()
        else:
            run_biomedbert()
    except Exception as exc:
        log(f"{args.model}: FAILED: {type(exc).__name__}: {exc}")
        write_state(args.model, status="failed", error=f"{type(exc).__name__}: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
