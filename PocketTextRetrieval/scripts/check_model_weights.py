#!/usr/bin/env python3
"""Show download state and prepared model-weight files."""

from __future__ import annotations

import json
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATE_ROOT = ROOT / "logs" / "model_weights_state"

EXPECTED = {
    "profsa": [
        ROOT / "checkpoints/profsa/profsa_last.ckpt",
        ROOT / "checkpoints/profsa/mol_pre_no_h_220816.pt",
    ],
    "unimol": [ROOT / "checkpoints/unimol/pocket_pre_220816.pt"],
    "biomedbert": [
        ROOT / "checkpoints/text/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext/pytorch_model.bin",
        ROOT / "checkpoints/text/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext/model.safetensors",
    ],
}


def human(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{size} B"


def pid_state(model: str) -> str:
    pid_file = ROOT / "logs" / f"model_weights_{model}.pid"
    if not pid_file.exists():
        return "no pid"
    try:
        pid = int(pid_file.read_text().strip())
        os.kill(pid, 0)
        return f"running (PID {pid})"
    except (ValueError, OSError):
        return "not running"


for model, candidates in EXPECTED.items():
    state_path = STATE_ROOT / f"{model}.json"
    state = {}
    if state_path.exists():
        try:
            state = json.loads(state_path.read_text())
        except json.JSONDecodeError:
            state = {"status": "invalid state file"}
    status = state.get("status", "not started")
    print(f"\n{model}: {status}; {pid_state(model)}")
    downloaded = state.get("bytes_downloaded")
    total = state.get("total_bytes")
    if isinstance(downloaded, int) and isinstance(total, int) and total > 0:
        print(f"  progress: {human(downloaded)} / {human(total)} ({downloaded * 100 / total:.2f}%)")
    elif isinstance(downloaded, int):
        print(f"  progress: {human(downloaded)}")
    if state.get("error"):
        print(f"  error: {state['error']}")
    for path in candidates:
        if path.exists():
            print(f"  [OK] {path.relative_to(ROOT)} ({human(path.stat().st_size)})")
    partial_roots = {
        "profsa": [ROOT / "checkpoints/profsa"],
        "unimol": [ROOT / "checkpoints/unimol"],
        "biomedbert": [ROOT / "checkpoints/text/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext"],
    }
    partials = []
    for partial_root in partial_roots[model]:
        partials.extend(partial_root.glob("*.part"))
    for partial in sorted(partials):
        print(f"  [part] {partial.relative_to(ROOT)} ({human(partial.stat().st_size)})")
    if model == "biomedbert":
        directory = candidates[0].parent
        if directory.exists():
            for path in sorted(directory.iterdir()):
                if path.is_file() and path.name not in {"pytorch_model.bin", "model.safetensors"}:
                    print(f"  [aux] {path.name} ({human(path.stat().st_size)})")

print("\nFull logs: logs/model_weights_<model>.log")
