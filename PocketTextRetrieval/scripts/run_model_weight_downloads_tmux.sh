#!/usr/bin/env bash
set -u

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="/mnt/HDD0/home/zf25/miniconda3/envs/protein/bin/python"

mkdir -p "$PROJECT_ROOT/logs/model_weights_state"

for model in profsa unimol biomedbert; do
    log_file="$PROJECT_ROOT/logs/model_weights_${model}.log"
    pid_file="$PROJECT_ROOT/logs/model_weights_${model}.pid"
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/download_model_weights.py" \
        --model "$model" >> "$log_file" 2>&1 &
    echo $! > "$pid_file"
    echo "started $model with PID $(<"$pid_file")"
done

wait
