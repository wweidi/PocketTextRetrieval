#!/usr/bin/env bash
set -u

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SESSION="pocket_weight_downloads"

if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "tmux session $SESSION is already running"
else
    tmux new-session -d -s "$SESSION" \
        "bash $PROJECT_ROOT/scripts/run_model_weight_downloads_tmux.sh"
    echo "started tmux session $SESSION"
fi

echo "Check: /mnt/HDD0/home/zf25/miniconda3/envs/protein/bin/python $PROJECT_ROOT/scripts/check_model_weights.py"
echo "Logs:  tail -f $PROJECT_ROOT/logs/model_weights_<profsa|unimol|biomedbert>.log"
