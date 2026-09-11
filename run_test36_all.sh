#!/usr/bin/env bash
# One-command AA result -> test36 CG training/generation/MD/analysis.
set -euo pipefail
root_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
py="${TEST36_PYTHON:-$root_dir/../.pixi/envs/default/bin/python}"
aa_run="${AA_RUN:-$root_dir/test36-aa/long-run/run}"
dataset="${DATASET_DIR:-$root_dir/test36-aa/long-run/dataset}"
train="${TRAIN_DIR:-$root_dir/test36-aa/long-run/train}"
generated="${GENERATED_DIR:-$root_dir/test36-aa/long-run/generated}"
cgmd="${CGMD_DIR:-$root_dir/test36-aa/long-run/cgmd}"
analysis="${ANALYSIS_DIR:-$root_dir/test36-aa/long-run/cgmd-analysis}"
trajectory="$aa_run/production.extxyz"
mapping="$aa_run/mapping.json"
checkpoint="$train/checkpoint.pt"

[[ -f "$trajectory" && -f "$mapping" ]] || { echo "Missing AA trajectory or mapping: $aa_run" >&2; exit 2; }
[[ -d "${DM2_ROOT:-$root_dir/DM2}/src/graphite" ]] || { echo "Set DM2_ROOT to a DM2 checkout" >&2; exit 2; }

if [[ ! -f "$dataset/metadata.json" ]]; then
  "$py" -u "$root_dir/test36.py" prepare --trajectory "$trajectory" --mapping "$mapping" --output "$dataset"
fi

if [[ ! -f "$checkpoint" ]]; then
  "$py" -u "$root_dir/test36.py" train --dataset "$dataset" --output "$train" \
    --device "${DEVICE:-mps}" --updates "${UPDATES:-30000}" \
    --batch-size "${BATCH_SIZE:-16}" --time-budget-hours "${TRAIN_TIME_HOURS:-5}"
fi

if [[ ! -f "$generated/generation.json" ]]; then
  "$py" -u "$root_dir/test36.py" generate --checkpoint "$checkpoint" --output "$generated" \
    --device "${DEVICE:-mps}" --time-budget-hours "${GENERATE_TIME_HOURS:-2}"
fi

if [[ ! -f "$cgmd/run_metrics.json" ]]; then
  "$py" -u "$root_dir/test36.py" md --checkpoint "$checkpoint" --output "$cgmd" \
    --device "${DEVICE:-mps}" --sigma-ref "${SIGMA_REF:-0.1}" \
    --steps "${CGMD_STEPS:-10000}" --time-budget-hours "${CGMD_TIME_HOURS:-2}"
fi

if [[ ! -f "$analysis/rdf.json" ]]; then
  "$py" -u "$root_dir/test36.py" analyze --reference "$dataset" \
    --samples "$cgmd/md.extxyz" --output "$analysis" --r-max "${R_MAX:-5}"
fi
echo "test36 pipeline complete: $analysis"
