#!/usr/bin/env bash
#SBATCH --job-name=fly_fst
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=fly_fst.%j.log
set -euo pipefail
# Edit paths and activate a Python environment containing numpy, scipy, matplotlib.
SEQ_DIR=/path/to/all_seq_files
META=/path/to/metadata.csv
OUT_DIR=/path/to/fst_results
SCRIPT=/path/to/fly_fst.py
python "$SCRIPT" \
  --seq-dir "$SEQ_DIR" \
  --metadata "$META" \
  --outdir "$OUT_DIR" \
  --chunk-size 100000 \
  --min-call-fraction 0.8 \
  --windows 5000:2500,10000:5000,50000:10000,100000:25000
