#!/usr/bin/env bash
# 9组实验样本与匹配no-primary对照的peak calling
#SBATCH --job-name=cuttag_peaks
#SBATCH --partition=cpu
#SBATCH --output=logs/peaks_%A_%a.out
#SBATCH --error=logs/peaks_%A_%a.err
#SBATCH --array=1-9%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"
module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print; exit}' "$PAIRS_TSV")
IFS=$'\t' read -r pair_id target replicate experiment no_primary <<< "$line"
[[ -n "$pair_id" ]] || { echo "No pair for array index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

TREATMENT="$OUT_ROOT/bam/$experiment/${experiment}.final.bam"
CONTROL="$OUT_ROOT/bam/$no_primary/${no_primary}.final.bam"
[[ -s "$TREATMENT" && -s "$CONTROL" ]] || { echo "Missing paired BAM: $TREATMENT $CONTROL" >&2; exit 3; }

P="$OUT_ROOT/peaks/$pair_id"
mkdir -p "$P"
macs3 callpeak -t "$TREATMENT" -c "$CONTROL" -f BAMPE -g hs \
  -n "$pair_id" --outdir "$P" --qvalue 0.05 --keep-dup all --call-summits

sort -k1,1 -k2,2n "$P/${pair_id}_peaks.narrowPeak" > "$P/${pair_id}.sorted.narrowPeak"
# 如果提供了hg19 blacklist则过滤；没有时沿用排序后的原始peak
if [[ -n "${BLACKLIST:-}" && -s "$BLACKLIST" ]]; then
  bedtools intersect -v -a "$P/${pair_id}.sorted.narrowPeak" -b "$BLACKLIST" \
    > "$P/${pair_id}.analysis.narrowPeak"
else
  cp "$P/${pair_id}.sorted.narrowPeak" "$P/${pair_id}.analysis.narrowPeak"
  echo "提示：未找到hg19 blacklist，本次未进行blacklist过滤" >&2
fi
