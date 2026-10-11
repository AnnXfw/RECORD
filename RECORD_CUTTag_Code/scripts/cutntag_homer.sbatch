#!/usr/bin/env bash
# 对三个靶点的2-of-3共识peak分别进行HOMER motif分析
#SBATCH --job-name=cuttag_homer
#SBATCH --partition=cpu
#SBATCH --output=logs/homer_%A_%a.out
#SBATCH --error=logs/homer_%A_%a.err
#SBATCH --array=1-3%1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --export=NONE

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"
module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

command -v findMotifsGenome.pl >/dev/null 2>&1 || {
  echo "当前cuttag环境尚未安装HOMER，安装并配置hg19后再提交此作业" >&2
  exit 2
}

targets=(HEY1 ING4 KLF16)
target=${targets[$((SLURM_ARRAY_TASK_ID-1))]}
PEAKS="$OUT_ROOT/consensus/$target/${target}.consensus.2of3.bed"
MOTIF_OUT="$OUT_ROOT/motif/${target}_hg19"
[[ -s "$PEAKS" ]] || { echo "缺少2-of-3共识peak：$PEAKS" >&2; exit 3; }
mkdir -p "$OUT_ROOT/motif"

findMotifsGenome.pl "$PEAKS" hg19 "$MOTIF_OUT" -size 200 -p "$SLURM_CPUS_PER_TASK"
