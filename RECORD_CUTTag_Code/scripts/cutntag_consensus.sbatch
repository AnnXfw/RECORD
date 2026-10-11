#!/usr/bin/env bash
# 每个靶点根据3个生物学重复生成1-of-3、2-of-3和3-of-3 consensus peaks
#SBATCH --job-name=cuttag_consensus
#SBATCH --partition=cpu
#SBATCH --output=logs/consensus_%A_%a.out
#SBATCH --error=logs/consensus_%A_%a.err
#SBATCH --array=1-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --export=NONE
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"
module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

targets=(HEY1 ING4 KLF16)
target=${targets[$((SLURM_ARRAY_TASK_ID-1))]}
C="$OUT_ROOT/consensus/$target"
mkdir -p "$C"

peaks=()
for rep in 1 2 3; do
  peak="$OUT_ROOT/peaks/${target}_rep${rep}/${target}_rep${rep}.analysis.narrowPeak"
  [[ -s "$peak" ]] || { echo "Missing peak file: $peak" >&2; exit 3; }
  peaks+=("$peak")
done

bedtools multiinter -i "${peaks[@]}" > "$C/${target}.multiinter.bed"
for support in 1 2 3; do
  awk -v n="$support" 'BEGIN{OFS="\t"} $4>=n {print $1,$2,$3}' "$C/${target}.multiinter.bed" \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$C/${target}.consensus.${support}of3.bed"
done
