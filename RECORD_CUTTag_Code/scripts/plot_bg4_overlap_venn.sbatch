#!/bin/bash
#SBATCH --job-name=bg4_venn
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/bg4_venn_%A_%a.out
#SBATCH --error=logs/bg4_venn_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
BG4_BED="$ROOT/reference/G4.bed"
ENV="${GENOMEPLOT_ENV:-genomeplot}"
OUT="$RESULTS/bg4_overlap_figures"

module purge
module load miniconda3/24.3.0
source activate "$ENV"
mkdir -p "$OUT"

case "${SLURM_ARRAY_TASK_ID}" in
    1)
        TARGET=KLF16
        LABEL=KLF16_strict_2of3
        ;;
    2)
        TARGET=HEY1
        LABEL=HEY1_pooled_default_2of3
        ;;
    *)
        echo "错误的数组任务编号" >&2
        exit 2
        ;;
esac

SUMMARY="$RESULTS/g4_enrichment/$LABEL/$LABEL.g4_enrichment_summary.tsv"

python "$SCRIPT_DIR/plot_bg4_overlap_venn.py" \
    --summary "$SUMMARY" \
    --bg4-bed "$BG4_BED" \
    --output "$OUT/${TARGET}.BG4_overlap.svg"

ls -lh "$OUT/${TARGET}.BG4_overlap.svg"
