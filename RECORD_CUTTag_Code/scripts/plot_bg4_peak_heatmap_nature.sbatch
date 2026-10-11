#!/bin/bash
#SBATCH --job-name=bg4_heatmap
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/bg4_heatmap_%A_%a.out
#SBATCH --error=logs/bg4_heatmap_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
ENV="${GENOMEPLOT_ENV:-genomeplot}"
OUT="$RESULTS/bg4_peak_heatmap_nature"

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

MATRIX="$RESULTS/g4_signal_redraw/$LABEL/${TARGET}.BG4.log2ratio.matrix.gz"

python "$SCRIPT_DIR/plot_bg4_peak_heatmap_nature.py" \
    --matrix "$MATRIX" \
    --target "$TARGET" \
    --output-dir "$OUT"

ls -lh \
    "$OUT/${TARGET}.BG4_positive_heatmap.svg" \
    "$OUT/${TARGET}.BG4_positive_heatmap.sorted_peaks.tsv"
