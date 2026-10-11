#!/bin/bash
#SBATCH --job-name=bg4_posneg
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/bg4_posneg_%A_%a.out
#SBATCH --error=logs/bg4_posneg_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
ENV="${GENOMEPLOT_ENV:-genomeplot}"
OUT="$RESULTS/bg4_positive_negative_signal"

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

ENRICH="$RESULTS/g4_enrichment/$LABEL"
SIGNAL="$RESULTS/g4_signal_redraw/$LABEL"

python "$SCRIPT_DIR/compare_bg4_positive_negative_signal.py" \
    --target "$TARGET" \
    --with-bg4 "$ENRICH/${LABEL}_with_BG4.bed" \
    --without-bg4 "$ENRICH/${LABEL}_without_BG4.bed" \
    --bigwig \
        "$SIGNAL/${TARGET}.Rep1_vs_Ctrl1.log2ratio.bw" \
        "$SIGNAL/${TARGET}.Rep2_vs_Ctrl2.log2ratio.bw" \
        "$SIGNAL/${TARGET}.Rep3_vs_Ctrl3.log2ratio.bw" \
    --output-dir "$OUT"

ls -lh \
    "$OUT/${TARGET}.BG4_positive_negative.signal.svg" \
    "$OUT/${TARGET}.BG4_positive_negative.per_peak.tsv" \
    "$OUT/${TARGET}.BG4_positive_negative.summary.tsv"
