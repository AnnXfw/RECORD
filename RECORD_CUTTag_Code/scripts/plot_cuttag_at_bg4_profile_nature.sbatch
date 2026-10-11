#!/bin/bash
#SBATCH --job-name=cuttag_bg4_profile
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/cuttag_bg4_profile_%A_%a.out
#SBATCH --error=logs/cuttag_bg4_profile_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
ENV="${GENOMEPLOT_ENV:-genomeplot}"
OUT="$RESULTS/cuttag_at_bg4_profile_nature"

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

MATRIX="$RESULTS/g4_signal_redraw/$LABEL/${TARGET}.BG4.raw_signal.matrix.gz"

python "$SCRIPT_DIR/plot_cuttag_at_bg4_profile_nature.py" \
    --matrix "$MATRIX" \
    --target "$TARGET" \
    --output "$OUT/${TARGET}.CUTTag_at_BG4.svg"

ls -lh "$OUT/${TARGET}.CUTTag_at_BG4.svg"
