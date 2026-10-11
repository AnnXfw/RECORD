#!/bin/bash
#SBATCH --job-name=bg4_svg
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/bg4_svg_%A_%a.out
#SBATCH --error=logs/bg4_svg_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="$ROOT/results/public_bg4_at_cuttag_peaks"
CONDA_ENV="${CUTTAG_ENV:-cuttag}"

module purge
module load miniconda3/24.3.0
source activate "$CONDA_ENV"

case "${SLURM_ARRAY_TASK_ID}" in
    1) TARGET=KLF16 ;;
    2) TARGET=HEY1 ;;
    *) echo "错误的数组任务编号" >&2; exit 2 ;;
esac

OUT="$OUT_ROOT/$TARGET"
MATRIX="$OUT/${TARGET}.public_BG4_at_peaks.matrix.gz"
[[ -s "$MATRIX" ]] || { echo "缺少矩阵：$MATRIX" >&2; exit 3; }

# 曲线图：SVG中曲线、文字、坐标轴和图例均为矢量对象
plotProfile \
    -m "$MATRIX" \
    -out "$OUT/${TARGET}.public_BG4_at_peaks.profile.svg" \
    --plotFileFormat svg \
    --plotTitle "Published K562 BG4 signal at ${TARGET} CUT&Tag peaks" \
    --refPointLabel "${TARGET} peak center" \
    --regionsLabel "${TARGET} peaks" \
    --plotHeight 7 \
    --plotWidth 11 \
    --colors '#08306B' '#2171B5' '#4292C6' '#6BAED6' '#9ECAE1' '#D7301F' '#FC8D59' \
    --legendLocation best

# 热图：输出SVG容器，标签和坐标轴便于编辑
plotHeatmap \
    -m "$MATRIX" \
    -out "$OUT/${TARGET}.public_BG4_at_peaks.heatmap.svg" \
    --plotFileFormat svg \
    --plotTitle "Published K562 BG4 signal at ${TARGET} CUT&Tag peaks" \
    --refPointLabel "${TARGET} peak center" \
    --regionsLabel "${TARGET} peaks" \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 14 \
    --heatmapWidth 5

ls -lh "$OUT/${TARGET}.public_BG4_at_peaks.profile.svg" \
       "$OUT/${TARGET}.public_BG4_at_peaks.heatmap.svg"
