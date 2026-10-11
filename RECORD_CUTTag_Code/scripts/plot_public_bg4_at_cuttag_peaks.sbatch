#!/bin/bash
#SBATCH --job-name=public_bg4
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/public_bg4_%A_%a.out
#SBATCH --error=logs/public_bg4_%A_%a.err

set -euo pipefail

# 项目路径
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="$ROOT/results"
BG4_DIR="$ROOT/reference/bg4_k562_gse107690"
OUT_ROOT="$RESULTS/public_bg4_at_cuttag_peaks"
CONDA_ENV="${CUTTAG_ENV:-cuttag}"

# 加载现有CUT&Tag分析环境
module purge
module load miniconda3/24.3.0
source activate "$CONDA_ENV"

# 数组任务1分析KLF16，任务2分析HEY1
case "${SLURM_ARRAY_TASK_ID}" in
    1)
        TARGET=KLF16
        PEAKS="$RESULTS/consensus/KLF16/KLF16.consensus.2of3.bed"
        ;;
    2)
        TARGET=HEY1
        PEAKS="$RESULTS/diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed"
        ;;
    *)
        echo "错误的数组任务编号：${SLURM_ARRAY_TASK_ID}" >&2
        exit 2
        ;;
esac

OUT="$OUT_ROOT/$TARGET"
mkdir -p "$OUT"

# 论文GSE107690中的5个BG4技术重复
BG4_TRACKS=(
    "$BG4_DIR/GSM2876090_B_REP1.SLX-12319.K562_asynch_a_701_504.clean.bw"
    "$BG4_DIR/GSM2876091_B_REP1.SLX-12319.K562_asynch_c_703_504.clean.bw"
    "$BG4_DIR/GSM2876093_B_REP2.SLX-12320.K562_P9_Async_c_701_503.rmdup.clean.bw"
    "$BG4_DIR/GSM2876094_B_REP2.SLX-12320.K562_P9_Async_b_701_502.rmdup.clean.bw"
    "$BG4_DIR/GSM2876095_B_REP2.SLX-12320.K562_P9_Async_a_701_517.rmdup.clean.bw"
)

# GEO文件名中写作IP，但GEO样本元数据标注为两个Input
INPUT_TRACKS=(
    "$BG4_DIR/GSM2876092_B_REP1.SLX-12319.K562_asynch_IP_704_504.clean.bw"
    "$BG4_DIR/GSM2876096_B_REP2.SLX-12320.K562_P9_Async_IP_702_504.rmdup.clean.bw"
)

[[ -s "$PEAKS" ]] || { echo "缺少或为空的峰文件：$PEAKS" >&2; exit 3; }
for f in "${BG4_TRACKS[@]}" "${INPUT_TRACKS[@]}"; do
    [[ -s "$f" ]] || { echo "缺少BigWig：$f" >&2; exit 4; }
done

# 去除可能存在的表头/空行，并按hg19自然染色体顺序排序
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print}' \
    "$PEAKS" | sort -k1,1V -k2,2n > "$OUT/${TARGET}.peaks.clean.bed"

PEAK_COUNT=$(wc -l < "$OUT/${TARGET}.peaks.clean.bed")
[[ "$PEAK_COUNT" -gt 0 ]] || { echo "清理后没有可分析峰" >&2; exit 5; }
echo "$TARGET用于绘图的峰数量：$PEAK_COUNT"

# 在蛋白峰中心上下游2 kb提取公开BG4与Input信号
computeMatrix reference-point \
    --referencePoint center \
    -R "$OUT/${TARGET}.peaks.clean.bed" \
    -S "${BG4_TRACKS[@]}" "${INPUT_TRACKS[@]}" \
    --samplesLabel BG4_R1a BG4_R1c BG4_R2c BG4_R2b BG4_R2a Input_R1 Input_R2 \
    -b 2000 -a 2000 \
    --binSize 25 \
    --missingDataAsZero \
    --skipZeros \
    --numberOfProcessors "$SLURM_CPUS_PER_TASK" \
    -o "$OUT/${TARGET}.public_BG4_at_peaks.matrix.gz" \
    --outFileSortedRegions "$OUT/${TARGET}.matrix_regions.bed"

# 每条线代表一个公开BG4技术重复或Input
plotProfile \
    -m "$OUT/${TARGET}.public_BG4_at_peaks.matrix.gz" \
    -out "$OUT/${TARGET}.public_BG4_at_peaks.profile.svg" \
    --plotFileFormat svg \
    --plotTitle "Published K562 BG4 signal at ${TARGET} CUT&Tag peaks" \
    --refPointLabel "${TARGET} peak center" \
    --regionsLabel "${TARGET} peaks" \
    --plotHeight 7 \
    --plotWidth 11 \
    --colors '#08306B' '#2171B5' '#4292C6' '#6BAED6' '#9ECAE1' '#D7301F' '#FC8D59' \
    --legendLocation best

# 热图直观看每一个蛋白峰处的BG4信号
plotHeatmap \
    -m "$OUT/${TARGET}.public_BG4_at_peaks.matrix.gz" \
    -out "$OUT/${TARGET}.public_BG4_at_peaks.heatmap.svg" \
    --plotFileFormat svg \
    --plotTitle "Published K562 BG4 signal at ${TARGET} CUT&Tag peaks" \
    --refPointLabel "${TARGET} peak center" \
    --regionsLabel "${TARGET} peaks" \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 14 \
    --heatmapWidth 5

printf "target\tpeak_count\tgenome\tBG4_source\n%s\t%s\thg19\tGSE107690_K562\n" \
    "$TARGET" "$PEAK_COUNT" > "$OUT/${TARGET}.public_BG4_plot_summary.tsv"

echo "分析完成：$OUT"
