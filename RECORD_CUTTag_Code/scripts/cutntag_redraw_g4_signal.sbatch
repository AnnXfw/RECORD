#!/usr/bin/env bash
# 重新绘制KLF16和HEY1在“与TF差异峰重叠的G4位点”附近的CUT&Tag信号
# 数组任务1：KLF16；数组任务2：HEY1

#SBATCH --job-name=cuttag_redraw_g4
#SBATCH --partition=cpu
#SBATCH --array=1-2%1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/redraw_g4_%A_%a.out
#SBATCH --error=logs/redraw_g4_%A_%a.err
#SBATCH --export=NONE

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"

module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
    target=KLF16
    label=KLF16_strict_2of3
else
    target=HEY1
    label=HEY1_pooled_default_2of3
fi

ENRICH_DIR="$OUT_ROOT/g4_enrichment/$label"
OUT="$OUT_ROOT/g4_signal_redraw/$label"
mkdir -p "$OUT"

# 这里使用的是被当前TF峰覆盖的G4位点，不再使用全基因组所有G4位点
BG4_REGIONS="$ENRICH_DIR/BG4_hit_by_${label}.bed"
POTENTIAL_REGIONS="$ENRICH_DIR/Potential_G4_hit_by_${label}.bed"

[[ -s "$BG4_REGIONS" ]] || { echo "缺少BG4交集区域：$BG4_REGIONS" >&2; exit 2; }
[[ -s "$POTENTIAL_REGIONS" ]] || { echo "缺少Potential G4交集区域：$POTENTIAL_REGIONS" >&2; exit 3; }

# 三个实验重复和三个相应的C对照RPKM bigWig
REP_TRACKS=()
CTRL_TRACKS=()
for rep in 1 2 3; do
    REP_TRACKS+=("$OUT_ROOT/bigWig/${target}_${rep}.bs50.RPKM.bw")
    CTRL_TRACKS+=("$OUT_ROOT/bigWig/${target}_C${rep}.bs50.RPKM.bw")
done

for bw in "${REP_TRACKS[@]}" "${CTRL_TRACKS[@]}"; do
    [[ -s "$bw" ]] || { echo "缺少bigWig：$bw" >&2; exit 4; }
done

# 分别计算三个重复的log2(实验组/对应C对照)轨迹
RATIO_TRACKS=()
for rep in 1 2 3; do
    ratio_bw="$OUT/${target}.Rep${rep}_vs_Ctrl${rep}.log2ratio.bw"
    bigwigCompare \
        -b1 "${REP_TRACKS[$((rep-1))]}" \
        -b2 "${CTRL_TRACKS[$((rep-1))]}" \
        --operation log2 \
        --pseudocount 1 \
        --binSize 25 \
        --numberOfProcessors "$SLURM_CPUS_PER_TASK" \
        -o "$ratio_bw"
    RATIO_TRACKS+=("$ratio_bw")
done

draw_one_set() {
    local set_name=$1
    local regions=$2

    # 图A：六个原始RPKM轨迹，用于透明展示实验组和C对照
    computeMatrix reference-point \
        --referencePoint center \
        -b 1000 -a 1000 --binSize 25 \
        --missingDataAsZero --skipZeros \
        -S "${REP_TRACKS[@]}" "${CTRL_TRACKS[@]}" \
        -R "$regions" \
        -p "$SLURM_CPUS_PER_TASK" \
        -o "$OUT/${target}.${set_name}.raw_signal.matrix.gz"

    plotProfile \
        -m "$OUT/${target}.${set_name}.raw_signal.matrix.gz" \
        -out "$OUT/${target}.${set_name}.raw_signal.profile.svg" \
        --samplesLabel Rep1 Rep2 Rep3 Ctrl1 Ctrl2 Ctrl3 \
        --regionsLabel "${set_name}_overlapping_${target}_peaks" \
        --refPointLabel "${set_name} center" \
        --plotTitle "${target}: raw RPKM signal at ${set_name} sites" \
        --plotType lines --perGroup --dpi 300 \
        --plotWidth 9 --plotHeight 7

    # 图B：三个log2(实验组/对照组)轨迹，这是主要profile图
    computeMatrix reference-point \
        --referencePoint center \
        -b 1000 -a 1000 --binSize 25 \
        --missingDataAsZero --skipZeros \
        -S "${RATIO_TRACKS[@]}" \
        -R "$regions" \
        -p "$SLURM_CPUS_PER_TASK" \
        -o "$OUT/${target}.${set_name}.log2ratio.matrix.gz"

    plotProfile \
        -m "$OUT/${target}.${set_name}.log2ratio.matrix.gz" \
        -out "$OUT/${target}.${set_name}.log2ratio.profile.svg" \
        --samplesLabel Rep1_vs_Ctrl1 Rep2_vs_Ctrl2 Rep3_vs_Ctrl3 \
        --regionsLabel "${set_name}_overlapping_${target}_peaks" \
        --refPointLabel "${set_name} center" \
        --yAxisLabel "log2(Treatment / Control)" \
        --plotTitle "${target}: treatment/control enrichment at ${set_name} sites" \
        --plotType lines --perGroup --dpi 300 \
        --plotWidth 9 --plotHeight 7

    # 图C：逐位点热图，检查富集是否由少量极强区域驱动
    plotHeatmap \
        -m "$OUT/${target}.${set_name}.log2ratio.matrix.gz" \
        -out "$OUT/${target}.${set_name}.log2ratio.heatmap.svg" \
        --samplesLabel Rep1_vs_Ctrl1 Rep2_vs_Ctrl2 Rep3_vs_Ctrl3 \
        --regionsLabel "${set_name}_overlapping_${target}_peaks" \
        --refPointLabel "${set_name} center" \
        --plotTitle "${target}: log2 treatment/control at ${set_name} sites" \
        --sortUsing mean --sortRegions descend \
        --dpi 300 --heatmapHeight 12 --heatmapWidth 8
}

# BG4作为主结果，Potential G4作为补充结果
draw_one_set BG4 "$BG4_REGIONS"
draw_one_set Potential_G4 "$POTENTIAL_REGIONS"

{
    printf 'target\tpeak_set\tregion_set\tregion_count\n'
    printf '%s\t%s\tBG4\t%s\n' "$target" "$label" "$(wc -l < "$BG4_REGIONS")"
    printf '%s\t%s\tPotential_G4\t%s\n' "$target" "$label" "$(wc -l < "$POTENTIAL_REGIONS")"
} > "$OUT/${target}.redraw_region_counts.tsv"

echo "$target 的G4信号重画完成：$OUT"
