#!/usr/bin/env bash
# 复用已经完成的Potential G4矩阵，补齐平均曲线和可选BG4图形
# 数组1为HEY1 pooled-default，数组2为KLF16正式结果
#SBATCH --job-name=cuttag_g4_finish
#SBATCH --partition=cpu
#SBATCH --output=logs/g4_finish_%A_%a.out
#SBATCH --error=logs/g4_finish_%A_%a.err
#SBATCH --array=1-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
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

if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
    target=HEY1
    ANALYSIS="$OUT_ROOT/diagnostic_g4/HEY1_pooled_default"
else
    target=KLF16
    ANALYSIS="$OUT_ROOT/analysis/KLF16"
fi

mkdir -p "$ANALYSIS"
POTENTIAL_MATRIX="$ANALYSIS/${target}_in_Potential_G4.mat.gz"
[[ -s "$POTENTIAL_MATRIX" ]] || {
    echo "缺少已完成的Potential G4矩阵：$POTENTIAL_MATRIX" >&2
    exit 2
}

labels=(Rep1 Rep2 Rep3 Ctrl1 Ctrl2 Ctrl3)

# 复用矩阵生成Potential G4平均曲线，不重新运行computeMatrix
plotProfile -m "$POTENTIAL_MATRIX" \
  -out "$ANALYSIS/${target}_in_Potential_G4.profile.svg" \
  --refPointLabel "Potential G4 Center" --regionsLabel "Potential G4" \
  --samplesLabel "${labels[@]}" --perGroup --dpi 300 --plotHeight 8 --plotWidth 8

# BG4为较小的实验性高可信区域，可生成平均曲线和热图
if [[ -n "${BG4_BED:-}" && -s "$BG4_BED" ]]; then
    TMP_BASE=${SLURM_TMPDIR:-/tmp/${USER}_g4_finish_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
    mkdir -p "$TMP_BASE"
    trap 'rm -rf -- "$TMP_BASE"' EXIT
    BG4_SORTED="$TMP_BASE/G4.sorted.bed"
    LC_ALL=C sort -k1,1 -k2,2n "$BG4_BED" > "$BG4_SORTED"

    tracks=()
    for rep in 1 2 3; do tracks+=("$OUT_ROOT/bigWig/${target}_${rep}.bs50.RPKM.bw"); done
    for rep in 1 2 3; do tracks+=("$OUT_ROOT/bigWig/${target}_C${rep}.bs50.RPKM.bw"); done
    for track in "${tracks[@]}"; do
        [[ -s "$track" ]] || { echo "缺少bigWig：$track" >&2; exit 3; }
    done

    computeMatrix reference-point --referencePoint center --missingDataAsZero \
      -b 400 -a 400 -S "${tracks[@]}" -R "$BG4_SORTED" --skipZeros \
      -o "$ANALYSIS/${target}_in_BG4.mat.gz" -p "$SLURM_CPUS_PER_TASK"
    plotProfile -m "$ANALYSIS/${target}_in_BG4.mat.gz" \
      -out "$ANALYSIS/${target}_in_BG4.profile.svg" \
      --refPointLabel "G4 Center" --regionsLabel BG4 --samplesLabel "${labels[@]}" \
      --perGroup --dpi 300 --plotHeight 8 --plotWidth 8
    plotHeatmap -m "$ANALYSIS/${target}_in_BG4.mat.gz" \
      -out "$ANALYSIS/${target}_in_BG4.heatmap.svg" \
      --refPointLabel "G4 Center" --regionsLabel BG4 --samplesLabel "${labels[@]}" \
      --dpi 300 --heatmapHeight 12 --heatmapWidth 10
fi

echo "$target 的G4绘图续跑完成"
