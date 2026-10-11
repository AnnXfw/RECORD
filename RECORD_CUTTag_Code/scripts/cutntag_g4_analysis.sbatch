#!/usr/bin/env bash
# 三个靶点分别完成Potential G4交集和deepTools图形；实验性BG4为可选分析
#SBATCH --job-name=cuttag_g4
#SBATCH --partition=cpu
#SBATCH --output=logs/g4_%A_%a.out
#SBATCH --error=logs/g4_%A_%a.err
# Potential G4区域较多，为避免内存和I/O拥堵，每次只运行一个靶点
#SBATCH --array=1-3%1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
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

targets=(HEY1 ING4 KLF16)
default_target=${targets[$((SLURM_ARRAY_TASK_ID-1))]}
# 默认读取正式共识峰；诊断性分析可在提交时指定替代靶点、峰集和输出目录
target=${TARGET_OVERRIDE:-$default_target}
CONSENSUS=${CONSENSUS_OVERRIDE:-$OUT_ROOT/consensus/$target/${target}.consensus.2of3.bed}
ANALYSIS=${ANALYSIS_OVERRIDE:-$OUT_ROOT/analysis/$target}
mkdir -p "$ANALYSIS"

[[ -s "$CONSENSUS" ]] || { echo "缺少2-of-3共识peak：$CONSENSUS" >&2; exit 2; }
[[ -s "$POTENTIAL_G4_BED" ]] || { echo "缺少Potential G4文件：$POTENTIAL_G4_BED" >&2; exit 3; }

# 在计算节点临时目录排序参考BED，避免三个数组任务同时覆盖同一文件
TMP_BASE=${SLURM_TMPDIR:-/tmp/${USER}_cuttag_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
mkdir -p "$TMP_BASE"
trap 'rm -rf -- "$TMP_BASE"' EXIT
POTENTIAL_SORTED="$TMP_BASE/Potential_G4.sorted.bed"
LC_ALL=C sort -k1,1 -k2,2n "$POTENTIAL_G4_BED" > "$POTENTIAL_SORTED"
BG4_SORTED=""
if [[ -n "${BG4_BED:-}" && -s "$BG4_BED" ]]; then
  BG4_SORTED="$TMP_BASE/G4.sorted.bed"
  LC_ALL=C sort -k1,1 -k2,2n "$BG4_BED" > "$BG4_SORTED"
else
  echo "提示：未找到实验性BG4 BED，只执行Potential G4主分析" >&2
fi

# 主分析：分别输出与Potential G4重叠和不重叠的2-of-3共识peak
bedtools intersect -a "$CONSENSUS" -b "$POTENTIAL_SORTED" -u > "$ANALYSIS/${target}_with_Potential_G4.bed"
bedtools intersect -a "$CONSENSUS" -b "$POTENTIAL_SORTED" -v > "$ANALYSIS/${target}_without_Potential_G4.bed"
# 反向统计：输出被当前靶点2-of-3 peak覆盖的Potential G4位点
bedtools intersect -a "$POTENTIAL_SORTED" -b "$CONSENSUS" -u > "$ANALYSIS/Potential_G4_with_${target}.bed"
# 如果提供实验性BG4 BED，则额外输出BG4交集结果
if [[ -n "$BG4_SORTED" ]]; then
  bedtools intersect -a "$CONSENSUS" -b "$BG4_SORTED" -u > "$ANALYSIS/${target}_with_BG4.bed"
  bedtools intersect -a "$CONSENSUS" -b "$BG4_SORTED" -v > "$ANALYSIS/${target}_without_BG4.bed"
fi

# 汇总区域数量和重叠比例，便于后续柱形图使用
consensus_n=$(wc -l < "$CONSENSUS")
potential_n=$(wc -l < "$POTENTIAL_SORTED")
potential_overlap_n=$(wc -l < "$ANALYSIS/${target}_with_Potential_G4.bed")
potential_hit_n=$(wc -l < "$ANALYSIS/Potential_G4_with_${target}.bed")
potential_overlap_pct=$(awk -v a="$potential_overlap_n" -v n="$consensus_n" 'BEGIN{if(n==0) print "NA"; else printf "%.4f",100*a/n}')
potential_hit_pct=$(awk -v a="$potential_hit_n" -v n="$potential_n" 'BEGIN{if(n==0) print "NA"; else printf "%.4f",100*a/n}')
{
  printf 'category\tcount\n'
  printf 'consensus_2of3\t%s\n' "$consensus_n"
  printf 'Potential_G4\t%s\n' "$potential_n"
  printf 'with_Potential_G4\t%s\n' "$potential_overlap_n"
  printf 'without_Potential_G4\t%s\n' "$(wc -l < "$ANALYSIS/${target}_without_Potential_G4.bed")"
  printf 'Potential_G4_overlapped_by_target\t%s\n' "$potential_hit_n"
  if [[ -n "$BG4_SORTED" ]]; then
    printf 'BG4\t%s\n' "$(wc -l < "$BG4_SORTED")"
    printf 'with_BG4\t%s\n' "$(wc -l < "$ANALYSIS/${target}_with_BG4.bed")"
    printf 'without_BG4\t%s\n' "$(wc -l < "$ANALYSIS/${target}_without_BG4.bed")"
  fi
} > "$ANALYSIS/${target}.overlap_counts.tsv"
{
  printf 'target\tconsensus_2of3\tpeaks_with_Potential_G4\tpercent_of_peaks\ttotal_Potential_G4\tPotential_G4_hit_by_target\tpercent_of_Potential_G4\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$target" "$consensus_n" "$potential_overlap_n" "$potential_overlap_pct" \
    "$potential_n" "$potential_hit_n" "$potential_hit_pct"
} > "$ANALYSIS/${target}.Potential_G4_summary.tsv"

# 顺序固定为三个实验重复，再接三个同编号No-primary对照
tracks=()
for rep in 1 2 3; do tracks+=("$OUT_ROOT/bigWig/${target}_${rep}.bs50.RPKM.bw"); done
for rep in 1 2 3; do tracks+=("$OUT_ROOT/bigWig/${target}_C${rep}.bs50.RPKM.bw"); done
for track in "${tracks[@]}"; do
  [[ -s "$track" ]] || { echo "缺少bigWig：$track" >&2; exit 5; }
done
labels=(Rep1 Rep2 Rep3 Ctrl1 Ctrl2 Ctrl3)

# 以Potential G4中心上下游400 bp绘制平均曲线
# Potential G4超过一百万行，不绘制逐行热图，以免产生不可解释且过大的图形
computeMatrix reference-point --referencePoint center --missingDataAsZero \
  -b 400 -a 400 -S "${tracks[@]}" -R "$POTENTIAL_SORTED" --skipZeros \
  -o "$ANALYSIS/${target}_in_Potential_G4.mat.gz" -p "$SLURM_CPUS_PER_TASK"
plotProfile -m "$ANALYSIS/${target}_in_Potential_G4.mat.gz" \
  -out "$ANALYSIS/${target}_in_Potential_G4.profile.svg" \
  --refPointLabel "Potential G4 Center" --regionsLabel "Potential G4" \
  --samplesLabel "${labels[@]}" --perGroup --dpi 300 --plotHeight 8 --plotWidth 8
# 可选：以实验性BG4中心上下游400 bp绘图
if [[ -n "$BG4_SORTED" ]]; then
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

echo "$target 的G4下游分析完成"
