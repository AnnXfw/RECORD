#!/usr/bin/env bash
# 汇总三个靶点与Potential G4的双向重叠统计
#SBATCH --job-name=cuttag_g4_collect
#SBATCH --partition=cpu
#SBATCH --output=logs/g4_collect_%j.out
#SBATCH --error=logs/g4_collect_%j.err
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
SUMMARY="$OUT_ROOT/analysis/Potential_G4_all_targets_summary.tsv"
MULTI="$OUT_ROOT/analysis/Potential_G4_target_multiinter.bed"
mkdir -p "$OUT_ROOT/analysis"

first=1
hit_beds=()
for target in "${targets[@]}"; do
  table="$OUT_ROOT/analysis/$target/${target}.Potential_G4_summary.tsv"
  hit_bed="$OUT_ROOT/analysis/$target/Potential_G4_with_${target}.bed"
  [[ -s "$table" ]] || { echo "缺少汇总表：$table" >&2; exit 2; }
  [[ -f "$hit_bed" ]] || { echo "缺少Potential G4交集BED：$hit_bed" >&2; exit 3; }
  if [[ "$first" -eq 1 ]]; then
    cat "$table" > "$SUMMARY"
    first=0
  else
    tail -n +2 "$table" >> "$SUMMARY"
  fi
  hit_beds+=("$hit_bed")
done

# 比较哪些Potential G4区域被一个、两个或三个靶点的peak覆盖
bedtools multiinter -header -names "${targets[@]}" -i "${hit_beds[@]}" > "$MULTI"
for support in 1 2 3; do
  awk -v n="$support" 'BEGIN{OFS="\t"} NR>1 && $4>=n {print $1,$2,$3}' "$MULTI" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "$OUT_ROOT/analysis/Potential_G4_hit_by_${support}of3_targets.bed"
done

echo "三个靶点的Potential G4比较已汇总：$SUMMARY"
