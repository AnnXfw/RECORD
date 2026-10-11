#!/usr/bin/env bash
# CUT&Tag显著富集共识peak与Potential G4/BG4的交集及随机置换检验
# 数组1：KLF16正式2-of-3共识峰；数组2：HEY1 pooled-default探索性共识峰
# 随机峰保持原峰长度和所在染色体，重复1000次
#SBATCH --job-name=cuttag_g4_enrich
#SBATCH --partition=cpu
#SBATCH --output=logs/g4_enrich_%A_%a.out
#SBATCH --error=logs/g4_enrich_%A_%a.err
#SBATCH --array=1-2%1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
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

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHROM_SIZES="$PROJECT/reference/hg19.chrom.sizes"
N_PERM=1000

if [[ "$SLURM_ARRAY_TASK_ID" -eq 1 ]]; then
    target=KLF16
    result_label=KLF16_strict_2of3
    PEAKS="$OUT_ROOT/consensus/KLF16/KLF16.consensus.2of3.bed"
else
    target=HEY1
    result_label=HEY1_pooled_default_2of3
    PEAKS="$OUT_ROOT/diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed"
fi

OUT="$OUT_ROOT/g4_enrichment/$result_label"
mkdir -p "$OUT"
[[ -s "$PEAKS" ]] || { echo "缺少共识peak：$PEAKS" >&2; exit 2; }
[[ -s "$POTENTIAL_G4_BED" ]] || { echo "缺少Potential G4：$POTENTIAL_G4_BED" >&2; exit 3; }
[[ -s "$BG4_BED" ]] || { echo "缺少BG4：$BG4_BED" >&2; exit 4; }
[[ -s "$CHROM_SIZES" ]] || { echo "缺少hg19染色体长度：$CHROM_SIZES" >&2; exit 5; }

TMP_BASE=${SLURM_TMPDIR:-/tmp/${USER}_g4_enrich_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}
mkdir -p "$TMP_BASE"
trap 'rm -rf -- "$TMP_BASE"' EXIT

# 所有BED统一排序；正式输出只保留前三列
PEAKS_SORTED="$TMP_BASE/peaks.sorted.bed"
POTENTIAL_SORTED="$TMP_BASE/Potential_G4.sorted.bed"
BG4_SORTED="$TMP_BASE/BG4.sorted.bed"
LC_ALL=C sort -k1,1 -k2,2n "$PEAKS" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > "$PEAKS_SORTED"
LC_ALL=C sort -k1,1 -k2,2n "$POTENTIAL_G4_BED" > "$POTENTIAL_SORTED"
LC_ALL=C sort -k1,1 -k2,2n "$BG4_BED" > "$BG4_SORTED"

peak_n=$(wc -l < "$PEAKS_SORTED")

# 真实观察值：以CUT&Tag peak为分母，每个peak只计数一次
bedtools intersect -u -a "$PEAKS_SORTED" -b "$POTENTIAL_SORTED" \
  > "$OUT/${result_label}_with_Potential_G4.bed"
bedtools intersect -v -a "$PEAKS_SORTED" -b "$POTENTIAL_SORTED" \
  > "$OUT/${result_label}_without_Potential_G4.bed"
bedtools intersect -u -a "$PEAKS_SORTED" -b "$BG4_SORTED" \
  > "$OUT/${result_label}_with_BG4.bed"
bedtools intersect -v -a "$PEAKS_SORTED" -b "$BG4_SORTED" \
  > "$OUT/${result_label}_without_BG4.bed"

observed_potential=$(wc -l < "$OUT/${result_label}_with_Potential_G4.bed")
observed_bg4=$(wc -l < "$OUT/${result_label}_with_BG4.bed")

# 反向输出：哪些G4位点被CUT&Tag peak覆盖
bedtools intersect -u -a "$POTENTIAL_SORTED" -b "$PEAKS_SORTED" \
  > "$OUT/Potential_G4_hit_by_${result_label}.bed"
bedtools intersect -u -a "$BG4_SORTED" -b "$PEAKS_SORTED" \
  > "$OUT/BG4_hit_by_${result_label}.bed"

NULL="$OUT/${result_label}.permutation_counts.tsv"
printf 'iteration\tPotential_G4_overlap_peaks\tBG4_overlap_peaks\n' > "$NULL"

# blacklist存在时禁止随机峰落入blacklist；不存在时仍可完成检验
shuffle_extra=()
if [[ -n "${BLACKLIST:-}" && -s "$BLACKLIST" ]]; then
    shuffle_extra=(-excl "$BLACKLIST")
fi

SHUFFLED="$TMP_BASE/shuffled.bed"
for i in $(seq 1 "$N_PERM"); do
    bedtools shuffle -i "$PEAKS_SORTED" -g "$CHROM_SIZES" \
      -chrom -noOverlapping -maxTries 10000 -seed "$i" "${shuffle_extra[@]}" \
      > "$SHUFFLED"
    random_potential=$(bedtools intersect -u -a "$SHUFFLED" -b "$POTENTIAL_SORTED" | wc -l)
    random_bg4=$(bedtools intersect -u -a "$SHUFFLED" -b "$BG4_SORTED" | wc -l)
    printf '%s\t%s\t%s\n' "$i" "$random_potential" "$random_bg4" >> "$NULL"
done

# 计算随机均值、标准差、fold enrichment和经验P值；加1校正避免P=0
SUMMARY="$OUT/${result_label}.g4_enrichment_summary.tsv"
printf 'target\tpeak_set\tg4_set\tpeak_count\tobserved_overlap\toverlap_percent\tnull_mean\tnull_sd\tfold_enrichment\tempirical_p\tpermutations\n' > "$SUMMARY"

summarize_null() {
    local column=$1
    local observed=$2
    local g4_name=$3
    awk -v col="$column" -v obs="$observed" -v npeak="$peak_n" \
        -v target="$target" -v label="$result_label" -v g4="$g4_name" -v nperm="$N_PERM" \
      'BEGIN{OFS="\t"}
       NR>1 {x=$col; sum+=x; sumsq+=x*x; if(x>=obs) ge++}
       END {
         n=NR-1;
         mean=(n>0?sum/n:0);
         variance=(n>1?(sumsq-sum*sum/n)/(n-1):0);
         sd=(variance>0?sqrt(variance):0);
         pct=(npeak>0?100*obs/npeak:0);
         fold=(mean>0?obs/mean:"NA");
         p=(ge+1)/(n+1);
         printf "%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%s\t%.6g\t%d\n", target,label,g4,npeak,obs,pct,mean,sd,fold,p,n;
       }' "$NULL" >> "$SUMMARY"
}

summarize_null 2 "$observed_potential" Potential_G4
summarize_null 3 "$observed_bg4" BG4

echo "$result_label 的G4富集检验完成"
cat "$SUMMARY"
