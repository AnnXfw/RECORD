#!/bin/bash
#SBATCH --job-name=cuttag_pub_qc
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cuttag_pub_qc_%j.out
#SBATCH --error=logs/cuttag_pub_qc_%j.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
OUT="$RESULTS/publication_qc"
SAMPLES="$SCRIPT_DIR/samples.tsv"
CUTTAG_ENV="${CUTTAG_ENV:-cuttag}"
GENOMEPLOT_ENV="${GENOMEPLOT_ENV:-genomeplot}"
BLACKLIST="$ROOT/reference/hg19-blacklist.v2.bed"

module purge
module load miniconda3/24.3.0
source activate "$CUTTAG_ENV"
mkdir -p "$OUT"

SUMMARY="$OUT/KLF16_HEY1.sample_qc.tsv"
printf "sample\ttarget\tcondition\treplicate\tmapped_reads\talignment_rate\tpeak_count\treads_in_peaks\tfrip_percent\n" > "$SUMMARY"

declare -a BAMS=()
declare -a LABELS=()

while IFS=$'\t' read -r sample target condition replicate batch r1 r2; do
    [[ "$sample" == "sample" ]] && continue
    [[ "$target" == "KLF16" || "$target" == "HEY1" ]] || continue

    bam="$RESULTS/bam/$sample/${sample}.final.bam"
    bowtie_log="$RESULTS/bam/$sample/${sample}.bowtie2.log"
    if [[ "$target" == "KLF16" ]]; then
        peaks="$RESULTS/consensus/KLF16/KLF16.consensus.2of3.bed"
    else
        peaks="$RESULTS/diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed"
    fi

    [[ -s "$bam" ]] || { echo "缺少BAM：$bam" >&2; exit 3; }
    [[ -s "${bam}.bai" ]] || { echo "缺少BAM索引：${bam}.bai" >&2; exit 4; }
    [[ -s "$peaks" ]] || { echo "缺少峰文件：$peaks" >&2; exit 5; }
    [[ -s "$bowtie_log" ]] || { echo "缺少Bowtie2日志：$bowtie_log" >&2; exit 6; }

    mapped=$(samtools view -c -F 260 "$bam")
    in_peaks=$(samtools view -c -F 260 -L "$peaks" "$bam")
    peak_count=$(awk 'NF>=3 && $1 !~ /^(#|track|browser)/ {n++} END{print n+0}' "$peaks")
    alignment=$(awk '/overall alignment rate/{gsub("%","",$1); value=$1} END{print value}' "$bowtie_log")
    frip=$(awk -v numerator="$in_peaks" -v denominator="$mapped" \
        'BEGIN{if(denominator>0) printf "%.6f",100*numerator/denominator; else print "NA"}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$target" "$condition" "$replicate" \
        "$mapped" "$alignment" "$peak_count" "$in_peaks" "$frip" >> "$SUMMARY"

    BAMS+=("$bam")
    LABELS+=("$sample")
done < "$SAMPLES"

[[ "${#BAMS[@]}" -eq 12 ]] || {
    echo "预期KLF16和HEY1共12个BAM，实际${#BAMS[@]}" >&2
    exit 7
}

# 12个正式分析样本的全基因组bin相关性与PCA
blacklist_args=()
if [[ -s "$BLACKLIST" ]]; then
    blacklist_args=(--blackListFileName "$BLACKLIST")
fi

multiBamSummary bins \
    --bamfiles "${BAMS[@]}" \
    --labels "${LABELS[@]}" \
    "${blacklist_args[@]}" \
    --binSize 10000 \
    --outFileName "$OUT/KLF16_HEY1.multiBamSummary.npz" \
    --outRawCounts "$OUT/KLF16_HEY1.bam_counts.tsv" \
    --numberOfProcessors "$SLURM_CPUS_PER_TASK"

plotCorrelation \
    -in "$OUT/KLF16_HEY1.multiBamSummary.npz" \
    --corMethod spearman \
    --whatToPlot heatmap \
    --plotNumbers \
    --removeOutliers \
    --outFileCorMatrix "$OUT/KLF16_HEY1.spearman_correlation.tsv" \
    -o "$OUT/KLF16_HEY1.spearman_correlation.svg"

plotPCA \
    -in "$OUT/KLF16_HEY1.multiBamSummary.npz" \
    -o "$OUT/KLF16_HEY1.PCA.svg" \
    --outFileNameData "$OUT/KLF16_HEY1.PCA_data.tsv"

# 使用绘图环境输出每个靶点独立QC概览SVG
source activate "$GENOMEPLOT_ENV"
python "$SCRIPT_DIR/plot_cuttag_qc_summary.py" \
    --input "$SUMMARY" \
    --output-dir "$OUT"

ls -lh \
    "$SUMMARY" \
    "$OUT/KLF16.QC_summary.svg" \
    "$OUT/HEY1.QC_summary.svg" \
    "$OUT/KLF16_HEY1.spearman_correlation.svg" \
    "$OUT/KLF16_HEY1.PCA.svg"
