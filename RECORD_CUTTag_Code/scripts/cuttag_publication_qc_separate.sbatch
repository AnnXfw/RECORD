#!/bin/bash
#SBATCH --job-name=cuttag_qc_sep
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/cuttag_qc_sep_%A_%a.out
#SBATCH --error=logs/cuttag_qc_sep_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
SAMPLES="$SCRIPT_DIR/samples.tsv"
CUTTAG_ENV="${CUTTAG_ENV:-cuttag}"
GENOMEPLOT_ENV="${GENOMEPLOT_ENV:-genomeplot}"
BLACKLIST="$ROOT/reference/hg19-blacklist.v2.bed"

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
        echo "错误的数组任务编号" >&2
        exit 2
        ;;
esac

OUT="$RESULTS/publication_qc_separate/$TARGET"
mkdir -p "$OUT"

module purge
module load miniconda3/24.3.0
source activate "$CUTTAG_ENV"

SUMMARY="$OUT/${TARGET}.sample_qc.tsv"
printf "sample\ttarget\tcondition\treplicate\tmapped_reads\talignment_rate\tpeak_count\treads_in_peaks\tfrip_percent\n" > "$SUMMARY"

declare -a BAMS=()
declare -a LABELS=()

while IFS=$'\t' read -r sample target condition replicate batch r1 r2; do
    [[ "$sample" == "sample" ]] && continue
    [[ "$target" == "$TARGET" ]] || continue

    bam="$RESULTS/bam/$sample/${sample}.final.bam"
    bowtie_log="$RESULTS/bam/$sample/${sample}.bowtie2.log"
    [[ -s "$bam" ]] || { echo "缺少BAM：$bam" >&2; exit 3; }
    [[ -s "${bam}.bai" ]] || { echo "缺少BAM索引：${bam}.bai" >&2; exit 4; }
    [[ -s "$PEAKS" ]] || { echo "缺少峰文件：$PEAKS" >&2; exit 5; }
    [[ -s "$bowtie_log" ]] || { echo "缺少Bowtie2日志：$bowtie_log" >&2; exit 6; }

    mapped=$(samtools view -c -F 260 "$bam")
    in_peaks=$(samtools view -c -F 260 -L "$PEAKS" "$bam")
    peak_count=$(awk 'NF>=3 && $1 !~ /^(#|track|browser)/ {n++} END{print n+0}' "$PEAKS")
    alignment=$(awk '/overall alignment rate/{gsub("%","",$1); value=$1} END{print value}' "$bowtie_log")
    frip=$(awk -v numerator="$in_peaks" -v denominator="$mapped" \
        'BEGIN{if(denominator>0) printf "%.6f",100*numerator/denominator; else print "NA"}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$target" "$condition" "$replicate" \
        "$mapped" "$alignment" "$peak_count" "$in_peaks" "$frip" >> "$SUMMARY"

    BAMS+=("$bam")
    LABELS+=("$sample")
done < "$SAMPLES"

[[ "${#BAMS[@]}" -eq 6 ]] || {
    echo "预期${TARGET}有6个BAM，实际${#BAMS[@]}" >&2
    exit 7
}

blacklist_args=()
if [[ -s "$BLACKLIST" ]]; then
    blacklist_args=(--blackListFileName "$BLACKLIST")
fi

multiBamSummary bins \
    --bamfiles "${BAMS[@]}" \
    --labels "${LABELS[@]}" \
    "${blacklist_args[@]}" \
    --binSize 10000 \
    --outFileName "$OUT/${TARGET}.multiBamSummary.npz" \
    --outRawCounts "$OUT/${TARGET}.bam_counts.tsv" \
    --numberOfProcessors "$SLURM_CPUS_PER_TASK"

plotCorrelation \
    -in "$OUT/${TARGET}.multiBamSummary.npz" \
    --corMethod spearman \
    --whatToPlot heatmap \
    --plotNumbers \
    --removeOutliers \
    --outFileCorMatrix "$OUT/${TARGET}.spearman_correlation.tsv" \
    -o "$OUT/${TARGET}.spearman_correlation.svg"

plotPCA \
    -in "$OUT/${TARGET}.multiBamSummary.npz" \
    -o "$OUT/${TARGET}.PCA.svg" \
    --outFileNameData "$OUT/${TARGET}.PCA_data.tsv"

source activate "$GENOMEPLOT_ENV"
python "$SCRIPT_DIR/plot_cuttag_qc_summary.py" \
    --input "$SUMMARY" \
    --output-dir "$OUT" \
    --target "$TARGET"

ls -lh "$OUT"
