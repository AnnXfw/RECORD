#!/usr/bin/env bash
# HEY1和ING4诊断性重新找峰：比较pooled control的两种MACS3缩放策略
# 该作业只读取现有final BAM，结果写入diagnostic_peak_calling，不覆盖正式结果
#SBATCH --job-name=cuttag_diag_peaks
#SBATCH --partition=cpu
#SBATCH --output=logs/diag_peaks_%A_%a.out
#SBATCH --error=logs/diag_peaks_%A_%a.err
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

targets=(HEY1 ING4)
target=${targets[$((SLURM_ARRAY_TASK_ID-1))]}
DIAG="$OUT_ROOT/diagnostic_peak_calling/$target"
mkdir -p "$DIAG"

# 合并同一靶点的三个No-primary对照，增加对照有效深度
controls=()
for rep in 1 2 3; do
    control="$OUT_ROOT/bam/${target}_C${rep}/${target}_C${rep}.final.bam"
    [[ -s "$control" ]] || { echo "缺少对照BAM：$control" >&2; exit 2; }
    controls+=("$control")
done

POOLED="$DIAG/${target}.pooled_no_primary.bam"
samtools merge -f -@ "$SLURM_CPUS_PER_TASK" "$POOLED" "${controls[@]}"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$POOLED"
samtools flagstat "$POOLED" > "$DIAG/${target}.pooled_no_primary.flagstat.txt"
pooled_fragments=$(awk 'NR==1{printf "%.0f",$1/2}' "$DIAG/${target}.pooled_no_primary.flagstat.txt")

SUMMARY="$DIAG/${target}.diagnostic_peak_counts.tsv"
CONSENSUS_SUMMARY="$DIAG/${target}.diagnostic_consensus_counts.tsv"
printf 'target\treplicate\tstrategy\ttreatment_fragments\tpooled_control_fragments\tpeak_count\n' > "$SUMMARY"
printf 'target\tstrategy\tnonempty_replicates\tconsensus_2of3_count\n' > "$CONSENSUS_SUMMARY"

for strategy in default scale_large; do
    strategy_dir="$DIAG/$strategy"
    mkdir -p "$strategy_dir"
    nonempty_peaks=()

    for rep in 1 2 3; do
        sample="${target}_${rep}"
        treatment="$OUT_ROOT/bam/$sample/${sample}.final.bam"
        flagstat="$OUT_ROOT/bam/$sample/${sample}.flagstat.txt"
        [[ -s "$treatment" ]] || { echo "缺少实验BAM：$treatment" >&2; exit 3; }
        [[ -s "$flagstat" ]] || { echo "缺少flagstat：$flagstat" >&2; exit 4; }

        name="${target}_rep${rep}_${strategy}"
        outdir="$strategy_dir/${target}_rep${rep}"
        mkdir -p "$outdir"
        scale_args=()
        if [[ "$strategy" == "scale_large" ]]; then
            scale_args=(--scale-to large)
        fi

        macs3 callpeak \
          -t "$treatment" -c "$POOLED" -f BAMPE -g hs \
          -n "$name" --outdir "$outdir" --qvalue 0.05 \
          --keep-dup all --call-summits "${scale_args[@]}"

        raw_peak="$outdir/${name}_peaks.narrowPeak"
        sorted_peak="$outdir/${name}.analysis.narrowPeak"
        [[ -f "$raw_peak" ]] || { echo "MACS3未生成peak文件：$raw_peak" >&2; exit 5; }
        LC_ALL=C sort -k1,1 -k2,2n "$raw_peak" > "$sorted_peak"
        peak_count=$(wc -l < "$sorted_peak")
        treatment_fragments=$(awk 'NR==1{printf "%.0f",$1/2}' "$flagstat")
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$target" "$rep" "$strategy" "$treatment_fragments" \
          "$pooled_fragments" "$peak_count" >> "$SUMMARY"

        if [[ -s "$sorted_peak" ]]; then
            nonempty_peaks+=("$sorted_peak")
        fi
    done

    # 如果至少两个重复产生peak，则计算任意两个重复共同支持的区域
    multiinter="$strategy_dir/${target}.${strategy}.multiinter.bed"
    consensus="$strategy_dir/${target}.${strategy}.consensus.2of3.bed"
    : > "$multiinter"
    : > "$consensus"
    if [[ "${#nonempty_peaks[@]}" -ge 2 ]]; then
        bedtools multiinter -i "${nonempty_peaks[@]}" > "$multiinter"
        awk 'BEGIN{OFS="\t"} $4>=2 {print $1,$2,$3}' "$multiinter" \
          | LC_ALL=C sort -k1,1 -k2,2n \
          | bedtools merge -i - > "$consensus"
    fi
    consensus_count=$(wc -l < "$consensus")
    printf '%s\t%s\t%s\t%s\n' \
      "$target" "$strategy" "${#nonempty_peaks[@]}" "$consensus_count" \
      >> "$CONSENSUS_SUMMARY"
done

echo "$target 诊断性peak calling完成"
column -t -s $'\t' "$SUMMARY" || cat "$SUMMARY"
column -t -s $'\t' "$CONSENSUS_SUMMARY" || cat "$CONSENSUS_SUMMARY"
