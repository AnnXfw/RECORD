#!/usr/bin/env bash
# 汇总质控作业
#SBATCH --job-name=cuttag_qc
# 使用可访问Lustre的cpu分区
#SBATCH --partition=cpu
# 日志使用绝对路径
#SBATCH --output=logs/cuttag_qc_%j.out
#SBATCH --error=logs/cuttag_qc_%j.err
# 汇总任务使用8个CPU
#SBATCH --cpus-per-task=8
# 申请24 GB内存
#SBATCH --mem=24G
# 最长运行4小时
#SBATCH --time=04:00:00
# 不继承提交节点环境
#SBATCH --export=NONE
set -euo pipefail
# 读取与数组任务相同的配置
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"
mkdir -p "$OUT_ROOT/qc" "$OUT_ROOT/logs"
module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
# 汇总原始FastQC、fastp、BAM统计和比对日志
multiqc "$OUT_ROOT/fastqc" "$OUT_ROOT/fastp" "$OUT_ROOT/bam" -o "$OUT_ROOT/qc" -n multiqc.html
mapfile -t bams < <(find "$OUT_ROOT/bam" -name '*.final.bam' | sort)
[[ "${#bams[@]}" -eq 18 ]] || { echo "Expected 18 final BAM files, found ${#bams[@]}" >&2; exit 3; }
labels=()
for bam in "${bams[@]}"; do labels+=("$(basename "$bam" .final.bam)"); done
printf '%s\n' "${bams[@]}" > "$OUT_ROOT/qc/bams.list"
# 计算样本全基因组bin覆盖度，用于相关性和PCA
blacklist_args=()
if [[ -n "${BLACKLIST:-}" && -s "$BLACKLIST" ]]; then
  blacklist_args=(--blackListFileName "$BLACKLIST")
fi
multiBamSummary bins --bamfiles "${bams[@]}" --labels "${labels[@]}" \
  "${blacklist_args[@]}" --outFileName "$OUT_ROOT/qc/multiBamSummary.npz" \
  --outRawCounts "$OUT_ROOT/qc/bam_counts.tsv" --numberOfProcessors "$SLURM_CPUS_PER_TASK"
# 输出样本相关性热图和PCA图
plotCorrelation -in "$OUT_ROOT/qc/multiBamSummary.npz" --corMethod spearman --whatToPlot heatmap \
  --plotNumbers --outFileCorMatrix "$OUT_ROOT/qc/sample_correlation.tsv" -o "$OUT_ROOT/qc/sample_correlation.svg"
plotPCA -in "$OUT_ROOT/qc/multiBamSummary.npz" -o "$OUT_ROOT/qc/sample_PCA.svg"
