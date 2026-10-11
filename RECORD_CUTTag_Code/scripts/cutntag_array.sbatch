#!/usr/bin/env bash
# 作业名称
#SBATCH --job-name=cuttag
# 使用可访问Lustre的cpu分区
#SBATCH --partition=cpu
# 日志使用绝对路径，避免从不同目录提交时找不到logs
#SBATCH --output=logs/cuttag_%A_%a.out
#SBATCH --error=logs/cuttag_%A_%a.err
# 18个样本，最多同时运行6个
#SBATCH --array=1-18%6
# 每个样本使用12个CPU线程
#SBATCH --cpus-per-task=12
# 每个样本申请32 GB内存
#SBATCH --mem=32G
# 单个样本最长运行24小时
#SBATCH --time=24:00:00
# 不继承提交节点环境，保证环境可控
#SBATCH --export=NONE
set -euo pipefail

# 读取超算路径和Conda环境配置
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="${ROOT}/scripts"
source "$SCRIPT_DIR/config.env"
# 创建输出目录：质控、BAM、peak、bigWig和日志
mkdir -p "$OUT_ROOT" "$OUT_ROOT/fastqc" "$OUT_ROOT/fastp" "$OUT_ROOT/bam" "$OUT_ROOT/bigWig" "$OUT_ROOT/logs"
# 显式加载超算Miniconda模块并激活CUT&Tag环境
module purge
module load miniconda3/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# 根据数组编号从样本表读取一个样本
line=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print; exit}' "$SAMPLES_TSV")
IFS=$'\t' read -r sample target condition replicate batch r1 r2 <<< "$line"
[[ -n "$sample" ]] || { echo "No sample for array index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }
# FASTQ在rawdata目录中平铺存放。样本表保留原始批次路径用于追溯，
# 运行时只取文件名；RE重测文件名会被完整保留并正确映射到C1/C3样本。
R1="$RAW_ROOT/$(basename "$r1")"
R2="$RAW_ROOT/$(basename "$r2")"
[[ -s "$R1" && -s "$R2" ]] || { echo "Missing FASTQ: $R1 $R2" >&2; exit 3; }

# 设置当前样本的原始质控、过滤和比对结果目录
Q="$OUT_ROOT/fastqc/$sample"; F="$OUT_ROOT/fastp/$sample"; B="$OUT_ROOT/bam/$sample"
mkdir -p "$Q" "$F" "$B"
# 对原始双端FASTQ执行FastQC
fastqc -t 2 -o "$Q" "$R1" "$R2"
# PE150双端接头和低质量过滤
fastp -i "$R1" -I "$R2" -o "$F/${sample}.R1.fq.gz" -O "$F/${sample}.R2.fq.gz" \
  --detect_adapter_for_pe --trim_poly_g --length_required 30 --thread "$SLURM_CPUS_PER_TASK" \
  --html "$F/${sample}.fastp.html" --json "$F/${sample}.fastp.json"

# 使用Bowtie2将双端reads比对到hg19；BT2_INDEX在config.env中指定
bowtie2 --very-sensitive --no-mixed --no-discordant --dovetail -I 10 -X 1000 -p "$SLURM_CPUS_PER_TASK" \
  -x "$BT2_INDEX" -1 "$F/${sample}.R1.fq.gz" -2 "$F/${sample}.R2.fq.gz" \
  2> "$B/${sample}.bowtie2.log" | samtools view -b -o "$B/${sample}.name.bam" -
# name sort、fixmate、坐标排序，为samtools markdup准备输入
samtools sort -n -@ "$SLURM_CPUS_PER_TASK" -o "$B/${sample}.namesort.bam" "$B/${sample}.name.bam"
samtools fixmate -m "$B/${sample}.namesort.bam" "$B/${sample}.fixmate.bam"
samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$B/${sample}.positionsort.bam" "$B/${sample}.fixmate.bam"
# 去除PCR重复，同时输出重复率统计
samtools markdup -r -s -f "$B/${sample}.markdup.stats.txt" -@ "$SLURM_CPUS_PER_TASK" "$B/${sample}.positionsort.bam" "$B/${sample}.dedup.bam"
# 保留properly paired、MAPQ>=30且排除异常flag的reads
samtools view -@ "$SLURM_CPUS_PER_TASK" -b -q 30 -f 2 -F 3844 "$B/${sample}.dedup.bam" -o "$B/${sample}.filtered.bam"
samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$B/${sample}.final.bam" "$B/${sample}.filtered.bam"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$B/${sample}.final.bam"
# 输出最终BAM统计信息
samtools flagstat "$B/${sample}.final.bam" > "$B/${sample}.flagstat.txt"
samtools stats "$B/${sample}.final.bam" > "$B/${sample}.stats.txt"
samtools idxstats "$B/${sample}.final.bam" > "$B/${sample}.idxstats.txt"

# 按原分析流程生成RPKM标准化bigWig，窗口大小为50 bp
bamCoverage -b "$B/${sample}.final.bam" -o "$OUT_ROOT/bigWig/${sample}.bs50.RPKM.bw" \
  --normalizeUsing RPKM --binSize 50 --numberOfProcessors "$SLURM_CPUS_PER_TASK"
# 删除中间BAM，仅保留final.bam、索引和统计结果以节省空间
rm -f "$B/${sample}.name.bam" "$B/${sample}.namesort.bam" "$B/${sample}.fixmate.bam" "$B/${sample}.positionsort.bam" "$B/${sample}.dedup.bam" "$B/${sample}.filtered.bam"
