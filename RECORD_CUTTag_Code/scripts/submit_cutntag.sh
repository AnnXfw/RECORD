#!/usr/bin/env bash
# 按依赖关系一次性提交CUT&Tag全流程
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT/scripts"
mkdir -p "$PROJECT/logs"

align_job=$(sbatch --parsable cutntag_array.sbatch)
qc_job=$(sbatch --parsable --dependency="afterok:${align_job}" cutntag_collect.sbatch)
peak_job=$(sbatch --parsable --dependency="afterok:${align_job}" cutntag_peaks.sbatch)
consensus_job=$(sbatch --parsable --dependency="afterok:${peak_job}" cutntag_consensus.sbatch)
g4_job=$(sbatch --parsable --dependency="afterok:${consensus_job}" cutntag_g4_analysis.sbatch)
g4_collect_job=$(sbatch --parsable --dependency="afterok:${g4_job}" cutntag_g4_collect.sbatch)

printf '比对数组作业: %s\n' "$align_job"
printf '质控汇总作业: %s\n' "$qc_job"
printf 'MACS3找峰作业: %s\n' "$peak_job"
printf '三重复共识作业: %s\n' "$consensus_job"
printf 'G4下游分析作业: %s\n' "$g4_job"
printf '三个靶点G4汇总作业: %s\n' "$g4_collect_job"
printf 'HOMER尚未自动提交；安装HOMER和hg19数据后再单独提交cutntag_homer.sbatch。\n'
