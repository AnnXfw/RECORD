#!/usr/bin/env bash
# 使用HEY1 pooled No-primary默认缩放得到的105个2-of-3共识峰进行G4分析
# 结果写入独立目录，不覆盖原始HEY1分析
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$PROJECT/scripts"
CONSENSUS="$PROJECT/results/diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed"
ANALYSIS="$PROJECT/results/diagnostic_g4/HEY1_pooled_default"

[[ -s "$CONSENSUS" ]] || { echo "缺少HEY1的105个共识峰：$CONSENSUS" >&2; exit 2; }
mkdir -p "$PROJECT/logs" "$ANALYSIS"
cd "$SCRIPT_DIR"

job_id=$(sbatch --parsable --array=1 \
  --export="ALL,TARGET_OVERRIDE=HEY1,CONSENSUS_OVERRIDE=$CONSENSUS,ANALYSIS_OVERRIDE=$ANALYSIS" \
  cutntag_g4_analysis.sbatch)

echo "HEY1 pooled-default G4作业：$job_id"
echo "输入共识峰：$CONSENSUS"
echo "输出目录：$ANALYSIS"
