#!/bin/bash
#SBATCH --job-name=g4_loci
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-2%2
#SBATCH --output=logs/g4_loci_%A_%a.out
#SBATCH --error=logs/g4_loci_%A_%a.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
GENOMEPLOT_ENV="${GENOMEPLOT_ENV:-genomeplot}"

module purge
module load miniconda3/24.3.0
source activate "$GENOMEPLOT_ENV"

case "${SLURM_ARRAY_TASK_ID}" in
    1) TARGET=KLF16 ;;
    2) TARGET=HEY1 ;;
    *) echo "错误的数组任务编号" >&2; exit 2 ;;
esac

python "$SCRIPT_DIR/select_and_plot_g4_loci.py" \
    --target "$TARGET" \
    --root "$ROOT" \
    --top-n 3 \
    --window 20000
