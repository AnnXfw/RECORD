#!/bin/bash
#SBATCH --job-name=peak_annotation
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/peak_annotation_%j.out
#SBATCH --error=logs/peak_annotation_%j.err

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT_DIR="$ROOT/scripts"
RESULTS="$ROOT/results"
REFERENCE="$ROOT/reference"
ENV="${GENOMEPLOT_ENV:-genomeplot}"

module purge
module load miniconda3/24.3.0
source activate "$ENV"

python "$SCRIPT_DIR/annotate_peak_genomic_distribution.py" \
    --gtf "$REFERENCE/hg19.refGene.gtf.gz" \
    --klf16-peaks "$RESULTS/consensus/KLF16/KLF16.consensus.2of3.bed" \
    --hey1-peaks "$RESULTS/diagnostic_peak_calling/HEY1/default/HEY1.default.consensus.2of3.bed" \
    --bg4-peaks "$REFERENCE/G4.bed" \
    --output-dir "$RESULTS/genomic_distribution"

ls -lh \
    "$RESULTS/genomic_distribution/KLF16.genomic_distribution.tsv" \
    "$RESULTS/genomic_distribution/KLF16.genomic_distribution.svg" \
    "$RESULTS/genomic_distribution/HEY1.genomic_distribution.tsv" \
    "$RESULTS/genomic_distribution/HEY1.genomic_distribution.svg"
