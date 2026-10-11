# RECORD CUT&Tag and BG4 analysis code

This package contains the code used for CUT&Tag processing and comparison with published BG4 and G4-seq datasets.

## Layout

- `scripts/`: shell, SLURM and Python scripts, together with `config.env`, `samples.tsv` and `peak_pairs.tsv`.
- `data/cuttag_bg4/`: raw-data links and reference resources.
- `results/cuttag_bg4/`: processed outputs and figures.
- `logs/`: SLURM logs.
- `environment/`: software-version record.

Submit SLURM jobs from the package root so that the relative `logs/` paths resolve correctly. The scripts derive the package root from their own location. Paths in `scripts/config.env` can be overridden with environment variables when data or references are stored elsewhere.

The default conda environment is `cuttag`; plotting scripts that use `genomeplot` can be redirected with `GENOMEPLOT_ENV`. Before release, complete the environment-version file and add the sequencing repository accession to this README.
