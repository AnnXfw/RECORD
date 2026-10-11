# RECORD analysis code

Custom analysis and figure-generation code for the RECORD study.

## Workflows

- [Proteomics](RECORD_Proteomics_Code/README.md): retron- and DRT10-based proteomics analyses and volcano plots.
- [CUT&Tag and BG4](RECORD_CUTTag_Code/README.md): read processing, peak calling, reproducible peaks, overlap, permutation and genomic visualization.

Run each workflow from its own package directory and follow its README. Raw experimental data and processed input tables are not included. Required inputs must be supplied separately; placeholder files mark their locations. Reference genome and BG4 datasets must match the assembly used in the scripts.

## Reproducibility status

Private machine paths have been removed from the packaged scripts. R and Python syntax checks passed; full workflows were not rerun during packaging. The environment files record available tool-version information, but R/Python package versions and the original HPC environment remain to be confirmed. Do not interpret missing version entries as verified dependencies.

## License

MIT; see [LICENSE](LICENSE). Third-party software and datasets retain their own terms.

## Archival release

The planned initial release is v1.0.0. A Zenodo DOI will be added after archival; no DOI has yet been assigned.
