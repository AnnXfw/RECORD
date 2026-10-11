# RECORD proteomics analysis code

This package contains the R scripts used for retron-based and DRT10-based RECORD proteomics analyses and volcano-plot generation.

## Contents

- `proteomics/retron_record/Retron1.9_analysis.R`: original retron RECORD analysis.
- `proteomics/retron_record/Retron_volcano_clean.R`: cleaned analysis workflow.
- `proteomics/retron_record/Retron_volcano_final.R`: final manuscript volcano-plot workflow.
- `proteomics/drt10_record/DRT10_analysis.R`: original DRT10 RECORD analysis.
- `proteomics/drt10_record/DRT10_volcano_final.R`: final DRT10 volcano-plot workflow.

## Directory layout

- Place retron DIA-NN input tables in `data/proteomics/retron_record/`.
- Place DRT10 DIA-NN input tables in `data/proteomics/drt10_record/`.
- Scripts write processed tables and figures to the matching directory under `results/proteomics/`.

## Running the scripts

Run scripts from the package root. By default, the scripts use the current working directory as the project root. Alternatively, set the `RECORD_CODE_ROOT` environment variable to the absolute path of this package.

Before deposition, add the R and package versions to `environment/versions_to_report.txt` and confirm that the deposited input filenames match those referenced by the scripts.


## Input sample headers

Retron scripts normalize Windows or Unix paths in DIA-NN column headers to raw filenames before matching WT-1.raw through WT-3.raw, TELO-1.raw through TELO-3.raw, and MUT-1.raw through MUT-3.raw. Input files do not need to retain the original computer paths. Duplicate normalized headers are rejected.
