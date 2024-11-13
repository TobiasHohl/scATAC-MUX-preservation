# Analysis

## 1. Create conda environment
- Use the environment yaml `bulk_analysis.yaml` to install all dependencies via `conda env create -f bulk_analysis.yml`.
- Alternatively, install the required packages (see yaml) manually. Most of them can also be installed via pip.

## 2. Run analysis scripts
- Activate the respective conda environment with `conda activate bulk_analysis`.
- Modify the required parameter in the script `02_analysis_runthis.sh`.
- If `pygenometracks` should be run, modify the `.ini` files under `./scripts/pygenometracks/` accordingly and un-comment line 44 of `02_analysis_runthis.sh`.
- Run `02_analysis_runthis.sh` via:
```
bash 02_analysis_runthis.sh
```

## 3. Plotting with R
Follow the `03_bulkFigures.Rmd` and modify paths accordingly.