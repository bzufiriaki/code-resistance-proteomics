# Input data

Data come from the Alzheimer's Disease Neuroimaging Initiative (ADNI) and
**cannot be redistributed** in this repository. Request access at
<https://adni.loni.usc.edu>. Everything in `data/` except this file is ignored by git.

Each script takes file paths as options, so you can keep the data anywhere.
The expected files and the columns the scripts use are listed below.

## 01_calculate_resistance.R

| File | Required columns |
|---|---|
| `ADNIMERGE_<date>.csv` (ADNI) | `RID`, `APOE4`, `PTETHCAT`, `PTRACCAT`, `DX_bl` |
| Subject-level dataset (e.g. `All_residuals_with_median.csv`) | `RID`, `AGE`, `PTEDUCAT`, `PTGENDER`, `Status` (CN/MCI/AD) and the outcome (`ABETA_adj`, `PTAU_adj`, …) |
