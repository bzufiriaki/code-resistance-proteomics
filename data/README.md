# Input data

Data for this project come from the Alzheimer's Disease Neuroimaging Initiative (ADNI) and
**cannot be redistributed** in this repository. Request access at
<https://adni.loni.usc.edu>. Everything in `data/` except this file is ignored by git.

Each script takes file paths as options, so you can keep the data anywhere.
The expected files and the columns the scripts use are listed below.

## 01_A_MORE.R

| File | Required columns |
|---|---|
| Subject-level dataset (e.g. `individuals.csv`) | `RID`, `AGE`, `PTEDUCAT`, `PTGENDER`, `Cognitive_Status` (CN/MCI/AD), 'APOE4' and the outcome (`ABETA`, `PTAU`, …) |
