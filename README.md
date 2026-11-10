# Proteomic signatures of resistance to Alzheimer's disease pathology

Regarding the manuscript "Molecular signatures of resistance to amyloid and tau pathology identify distinct pathways of protection and vulnerability in Alzheimer’s disease", here we provide the code for quantifying individual **resistance** to amyloid and tau pathology.

We developed a-MoRE (a Molecular Resistance Estimator), a computational framework that quantifies individual resistance to Aβ and pTau pathology as the difference between observed and expected CSF biomarker levels. Expected values were estimated using iterative XGBoost models incorporating age, sex, education, cognitive status (CN, MCI or AD) and APOE ε4 allele count (0, 1 or 2), given their established associations with AD risk and clinical progression.

Pereira Connectomics Lab, Karolinska Institutet.

## Pipeline

```
 01_A_MORE.R          ->  per-subject resistance residuals
```

## Repository structure

```
.
├── R/
│   ├── utils.R               # argument parsing, small helpers
├── scripts/
│   ├── 01_A_MORE.R
├── data/                     # data not included, see data/README.md
├── results/                  # output (not tracked)
├── install_packages.R
└── LICENSE
```

## Requirements

R ≥ 4.2. Install the dependencies with:

```bash
Rscript install_packages.R
```

The scripts were developed with **xgboost 1.7.x**. xgboost ≥ 2.0 changed its
R interface and is not guaranteed to reproduce the published results.
