# Proteomic signatures of resistance to Alzheimer's disease pathology

Code for quantifying individual **resistance** to amyloid and tau pathology
(ADNPC: Alzheimer's disease neuropathological change).

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
├── data/                     # ADNI data (not included, see data/README.md)
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
