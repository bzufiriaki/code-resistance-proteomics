# Proteomic signatures of resistance to Alzheimer's disease pathology

Code for quantifying individual **resistance** to amyloid and tau pathology
(ADNPC: Alzheimer's disease neuropathological change), selecting CSF protein
panels that capture it, and relating those proteins to clinical, imaging,
fluid-biomarker and cognitive trajectories in ADNI.

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
│   ├── xgb_helpers.R         # repeated-CV XGBoost + importance normalisation
│   └── lmm_helpers.R         # mixed-model grid, p-value correction, heatmaps
├── scripts/
│   ├── 01_calculate_resistance.R
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