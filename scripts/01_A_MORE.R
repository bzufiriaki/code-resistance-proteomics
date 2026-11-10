#!/usr/bin/env Rscript
# =============================================================================
# 01_A_MORE.R
#
# Iterative framework that quantifies each individual’s resistance 
# to Aβ and pTau pathology as the gap between observed CSF levels and those 
# expected from well-established AD risk factors, namely age, sex, education, 
# APOE ε4 genotype and clinical status
#
# Resistance = observed pathology - expected pathology
# 
# The expected pathology is predicted out-of-fold by XGBoost from the AD risk 
# factors (age, sex, education, cognitive status and APOE e4 by
# default). The 10-fold CV is repeated 'n_iter' times (default 1000) with
# random generated seeds, so each subject gets one out-of-fold prediction per
# iteration.
#
# Outputs (in <output_dir>/<label>/):
#   residuals_iter_<category>.csv    n_iter x subjects matrix of OOF residuals
#   predictions_iter_<category>.csv  n_iter x subjects matrix of OOF predictions
#                                    (both written with write.csv2: ';' separated)
#   Dataset_<label>_Residual.csv     mean OOF residual per subject
#   Summary_<label>_Residuals.csv    R2 and residual-outcome correlation
#   ImportancePlot_<category>.svg/.tiff
#
# Example (Abeta42 resistance):
#   Rscript scripts/01_A_MORE.R \
#     --adnimerge=data/ADNIMERGE.csv \
#     --dataset=data/individuals.csv \
#     --outcome=ABETA_adj --label=ABETA_adj \
#     --output_dir=results/XGboost_ABETA42
# =============================================================================

suppressPackageStartupMessages({
  library(xgboost)
  library(cvTools)
  library(dplyr)
  library(ggplot2)
  library(svglite)
})

# Locate repo root so helpers can be sourced from anywhere
script_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
                       error = function(e) "scripts")
source(file.path(script_dir, "..", "R", "utils.R"))

# ---- Options -----------------------------------------------------------------
opts <- parse_args(list(
  adnimerge  = NULL,                        # ADNIMERGE csv (APOE4, demographics)
  dataset    = NULL,                        # subject-level dataset with outcome + covariates
  outcome    = "ABETA_adj",                 # column to model (e.g. ABETA_adj, PTAU_adj)
  label      = "ABETA_adj",                 # name used for output folder/files
  predictors = "PTEDUCAT,APOE4,PTGENDER,AGE,Status",
  category   = "State",                     # name of the predictor set
  output_dir = "results/XGboost_ABETA42",
  n_iter     = "1000",
  k_folds    = "10"
))
predictors <- split_csv(opts$predictors)
n_iter     <- as.integer(opts$n_iter)
k_folds    <- as.integer(opts$k_folds)
out_dir    <- ensure_dir(file.path(opts$output_dir, opts$label))

# ---- XGBoost parameters ------------------------------------------------------
xgb_params <- list(
  max_depth        = 3,
  eta              = 0.01,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  lambda           = 10,
  alpha            = 0,
  objective        = "reg:squarederror",
  eval_metric      = "rmse",
  nthread          = -1
)
nrounds               <- 1000
early_stopping_rounds <- 10

# ---- Load and assemble data --------------------------------------------------
adnimerge <- read.csv(require_file(opts$adnimerge))
dataset   <- read.csv(require_file(opts$dataset))

variables_of_interest <- c("RID", "AGE", "PTEDUCAT", "PTGENDER",
               "ABETA_adj", "AV45", "PTAU_adj",
               "FDG", "Ventricles_bl", "Hippocampus_bl", "WholeBrain_bl",
               "Entorhinal_bl", "Fusiform_bl", "MidTemp_bl", "ICV_bl",
               "MMSE", "ADAS13", "RAVLT_perc_forgetting", "MOCA", "FAQ",
               "Status")
variables_of_interest <- union(variables_of_interest, setdiff(c(opts$outcome, predictors), "APOE4"))
dataset <- dataset[, intersect(variables_of_interest, colnames(dataset))]

# Add APOE4 and demographic info from ADNIMERGE (merge() also sorts by RID,
# which determines the fold assignment — keep this order for reproducibility)
dataset <- merge(dataset, unique(adnimerge[, c("APOE4", "RID")]), by = "RID", all.x = TRUE)
dataset <- merge(dataset, unique(adnimerge[, c("PTETHCAT", "RID", "PTRACCAT", "DX_bl")]),
                 by = "RID", all.x = TRUE)

missing_pred <- setdiff(c(predictors, opts$outcome), colnames(dataset))
if (length(missing_pred) > 0) stop("Columns not found in data: ", paste(missing_pred, collapse = ", "))

# ---- Iterative out-of-fold XGBoost ------------------------------------------
#' @return list with mean residuals, R2, correlation and importance
compute_oof_residuals <- function(X_df, y, RID, category_name) {
  # Encode character/factor columns as integer codes
  X <- as.matrix(as.data.frame(lapply(X_df, function(col) {
    if (is.character(col) || is.factor(col)) as.numeric(factor(col)) else as.numeric(col)
  })))

  pred_iter  <- vector("list", n_iter)
  y_iter     <- vector("list", n_iter)
  rid_iter   <- vector("list", n_iter)
  last_model <- NULL

  for (j in seq_len(n_iter)) {
    set.seed(j)
    message(category_name, " iteration: ", j)
    folds <- cvFolds(n = nrow(X), K = k_folds)

    preds <- ys <- rids <- vector("list", k_folds)
    for (i in seq_len(k_folds)) {
      train_idx <- folds$subsets[folds$which != i]
      test_idx  <- folds$subsets[folds$which == i]

      dtrain <- xgb.DMatrix(data = X[train_idx, , drop = FALSE], label = y[train_idx])
      dtest  <- xgb.DMatrix(data = X[test_idx, , drop = FALSE],  label = y[test_idx])

      # Equivalent to the legacy xgboost() wrapper: early stopping monitored
      # on the training set
      model <- xgb.train(params = xgb_params, data = dtrain, nrounds = nrounds,
                         watchlist = list(train = dtrain),
                         early_stopping_rounds = early_stopping_rounds, verbose = 0)

      preds[[i]] <- predict(model, dtest)
      ys[[i]]    <- y[test_idx]
      rids[[i]]  <- RID[test_idx]
      last_model <- model
    }

    rid_j <- unlist(rids)
    ord   <- order(rid_j)
    rid_iter[[j]]  <- rid_j[ord]
    pred_iter[[j]] <- unlist(preds)[ord]
    y_iter[[j]]    <- unlist(ys)[ord]
  }

  resid_iter <- Map(`-`, y_iter, pred_iter)

  # Save per-iteration matrices (rows = iterations, columns = subjects)
  to_matrix <- function(lst) {
    m <- as.data.frame(do.call(rbind, lst))
    colnames(m) <- rid_iter[[1]]
    m
  }
  write.csv2(to_matrix(resid_iter), file.path(out_dir, paste0("residuals_iter_", category_name, ".csv")),   row.names = FALSE)
  write.csv2(to_matrix(pred_iter),  file.path(out_dir, paste0("predictions_iter_", category_name, ".csv")), row.names = FALSE)

  # Summaries across iterations
  mean_resid <- colMeans(do.call(rbind, resid_iter))
  mean_pred  <- colMeans(do.call(rbind, pred_iter))
  y_sorted   <- y_iter[[1]]

  r2      <- 1 - sum((y_sorted - mean_pred)^2) / sum((y_sorted - mean(y_sorted))^2)
  ct      <- cor.test(y_sorted - mean_pred, y_sorted)

  # Feature importance of the last fitted model (last fold, last iteration)
  importance <- xgb.importance(model = last_model)
  p_imp <- ggplot(importance, aes(x = reorder(Feature, Gain), y = Gain)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Feature Importance for", category_name), x = "Feature", y = "Gain") +
    theme_classic()
  ggsave(file.path(out_dir, paste0("ImportancePlot_", category_name, ".svg")),  p_imp, width = 8, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0("ImportancePlot_", category_name, ".tiff")), p_imp, width = 8, height = 6, dpi = 300)

  list(Residuals = mean_resid, RID = rid_iter[[1]], R2 = r2,
       Correlation = unname(ct$estimate), P_value = ct$p.value, Importance = importance)
}

# ---- Run ---------------------------------------------------------------------
message("Analysing category: ", opts$category)
res <- compute_oof_residuals(dataset[, predictors, drop = FALSE],
                             y = dataset[[opts$outcome]],
                             RID = dataset$RID,
                             category_name = opts$category)

# NB: residuals are ordered by RID; dataset is already sorted by RID (merge)
stopifnot(identical(as.numeric(res$RID), as.numeric(dataset$RID)))
residuals_dataset <- data.frame(RID = dataset$RID, PTAU_adj = dataset$PTAU_adj, Status = dataset$Status)
residuals_dataset[[paste0("Residuals_", opts$label, "_", opts$category)]] <- res$Residuals

summary_data <- data.frame(Category = opts$category, R2 = res$R2,
                           Correlation = res$Correlation, P_value = res$P_value)

write.csv(residuals_dataset, file.path(out_dir, paste0("Dataset_", opts$label, "_Residual.csv")),   row.names = FALSE)
write.csv(summary_data,      file.path(out_dir, paste0("Summary_", opts$label, "_Residuals.csv")), row.names = FALSE)
message("Done. Results in ", out_dir)
