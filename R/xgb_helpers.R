# -----------------------------------------------------------------------------
# xgb_helpers.R — XGBoost cross-validation and feature-importance helpers
# used by scripts/02_proteomic_biomarker_selection.R
# -----------------------------------------------------------------------------

#' Default XGBoost parameters for the protein-selection steps
xgb_params_selection <- function() {
  list(
    objective        = "reg:squarederror",
    eval_metric      = "rmse",
    booster          = "gbtree",
    max_depth        = 6,
    eta              = 0.01,
    subsample        = 0.7,
    colsample_bytree = 0.7,
    nthread          = 1
  )
}

#' One iteration of 10-fold cross-validated XGBoost
#'
#' The seed is fixed to 123 + iteration so that every iteration is
#' reproducible independently of the parallel backend.
#'
#' @return list with per-fold Gain / Cover importance and test RMSE, R2, r.
run_xgb_cv <- function(iteration, X, y, params, nrounds = 100, k = 10) {
  set.seed(123 + iteration)
  folds <- caret::createFolds(y, k = k, list = TRUE, returnTrain = TRUE)

  gain <- cover <- rmse <- r2 <- r <- vector("list", length(folds))

  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(X)), train_idx)

    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X[train_idx, , drop = FALSE]), label = y[train_idx])
    dtest  <- xgboost::xgb.DMatrix(data = as.matrix(X[test_idx, , drop = FALSE]),  label = y[test_idx])

    model <- xgboost::xgb.train(params, dtrain, nrounds = nrounds, verbose = 0)

    # Feature importance (training fold)
    imp <- xgboost::xgb.importance(feature_names = colnames(X), model = model)
    gain[[i]]  <- setNames(imp$Gain,  imp$Feature)
    cover[[i]] <- setNames(imp$Cover, imp$Feature)

    # Performance (held-out fold)
    y_test <- y[test_idx]
    pred   <- predict(model, dtest)
    rmse[[i]] <- sqrt(mean((pred - y_test)^2))
    r2[[i]]   <- 1 - sum((y_test - pred)^2) / sum((y_test - mean(y_test))^2)
    r[[i]]    <- unname(cor.test(pred, y_test)$estimate)
  }

  list(Gain = gain, Cover = cover, RMSE = rmse, R2 = r2, R = r)
}

#' Run `n_iter` iterations of run_xgb_cv (in parallel via future.apply)
run_xgb_iterations <- function(X, y, params, n_iter = 100) {
  future.apply::future_lapply(
    seq_len(n_iter),
    function(it) run_xgb_cv(it, X, y, params),
    future.seed = TRUE
  )
}

#' Normalise importance to 100 within each (iteration, fold) and average
#'
#' Proteins that are not used in a given fold are absent (not zero) for that
#' fold, so the mean is taken over the folds in which the protein appears.
#'
#' @param imp_list list over iterations of lists over folds of named vectors.
#' @return data.frame with Protein, Mean_Value, Standard_Deviation (sorted).
normalize_importance_per_fold <- function(imp_list) {
  rows <- list()
  for (it in seq_along(imp_list)) {
    for (fold in seq_along(imp_list[[it]])) {
      v <- imp_list[[it]][[fold]]
      if (!is.null(v) && length(v) > 0) {
        rows[[length(rows) + 1]] <- data.frame(
          Iteration = it, Fold = fold, Protein = names(v), Value = as.numeric(v)
        )
      }
    }
  }
  if (length(rows) == 0) stop("No importance values found: check the input structure.")

  dplyr::bind_rows(rows) %>%
    dplyr::group_by(Iteration, Fold) %>%
    dplyr::mutate(Normalized_Value = 100 * Value / sum(Value, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(
      Mean_Value         = mean(Normalized_Value, na.rm = TRUE),
      Standard_Deviation = sd(Normalized_Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Mean_Value))
}

#' Summarise iterations: write performance tables and return merged importance
#'
#' @param results output of run_xgb_iterations()
#' @param out_prefix path prefix for the RMSE / R2 / r CSV files
#' @return data.frame with Gain, Cover and Ensemble (= mean of both) importance
summarise_xgb_iterations <- function(results, out_dir, step_label) {
  collect <- function(field) do.call(rbind, lapply(results, function(res) unlist(res[[field]])))
  write.csv(collect("RMSE"), file.path(out_dir, paste0("XGBoost1_rmse_", step_label, ".csv")),      row.names = FALSE)
  write.csv(collect("R2"),   file.path(out_dir, paste0("XGBoost1_r_squared_", step_label, ".csv")), row.names = FALSE)
  write.csv(collect("R"),    file.path(out_dir, paste0("XGBoost1_r_", step_label, ".csv")),         row.names = FALSE)

  gain_df  <- normalize_importance_per_fold(lapply(results, `[[`, "Gain"))
  cover_df <- normalize_importance_per_fold(lapply(results, `[[`, "Cover"))

  gain_df %>%
    dplyr::full_join(cover_df, by = "Protein", suffix = c("_Gain", "_Cover")) %>%
    dplyr::mutate(Ensemble = (Mean_Value_Gain + Mean_Value_Cover) / 2)
}
