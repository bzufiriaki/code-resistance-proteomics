# -----------------------------------------------------------------------------
# lmm_helpers.R — linear mixed models of protein level vs. clinical outcomes
# used by scripts/03_clinical_analyses.R
# -----------------------------------------------------------------------------

# Outcome variables grouped by modality --------------------------------------
outcome_sets <- list(
  CSF       = c("ABETA_adj_log_z", "PTAU_adj_log_z", "TAU_log_z"),
  PET       = c("AV45_z", "FDG_z"),
  MRI       = c("Ventricles_z", "Hippocampus_z", "WholeBrain_z",
                "Entorhinal_z", "Fusiform_z", "MidTemp_z"),
  `PET Tau` = c("Braak_I_SUVR_z", "Braak_III_IV_SUVR_z",
                "Braak_V_VI_SUVR_z", "META_TEMPORAL_SUVR_z"),
  Plasma    = c("Plasma_pTau217_log_z", "Plasma_AB42_log_z", "Plasma_AB40_log_z",
                "Plasma_AB42_AB40_Ratio_log_z", "Plasma_pTau217_AB42_Ratio_log_z",
                "Plasma_NfL_log_z", "Plasma_GFAP_log_z"),
  Cognitive = c("ADNI_EF", "ADNI_EF2", "ADNI_MEM", "PHC_MEM", "PHC_EXF",
                "PHC_LAN", "PHC_VSP", "MMSCORE_z", "PACC", "mPACCdigit", "mPACCtrailsB")
)
# Order in which outcomes are modelled (kept identical to the original script)
all_outcomes <- unlist(outcome_sets[c("CSF", "PET", "MRI", "PET Tau", "Plasma", "Cognitive")],
                       use.names = FALSE)
category_levels <- c("CSF", "Plasma", "PET", "PET Tau", "MRI", "Cognitive")

# Transformed variable name -> original variable name
transform_map <- c(
  ABETA_adj_log_z = "ABETA_adj", PTAU_adj_log_z = "PTAU_adj", TAU_log_z = "TAU",
  AV45_z = "AV45", FDG_z = "FDG",
  Braak_I_SUVR_z = "Braak_I_SUVR", Braak_III_IV_SUVR_z = "Braak_III_IV_SUVR",
  Braak_V_VI_SUVR_z = "Braak_V_VI_SUVR", META_TEMPORAL_SUVR_z = "META_TEMPORAL_SUVR",
  Ventricles_z = "Ventricles", Hippocampus_z = "Hippocampus", WholeBrain_z = "WholeBrain",
  Entorhinal_z = "Entorhinal", Fusiform_z = "Fusiform", MidTemp_z = "MidTemp",
  Plasma_pTau217_log_z = "Plasma_pTau217", Plasma_AB42_log_z = "Plasma_AB42",
  Plasma_AB40_log_z = "Plasma_AB40", Plasma_AB42_AB40_Ratio_log_z = "Plasma_AB42_AB40_Ratio",
  Plasma_pTau217_AB42_Ratio_log_z = "Plasma_pTau217_AB42_Ratio",
  Plasma_NfL_log_z = "Plasma_NfL", Plasma_GFAP_log_z = "Plasma_GFAP",
  ADNI_EF = "ADNI_EF", ADNI_EF2 = "ADNI_EF2", ADNI_MEM = "ADNI_MEM",
  PHC_MEM = "PHC_MEM", PHC_EXF = "PHC_EXF", PHC_LAN = "PHC_LAN", PHC_VSP = "PHC_VSP",
  MOCA = "MOCA", MMSCORE_z = "MMSCORE_z", PACC = "PACC",
  mPACCdigit = "mPACCdigit", mPACCtrailsB = "mPACCtrailsB"
)

# Original variable name -> label used in plots (order = plot order)
var_labels <- c(
  ABETA_adj = "Aβ42", PTAU_adj = "pTau181", TAU = "tTau",
  Plasma_AB40 = "Aβ40", Plasma_AB42 = " Aβ42", Plasma_AB42_AB40_Ratio = "Aβ42/40",
  Plasma_GFAP = "GFAP", Plasma_NfL = "NfL", Plasma_pTau217 = "pTau217",
  Plasma_pTau217_AB42_Ratio = "pTau217/Aβ42",
  AV45 = "AV45 (Amyloid)", FDG = "FDG (Glucose)",
  Braak_I_SUVR = "Braak I", Braak_III_IV_SUVR = "Braak III–IV",
  Braak_V_VI_SUVR = "Braak V–VI", META_TEMPORAL_SUVR = "Meta-Temporal",
  Entorhinal = "Entorhinal", Fusiform = "Fusiform", Hippocampus = "Hippocampus",
  MidTemp = "Mid Temporal", Ventricles = "Ventricles", WholeBrain = "Whole Brain",
  ADNI_EF = "Exec. Function", ADNI_EF2 = "Exec. Function 2", ADNI_MEM = "Memory",
  PHC_MEM = "Memory (PHC)", PHC_EXF = "Exec. Function (PHC)", PHC_LAN = "Language (PHC)",
  PHC_VSP = "Visuospatial (PHC)", MOCA = "MoCA", MMSCORE_z = "MMSE", PACC = "PACC",
  mPACCdigit = "mPACC Digit", mPACCtrailsB = "mPACC Trails-B"
)

outcome_category <- function(v) {
  out <- rep("Other", length(v))
  for (cat_name in names(outcome_sets)) out[v %in% outcome_sets[[cat_name]]] <- cat_name
  out
}

#' Fit one LMM per (group x protein x outcome)
#'
#' Model: outcome ~ Protein_Level * Years + APOE4 + sex + age + education
#'                  [+ extra covariate] [+ scale(ICV) for MRI] + (1 | RID)
#' Protein_Level is z-scored within each group.
#'
#' @param data long-format data (one row per visit)
#' @param group_col column defining the groups analysed separately
#' @param proteins protein column names
#' @param extra_covariate optional additional covariate (e.g. baseline biomarker)
fit_lmm_grid <- function(data, group_col, proteins, extra_covariate = NULL) {
  all_results <- list()

  for (grp in unique(data[[group_col]])) {
    message("Group: ", grp)
    grp_data <- data[data[[group_col]] == grp, , drop = FALSE]

    for (protein in proteins) {
      for (var in all_outcomes) {
        tmp <- grp_data
        tmp$Outcome_Var   <- tmp[[var]]
        tmp$Protein_Level <- scale(tmp[[protein]])[, 1]
        if (sum(!is.na(tmp$Outcome_Var)) == 0) next

        rhs <- paste(c("Protein_Level * Years", "as.factor(APOE4)", "as.factor(PTGENDER)",
                       "AGE", "PTEDUCAT", extra_covariate,
                       if (var %in% outcome_sets$MRI) "scale(ICV)"),
                     collapse = " + ")
        fml <- as.formula(paste("as.numeric(Outcome_Var) ~", rhs, "+ (1 | RID)"))

        fit <- tryCatch(lmerTest::lmer(fml, data = tmp), error = function(e) {
          message("  Skipping ", protein, " / ", var, ": ", conditionMessage(e))
          NULL
        })
        if (is.null(fit)) next

        coefs <- as.data.frame(coef(summary(fit)))
        all_results[[length(all_results) + 1]] <- coefs %>%
          tibble::rownames_to_column("term") %>%
          dplyr::rename(estimate = Estimate, std.error = `Std. Error`,
                        statistic = `t value`, p.value = `Pr(>|t|)`) %>%
          dplyr::mutate(effect = "fixed", Outcome_Variable = var,
                        Protein = protein, Status = grp)
      }
    }
  }
  dplyr::bind_rows(all_results)
}

signif_stars <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

#' Keep the protein main effect and protein x time interaction and correct p-values
#'
#' p_adj  : Bonferroni across all retained tests
#' p_bonf : Bonferroni within the groups given in `bonf_by`
postprocess_lmm <- function(all_results, bonf_by = "Category") {
  all_results %>%
    dplyr::filter(term %in% c("Protein_Level", "Protein_Level:Years")) %>%
    dplyr::mutate(
      Type          = ifelse(term == "Protein_Level", "Cross-sectional", "Longitudinal"),
      Base_Var      = transform_map[Outcome_Variable],
      Category      = outcome_category(Outcome_Variable),
      Outcome_Label = factor(var_labels[Base_Var], levels = var_labels),
      p_adj         = p.adjust(p.value, method = "bonferroni"),
      Significance  = signif_stars(p_adj)
    ) %>%
    dplyr::filter(!is.na(Category) & !is.na(Outcome_Label)) %>%
    dplyr::mutate(Category = factor(Category, levels = category_levels)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(bonf_by))) %>%
    dplyr::mutate(p_bonf = p.adjust(p.value, method = "bonferroni")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Significance_bonf = signif_stars(p_bonf))
}

#' Heatmap of protein effects (one PDF per group)
plot_lmm_heatmaps <- function(results, out_dir) {
  for (st in unique(results$Status)) {
    tmp <- results %>%
      dplyr::filter(Status == st) %>%
      dplyr::mutate(alpha_val = ifelse(p_bonf <= 0.05, -log10(p_bonf), 0))
    tmp_cross <- dplyr::filter(tmp, Type == "Cross-sectional")
    tmp_long  <- dplyr::filter(tmp, Type == "Longitudinal")

    p <- ggplot() +
      geom_tile(data = tmp_cross, aes(x = Protein, y = Outcome_Label, fill = estimate, alpha = alpha_val),
                color = "white") +
      scale_fill_gradient2(low = "#00AFBB", mid = "white", high = "#bb0c00",
                           limits = c(-1, 1), midpoint = 0, name = "β") +
      ggnewscale::new_scale_fill() +
      geom_tile(data = tmp_long, aes(x = Protein, y = Outcome_Label, fill = estimate, alpha = alpha_val),
                color = "white") +
      scale_fill_gradient2(low = "#00AFBB", mid = "white", high = "#bb0c00",
                           limits = c(-0.2, 0.2), midpoint = 0, name = "β*time") +
      geom_text(data = tmp, aes(x = Protein, y = Outcome_Label, label = Significance_bonf), size = 4) +
      facet_grid(Category ~ Type, scales = "free_y", space = "free_y") +
      scale_alpha(range = c(0.1, 1), limits = c(0, 1.301), name = NULL) +
      labs(title = paste("Protein effects on outcomes in:", st), x = "", y = "") +
      theme_classic() +
      guides(alpha = "none") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y  = element_text(size = 9),
            strip.text.y = element_text(face = "bold"),
            plot.title   = element_text(face = "bold", hjust = 0.5))

    ggsave(file.path(out_dir, paste0("Plot_", st, ".pdf")), plot = p,
           width = 10, height = 8, dpi = 300, device = cairo_pdf)
  }
}

#' Fit, post-process, save CSV and heatmaps for one analysis
run_lmm_analysis <- function(data, group_col, proteins, out_dir, csv_name,
                             extra_covariate = NULL, bonf_by = "Category") {
  ensure_dir(out_dir)
  res <- fit_lmm_grid(data, group_col, proteins, extra_covariate) %>%
    postprocess_lmm(bonf_by = bonf_by)
  write.csv(res, file.path(out_dir, csv_name), row.names = FALSE)
  plot_lmm_heatmaps(res, out_dir)
  invisible(res)
}
