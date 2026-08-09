rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================
## PRDX4_RECALC_06B_M070_SEQC_CohortZ_External.R
##
## PURPOSE
##   Technical correction of the SEQC external-prediction preprocessing only.
##
## WHAT IS FROZEN AND MUST NOT CHANGE
##   - M070 model ID: SuperPC -> RSF
##   - the 20 selected genes
##   - the already fitted RSF forest from RECALC 06
##   - risk-score orientation
##   - training risk-score mean/SD
##   - E-MTAB-derived locked risk-group cutoff
##
## WHY 06B EXISTS
##   RECALC 06 correctly locked the model before SEQC, but applying the raw
##   E-MTAB expression center/scale directly to the different SEQC platform
##   produced one identical RSF prediction for all 498 SEQC patients. Thus the
##   external C-index/AUC/HR were not estimable. This is a cross-platform
##   feature-scale problem, not evidence of poor external discrimination.
##
## TECHNICAL CORRECTION
##   SEQC expression is standardized gene-by-gene within the SEQC cohort
##   (mean 0, SD 1) using EXPRESSION VALUES ONLY. No SEQC survival outcome is
##   used for this harmonization. The frozen E-MTAB-trained RSF is then applied
##   without refitting. The original E-MTAB-scaling prediction is also repeated
##   before outcomes are accessed to document the constant-score failure.
##
## IMPORTANT
##   This correction was introduced after observing constant external scores.
##   It must NOT be described as a prespecified preprocessing rule.
## ============================================================

seed_main <- 20260808L
root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"

locked_model_file <- file.path(
  root_dir,
  "SC28_RECALC_06_LockM070_SEQC_External_20260808",
  "Models", "M070_locked_model_BEFORE_SEQC.rds"
)
seqc_ready_file <- file.path(
  root_dir,
  "01_Preprocess", "SEQC_GSE49711_ready_FIX2",
  "SEQC_external_validation_ready_input_FIX2.rds"
)

out_dir <- file.path(
  root_dir,
  "SC28_RECALC_06B_M070_SEQC_CohortZ_External_20260808"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

warning_log <- data.frame(
  Context = character(),
  Message = character(),
  stringsAsFactors = FALSE
)

with_warnlog <- function(context, expr) {
  withCallingHandlers(
    force(expr),
    warning = function(w) {
      warning_log <<- rbind(
        warning_log,
        data.frame(
          Context = as.character(context),
          Message = conditionMessage(w),
          stringsAsFactors = FALSE
        )
      )
      invokeRestart("muffleWarning")
    }
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_event <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[y %in% c(
    "1", "dead", "deceased", "death", "event", "yes", "true", "died",
    "deceased due to disease"
  )] <- 1
  out[y %in% c(
    "0", "alive", "living", "censored", "no", "false",
    "alive without event", "alive with relapse/primary tumor"
  )] <- 0
  z <- suppressWarnings(as.numeric(y))
  use_z <- is.na(out) & !is.na(z)
  out[use_z] <- z[use_z]
  out
}

as_expr_matrix <- function(x) {
  if (is.matrix(x)) {
    mat <- x
  } else if (is.data.frame(x)) {
    df <- x
    first_num <- suppressWarnings(as.numeric(as.character(df[[1]])))
    if (sum(is.na(first_num)) > length(first_num) * 0.5) {
      rownames(df) <- make.unique(as.character(df[[1]]))
      df <- df[, -1, drop = FALSE]
    }
    mat <- as.matrix(as.data.frame(lapply(df, safe_num), check.names = FALSE))
    rownames(mat) <- rownames(df)
  } else {
    stop("Unsupported expression object")
  }
  mode(mat) <- "numeric"
  mat
}

collapse_duplicate_genes <- function(mat) {
  rn <- toupper(trimws(as.character(rownames(mat))))
  keep <- !is.na(rn) & rn != ""
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  rownames(mat) <- rn
  if (!anyDuplicated(rn)) return(mat)
  sp <- split(seq_along(rn), rn)
  z <- do.call(
    rbind,
    lapply(sp, function(i) colMeans(mat[i, , drop = FALSE], na.rm = TRUE))
  )
  rownames(z) <- names(sp)
  z
}

find_time_event_cols <- function(clin) {
  cn <- colnames(clin)
  low <- tolower(cn)
  tp <- c(
    "OS_time", "OS.time", "os_time", "os.time", "OS_time_days",
    "os_time_days", "OVERALLSURVIVAL", "overall_survival", "survival_time"
  )
  ep <- c(
    "OS_event", "OS.event", "os_event", "os.event", "DEATHOFDISEASE",
    "death_of_disease", "overall_survival_event", "vital_status", "event", "status"
  )
  tc_hits <- tp[tp %in% cn]
  ec_hits <- ep[ep %in% cn]
  tc <- if (length(tc_hits)) tc_hits[1] else NA_character_
  ec <- if (length(ec_hits)) ec_hits[1] else NA_character_

  if (is.na(tc)) {
    for (p in c("os.*time", "overall.*survival", "survival.*time", "follow.*up", "time")) {
      h <- cn[grepl(p, low)]
      if (length(h)) {
        for (q in h) {
          if (sum(!is.na(safe_num(clin[[q]]))) >= 30) {
            tc <- q
            break
          }
        }
      }
      if (!is.na(tc)) break
    }
  }

  if (is.na(ec)) {
    for (p in c("os.*event", "death.*disease", "death", "vital.*status", "event", "status")) {
      h <- cn[grepl(p, low)]
      if (length(h)) {
        for (q in h) {
          e <- clean_event(clin[[q]])
          if (sum(!is.na(e)) >= 30 && length(unique(e[!is.na(e)])) >= 2) {
            ec <- q
            break
          }
        }
      }
      if (!is.na(ec)) break
    }
  }
  list(time = tc, event = ec)
}

apply_locked_emtab_scaling <- function(mat, center, scale) {
  genes <- names(center)
  if (!all(genes %in% rownames(mat))) {
    stop(
      "SEQC expression lacks locked candidate gene(s): ",
      paste(setdiff(genes, rownames(mat)), collapse = ", ")
    )
  }
  z <- mat[genes, , drop = FALSE]
  z <- sweep(z, 1, center[genes], "-")
  z <- sweep(z, 1, scale[genes], "/")
  z
}

scale_within_external_cohort <- function(mat, genes) {
  x <- mat[genes, , drop = FALSE]
  mu <- rowMeans(x, na.rm = TRUE)
  sdv <- apply(x, 1, sd, na.rm = TRUE)
  if (any(!is.finite(mu))) {
    stop("Non-finite SEQC expression mean detected")
  }
  bad_sd <- !is.finite(sdv) | sdv <= 0
  if (any(bad_sd)) {
    stop(
      "SEQC gene(s) have zero/non-finite SD: ",
      paste(names(sdv)[bad_sd], collapse = ", ")
    )
  }
  z <- sweep(sweep(x, 1, mu, "-"), 1, sdv, "/")
  if (any(!is.finite(z))) {
    stop("Non-finite values remain after SEQC cohort-wise z standardization")
  }
  list(z = z, center = mu, scale = sdv)
}

predict_locked_rsf <- function(fit, zmat, selected_genes, context) {
  if (!all(selected_genes %in% rownames(zmat))) {
    stop("Prediction matrix lacks selected gene(s) in ", context)
  }
  nd <- as.data.frame(
    t(zmat[selected_genes, , drop = FALSE]),
    check.names = FALSE
  )
  pr <- with_warnlog(context, predict(fit, newdata = nd))
  score <- safe_num(pr$predicted)
  if (length(score) != nrow(nd)) {
    stop("RSF prediction length mismatch in ", context)
  }
  names(score) <- rownames(nd)
  score
}

prediction_diagnostic <- function(method, raw_prediction) {
  x <- raw_prediction[is.finite(raw_prediction)]
  data.frame(
    ScalingMethod = method,
    N = length(raw_prediction),
    FiniteN = length(x),
    UniquePredictionN = length(unique(x)),
    PredictionSD = if (length(x) > 1) sd(x) else NA_real_,
    PredictionMin = if (length(x)) min(x) else NA_real_,
    PredictionMedian = if (length(x)) median(x) else NA_real_,
    PredictionMax = if (length(x)) max(x) else NA_real_,
    ConstantPrediction = length(unique(x)) <= 1,
    stringsAsFactors = FALSE
  )
}

calc_metrics_locked <- function(time, event, score, cutoff, dataset_label) {
  d <- data.frame(
    time = safe_num(time),
    event = safe_num(event),
    score = safe_num(score)
  )
  d <- d[
    complete.cases(d) & d$time > 0 & d$event %in% c(0, 1),
    , drop = FALSE
  ]

  if (nrow(d) < 30 || sum(d$event == 1) < 5 || length(unique(d$score)) < 2) {
    return(data.frame(
      Dataset = dataset_label,
      N = nrow(d), Events = sum(d$event == 1),
      HR_per_trainingSD = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_,
      CoxP = NA_real_, Cindex = NA_real_,
      HighN = NA_integer_, LowN = NA_integer_,
      KM_HR_High_vs_Low = NA_real_, LogRankP = NA_real_,
      AUC1 = NA_real_, AUC3 = NA_real_, AUC5 = NA_real_, MeanAUC = NA_real_
    ))
  }

  cox_fit <- with_warnlog(
    paste0(dataset_label, ": continuous Cox"),
    survival::coxph(survival::Surv(time, event) ~ score, data = d)
  )
  cox_sum <- summary(cox_fit)
  cindex <- tryCatch(
    as.numeric(
      survival::concordance(
        survival::Surv(time, event) ~ score,
        data = d,
        reverse = TRUE
      )$concordance
    ),
    error = function(e) NA_real_
  )

  d$group <- factor(
    ifelse(d$score > cutoff, "High", "Low"),
    levels = c("Low", "High")
  )
  high_n <- sum(d$group == "High")
  low_n <- sum(d$group == "Low")

  if (high_n > 0 && low_n > 0) {
    group_fit <- tryCatch(
      with_warnlog(
        paste0(dataset_label, ": locked-cutoff Cox"),
        survival::coxph(survival::Surv(time, event) ~ group, data = d)
      ),
      error = function(e) NULL
    )
    lr_fit <- tryCatch(
      with_warnlog(
        paste0(dataset_label, ": locked-cutoff log-rank"),
        survival::survdiff(survival::Surv(time, event) ~ group, data = d)
      ),
      error = function(e) NULL
    )
    km_hr <- if (is.null(group_fit)) NA_real_ else exp(as.numeric(coef(group_fit)[1]))
    logrank_p <- if (is.null(lr_fit)) {
      NA_real_
    } else {
      1 - pchisq(lr_fit$chisq, length(lr_fit$n) - 1)
    }
  } else {
    km_hr <- NA_real_
    logrank_p <- NA_real_
  }

  auc <- rep(NA_real_, 3)
  roc_obj <- tryCatch(
    with_warnlog(
      paste0(dataset_label, ": timeROC"),
      timeROC::timeROC(
        T = d$time,
        delta = d$event,
        marker = d$score,
        cause = 1,
        weighting = "marginal",
        times = c(365, 1095, 1825),
        iid = FALSE
      )
    ),
    error = function(e) NULL
  )
  if (!is.null(roc_obj)) auc <- as.numeric(roc_obj$AUC)

  data.frame(
    Dataset = dataset_label,
    N = nrow(d),
    Events = sum(d$event == 1),
    HR_per_trainingSD = cox_sum$coefficients[1, "exp(coef)"],
    Lower95 = cox_sum$conf.int[1, "lower .95"],
    Upper95 = cox_sum$conf.int[1, "upper .95"],
    CoxP = cox_sum$coefficients[1, "Pr(>|z|)"],
    Cindex = cindex,
    HighN = high_n,
    LowN = low_n,
    KM_HR_High_vs_Low = km_hr,
    LogRankP = logrank_p,
    AUC1 = auc[1],
    AUC3 = auc[2],
    AUC5 = auc[3],
    MeanAUC = if (all(is.na(auc))) NA_real_ else mean(auc, na.rm = TRUE)
  )
}

main <- function() {
  on.exit({
    write.csv(
      warning_log,
      file.path(out_dir, "Logs", "01_warnings_captured.csv"),
      row.names = FALSE
    )
    writeLines(
      capture.output(sessionInfo()),
      file.path(out_dir, "Logs", "02_sessionInfo.txt")
    )
  }, add = TRUE)

  ## ---------------------------- package gate ----------------------------
  required_packages <- c("survival", "randomForestSRC", "timeROC")
  package_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  write.csv(
    data.frame(Package = required_packages, Available = package_ok),
    file.path(out_dir, "Tables", "00_package_status.csv"),
    row.names = FALSE
  )
  if (any(!package_ok)) {
    missing_packages <- required_packages[!package_ok]
    stop(
      "Missing required package(s): ", paste(missing_packages, collapse = ", "),
      ". Install once with: install.packages(c(\"",
      paste(missing_packages, collapse = "\",\""),
      "\"), repos=\"https://cloud.r-project.org\") ; then restart R and rerun this WHOLE script."
    )
  }

  if (!file.exists(locked_model_file)) {
    stop("Locked RECALC 06 model not found: ", locked_model_file)
  }
  if (!file.exists(seqc_ready_file)) {
    stop("SEQC ready input not found: ", seqc_ready_file)
  }

  ## ---------------------------- load FROZEN model first ----------------------------
  ## This file was created by RECALC 06 BEFORE SEQC was opened.
  locked <- readRDS(locked_model_file)
  required_fields <- c(
    "ModelID", "Selector", "Learner", "CandidateGenes", "SelectedGenes",
    "ExpressionCenter", "ExpressionScale", "RSFModel",
    "RiskOrientationSign", "TrainingRiskMean", "TrainingRiskSD",
    "LockedRiskCutoff", "SEQCUsedForModelSelection", "SEQCUsedForCutoff"
  )
  missing_fields <- setdiff(required_fields, names(locked))
  if (length(missing_fields)) {
    stop("Locked model object lacks field(s): ", paste(missing_fields, collapse = ", "))
  }
  if (!identical(as.character(locked$ModelID), "M070_SuperPC_RSF")) {
    stop("Unexpected locked model ID: ", as.character(locked$ModelID))
  }
  if (!identical(as.character(locked$Selector), "SuperPC") ||
      !identical(as.character(locked$Learner), "RSF")) {
    stop("Locked pipeline is not SuperPC -> RSF")
  }
  candidate_genes <- as.character(locked$CandidateGenes)
  selected_genes <- as.character(locked$SelectedGenes)
  if (length(candidate_genes) != 24L || anyDuplicated(candidate_genes)) {
    stop("Locked candidate panel is not the expected unique 24-gene panel")
  }
  if (length(selected_genes) != 20L || anyDuplicated(selected_genes)) {
    stop("Locked selected set is not the expected unique 20-gene set")
  }
  if (!all(selected_genes %in% candidate_genes)) {
    stop("Locked selected genes are not a subset of the candidate panel")
  }
  if (!identical(locked$SEQCUsedForModelSelection, FALSE) ||
      !identical(locked$SEQCUsedForCutoff, FALSE)) {
    stop("Locked model provenance indicates SEQC leakage")
  }

  orientation_sign <- safe_num(locked$RiskOrientationSign)[1]
  train_risk_mean <- safe_num(locked$TrainingRiskMean)[1]
  train_risk_sd <- safe_num(locked$TrainingRiskSD)[1]
  locked_cutoff <- safe_num(locked$LockedRiskCutoff)[1]
  if (!all(is.finite(c(orientation_sign, train_risk_mean, train_risk_sd, locked_cutoff))) ||
      train_risk_sd <= 0 || !orientation_sign %in% c(-1, 1)) {
    stop("Invalid locked risk-score transformation parameters")
  }

  model_md5 <- unname(tools::md5sum(locked_model_file))
  model_copy_file <- file.path(
    out_dir, "Models", "M070_locked_model_USED_UNCHANGED.rds"
  )
  copied <- file.copy(locked_model_file, model_copy_file, overwrite = TRUE)
  if (!isTRUE(copied)) stop("Could not copy frozen M070 model into 06B output")
  copied_md5 <- unname(tools::md5sum(model_copy_file))
  if (!identical(model_md5, copied_md5)) {
    stop("Frozen-model MD5 mismatch after copy")
  }

  model_integrity <- data.frame(
    Item = c(
      "Model ID", "Selector", "Learner", "Candidate gene N", "Selected gene N",
      "Selected genes", "PRDX4 selected", "Frozen model MD5",
      "Risk orientation sign", "Training risk mean", "Training risk SD",
      "Locked cutoff", "Model refit in 06B", "Genes reselected in 06B",
      "Cutoff changed in 06B", "SEQC outcomes used for model selection",
      "SEQC outcomes used for harmonization"
    ),
    Value = c(
      locked$ModelID, locked$Selector, locked$Learner,
      length(candidate_genes), length(selected_genes),
      paste(selected_genes, collapse = ";"), "PRDX4" %in% selected_genes,
      model_md5, orientation_sign, train_risk_mean, train_risk_sd, locked_cutoff,
      "FALSE", "FALSE", "FALSE", "FALSE", "FALSE"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    model_integrity,
    file.path(out_dir, "Tables", "01_locked_model_integrity_BEFORE_SEQC.csv"),
    row.names = FALSE
  )

  method_note <- data.frame(
    Item = c(
      "Why 06B was required",
      "Original 06 external preprocessing",
      "Original 06 failure mode",
      "06B technical correction",
      "Outcome information used for correction",
      "Model refitted",
      "Feature set changed",
      "Risk cutoff changed",
      "How this correction should be reported"
    ),
    Value = c(
      "Original RECALC 06 produced a constant RSF external prediction in SEQC",
      "E-MTAB expression center/scale applied directly to SEQC",
      "All SEQC patients received the same RSF risk prediction; external performance was therefore not estimable",
      "Gene-wise z standardization within SEQC using expression values only; apply the already frozen M070 RSF unchanged",
      "None",
      "No",
      "No",
      "No; E-MTAB-derived cutoff remains locked",
      "Post-audit cross-platform preprocessing correction; do not describe as prespecified"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    method_note,
    file.path(out_dir, "Tables", "02_method_correction_note.csv"),
    row.names = FALSE
  )

  cat("Frozen model loaded and verified BEFORE SEQC.\n")
  cat("Model MD5: ", model_md5, "\n", sep = "")
  cat("Selected genes: ", paste(selected_genes, collapse = ", "), "\n", sep = "")
  cat("PRDX4 selected: ", "PRDX4" %in% selected_genes, "\n", sep = "")

  ## ---------------------------- open SEQC expression ----------------------------
  ## The model is already frozen. Survival outcomes are NOT accessed yet.
  seqc_obj <- readRDS(seqc_ready_file)
  if (is.null(seqc_obj$expr)) {
    stop("SEQC ready object must contain $expr")
  }
  seqc_expr <- collapse_duplicate_genes(as_expr_matrix(seqc_obj$expr))
  missing_seqc <- setdiff(candidate_genes, rownames(seqc_expr))
  if (length(missing_seqc)) {
    stop(
      "Frozen candidate gene(s) missing from SEQC: ",
      paste(missing_seqc, collapse = ", ")
    )
  }

  emtab_center <- safe_num(locked$ExpressionCenter)
  emtab_scale <- safe_num(locked$ExpressionScale)
  names(emtab_center) <- names(locked$ExpressionCenter)
  names(emtab_scale) <- names(locked$ExpressionScale)
  if (!all(candidate_genes %in% names(emtab_center)) ||
      !all(candidate_genes %in% names(emtab_scale))) {
    stop("Locked E-MTAB scaling parameters do not cover all 24 candidate genes")
  }

  ## Reproduce original RECALC 06 scaling WITHOUT outcomes.
  seqc_emtab_z <- apply_locked_emtab_scaling(
    seqc_expr,
    center = emtab_center[candidate_genes],
    scale = emtab_scale[candidate_genes]
  )
  original_raw_prediction <- predict_locked_rsf(
    locked$RSFModel,
    seqc_emtab_z,
    selected_genes,
    "Original 06 E-MTAB-scaling prediction reproduced before outcomes"
  )

  ## Corrected expression-only cross-platform harmonization.
  seqc_cohort_scale <- scale_within_external_cohort(seqc_expr, candidate_genes)
  seqc_cohort_z <- seqc_cohort_scale$z
  corrected_raw_prediction <- predict_locked_rsf(
    locked$RSFModel,
    seqc_cohort_z,
    selected_genes,
    "06B SEQC cohort-z prediction before outcomes"
  )

  if (!identical(names(original_raw_prediction), names(corrected_raw_prediction))) {
    stop("Sample order mismatch between original and corrected SEQC predictions")
  }

  scaling_diag <- data.frame(
    Gene = candidate_genes,
    E_MTAB_RawCenter = as.numeric(emtab_center[candidate_genes]),
    E_MTAB_RawSD = as.numeric(emtab_scale[candidate_genes]),
    SEQC_RawMean = as.numeric(seqc_cohort_scale$center[candidate_genes]),
    SEQC_RawSD = as.numeric(seqc_cohort_scale$scale[candidate_genes]),
    SEQC_MeanShift_in_EMTAB_SD = as.numeric(
      (seqc_cohort_scale$center[candidate_genes] - emtab_center[candidate_genes]) /
        emtab_scale[candidate_genes]
    ),
    SEQC_SD_to_EMTAB_SD_Ratio = as.numeric(
      seqc_cohort_scale$scale[candidate_genes] / emtab_scale[candidate_genes]
    ),
    SEQC_CohortZ_Mean = rowMeans(seqc_cohort_z[candidate_genes, , drop = FALSE]),
    SEQC_CohortZ_SD = apply(seqc_cohort_z[candidate_genes, , drop = FALSE], 1, sd),
    SelectedInM070 = candidate_genes %in% selected_genes,
    stringsAsFactors = FALSE
  )
  write.csv(
    scaling_diag,
    file.path(out_dir, "Tables", "03_cross_platform_scaling_diagnostics_BEFORE_OUTCOMES.csv"),
    row.names = FALSE
  )

  prediction_diag <- rbind(
    prediction_diagnostic(
      "Original06_apply_EMTAB_centerSD_to_SEQC",
      original_raw_prediction
    ),
    prediction_diagnostic(
      "Corrected06B_SEQC_cohort_gene_zscore",
      corrected_raw_prediction
    )
  )
  write.csv(
    prediction_diag,
    file.path(out_dir, "Tables", "04_prediction_variation_diagnostics_BEFORE_OUTCOMES.csv"),
    row.names = FALSE
  )

  corrected_oriented <- orientation_sign * corrected_raw_prediction
  corrected_score <- (corrected_oriented - train_risk_mean) / train_risk_sd
  corrected_group <- ifelse(corrected_score > locked_cutoff, "High", "Low")

  preoutcome_scores <- data.frame(
    Sample = names(corrected_raw_prediction),
    Original06_RawPrediction = as.numeric(original_raw_prediction),
    Corrected06B_RawPrediction = as.numeric(corrected_raw_prediction),
    Corrected06B_M070_RiskScore = as.numeric(corrected_score),
    RiskGroup_Locked_EMTAB_Cutoff = corrected_group,
    stringsAsFactors = FALSE
  )
  write.csv(
    preoutcome_scores,
    file.path(out_dir, "Tables", "05_SEQC_predictions_BEFORE_OUTCOMES.csv"),
    row.names = FALSE
  )

  corrected_unique_n <- length(unique(corrected_raw_prediction[is.finite(corrected_raw_prediction)]))
  corrected_sd <- sd(corrected_raw_prediction, na.rm = TRUE)
  if (corrected_unique_n <= 1 || !is.finite(corrected_sd) || corrected_sd <= 0) {
    stop(
      "06B cohort-z prediction is still constant. STOPPED BEFORE SEQC outcomes. ",
      "Upload the 06B output folder for diagnosis."
    )
  }

  cat("\nPrediction variation restored BEFORE survival outcomes are accessed.\n")
  print(prediction_diag)
  cat("Locked-cutoff groups before outcomes: High=",
      sum(corrected_group == "High"),
      " Low=", sum(corrected_group == "Low"), "\n", sep = "")

  ## ---------------------------- NOW access SEQC outcomes ----------------------------
  if (is.null(seqc_obj$clin)) {
    stop("SEQC ready object must contain $clin")
  }
  seqc_clin <- as.data.frame(seqc_obj$clin)
  if (!"Sample" %in% colnames(seqc_clin)) seqc_clin$Sample <- rownames(seqc_clin)
  seqc_clin$Sample <- as.character(seqc_clin$Sample)
  te <- find_time_event_cols(seqc_clin)
  if (is.na(te$time) || is.na(te$event)) {
    stop("Cannot identify SEQC survival time/event columns")
  }
  seqc_time_all <- safe_num(seqc_clin[[te$time]])
  seqc_event_all <- clean_event(seqc_clin[[te$event]])

  result <- preoutcome_scores
  idx <- match(result$Sample, seqc_clin$Sample)
  result$OS_time <- seqc_time_all[idx]
  result$OS_event <- seqc_event_all[idx]
  result <- result[
    complete.cases(result[, c("OS_time", "OS_event", "Corrected06B_M070_RiskScore")]) &
      result$OS_time > 0 & result$OS_event %in% c(0, 1),
    , drop = FALSE
  ]
  if (nrow(result) < 30) stop("Too few SEQC patients after survival-data matching")
  result <- result[, c(
    "Sample", "OS_time", "OS_event",
    "Corrected06B_RawPrediction", "Corrected06B_M070_RiskScore",
    "RiskGroup_Locked_EMTAB_Cutoff"
  )]
  write.csv(
    result,
    file.path(out_dir, "Tables", "06_SEQC_EXTERNAL_M070_scores_with_outcomes.csv"),
    row.names = FALSE
  )

  perf <- calc_metrics_locked(
    result$OS_time,
    result$OS_event,
    result$Corrected06B_M070_RiskScore,
    cutoff = locked_cutoff,
    dataset_label = "SEQC external validation; frozen M070; cohort-wise expression z harmonization"
  )
  write.csv(
    perf,
    file.path(out_dir, "Tables", "07_SEQC_EXTERNAL_M070_performance.csv"),
    row.names = FALSE
  )

  final_summary <- data.frame(
    Item = c(
      "Frozen model", "Model refitted in 06B", "Genes reselected in 06B",
      "PRDX4 selected", "Frozen model MD5",
      "Original06 unique SEQC predictions",
      "Corrected06B unique SEQC predictions",
      "Corrected06B SEQC prediction SD",
      "SEQC N", "SEQC events",
      "Locked cutoff source", "Locked cutoff",
      "SEQC High N by locked cutoff", "SEQC Low N by locked cutoff",
      "SEQC survival time column", "SEQC survival event column",
      "SEQC outcomes used for harmonization", "Warnings captured"
    ),
    Value = c(
      locked$ModelID, "FALSE", "FALSE", "PRDX4" %in% selected_genes,
      model_md5,
      prediction_diag$UniquePredictionN[1], prediction_diag$UniquePredictionN[2],
      prediction_diag$PredictionSD[2],
      nrow(result), sum(result$OS_event == 1),
      "E-MTAB full-refit median; inherited unchanged from RECALC 06",
      locked_cutoff,
      sum(result$RiskGroup_Locked_EMTAB_Cutoff == "High"),
      sum(result$RiskGroup_Locked_EMTAB_Cutoff == "Low"),
      te$time, te$event, "FALSE", nrow(warning_log)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    final_summary,
    file.path(out_dir, "Tables", "08_FINAL_06B_SUMMARY.csv"),
    row.names = FALSE
  )

  summary_lines <- c(
    "SC28 RECALC 06B completed successfully.",
    paste0("End: ", as.character(Sys.time())),
    "Frozen model used unchanged: M070 = SuperPC -> RSF",
    "Model refit: FALSE",
    "Gene reselection: FALSE",
    paste0("PRDX4 selected in frozen model: ", "PRDX4" %in% selected_genes),
    paste0("Frozen model MD5: ", model_md5),
    paste0("Original06 unique SEQC predictions: ", prediction_diag$UniquePredictionN[1]),
    paste0("Corrected06B unique SEQC predictions: ", prediction_diag$UniquePredictionN[2]),
    paste0("SEQC N/events: ", nrow(result), "/", sum(result$OS_event == 1)),
    paste0(
      "Locked-cutoff High/Low: ",
      sum(result$RiskGroup_Locked_EMTAB_Cutoff == "High"), "/",
      sum(result$RiskGroup_Locked_EMTAB_Cutoff == "Low")
    ),
    paste0("Warnings captured: ", nrow(warning_log)),
    "SEQC outcomes were not used for cross-platform harmonization.",
    "Primary external-validation results: Tables/07_SEQC_EXTERNAL_M070_performance.csv"
  )
  writeLines(summary_lines, file.path(out_dir, "Logs", "03_run_summary.txt"))

  cat("\n============================================================\n")
  cat(paste(summary_lines, collapse = "\n"), "\n")
  cat("Output folder:\n", out_dir, "\n", sep = "")
  cat("============================================================\n")
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(
      c(
        paste0("Time: ", as.character(Sys.time())),
        paste0("ERROR: ", conditionMessage(e)),
        "If diagnostic tables were created, keep them and upload the entire 06B output folder."
      ),
      file.path(out_dir, "Logs", "99_ERROR.txt")
    )
    stop(e)
  }
)
