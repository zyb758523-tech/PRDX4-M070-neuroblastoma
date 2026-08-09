rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================================
## PRDX4_RECALC_07A_M070_Prognostic_ClinicalIncremental.R
##
## Purpose
##   Final prognostic and clinical-incremental-value refresh after locking M070.
##
## Frozen design
##   1. Model selection evidence comes from SC28 RECALC 05 repeated CV in
##      E-MTAB-179 only.
##   2. The final model is M070 = SuperPC -> RSF, frozen in RECALC 06 before
##      SEQC outcomes were accessed.
##   3. SEQC uses the audited 06B cohort-wise gene z harmonization because the
##      two expression platforms are on incompatible raw scales. No model
##      fitting, feature selection, cutoff tuning, or outcome-guided scaling is
##      performed in SEQC.
##   4. Primary displayed risk groups use the E-MTAB-179 locked cutoff.
##   5. Clinical incremental analyses standardize the already-fixed M070 score
##      within each cohort (M070_z) only for interpretable Cox coefficients.
##      This affine transformation does not refit the prognostic model.
##
## Clinical covariates (kept identical to the earlier curated analysis)
##   E-MTAB-179: Age_months + Male + INSS_stage4
##   SEQC:       Age_months + Male
##
## Expected complete-case counts from the audited prior clinical merge
##   E-MTAB-179: N = 389, events = 57
##   SEQC:       N = 498, events = 105
## ============================================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
out_dir <- file.path(root_dir, "SC28_RECALC_07A_M070_Prognostic_ClinicalIncremental_20260808")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "timeROC", "ggplot2", "rms", "dcurves")
pkg_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
write.csv(
  data.frame(Package = required_packages, Available = unname(pkg_ok)),
  file.path(out_dir, "Tables", "00_package_status.csv"), row.names = FALSE
)
if (any(!pkg_ok)) {
  missing_pkgs <- required_packages[!pkg_ok]
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install once with:\ninstall.packages(c(",
    paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "),
    "), repos = \"https://cloud.r-project.org\")\n",
    "Then restart R and rerun this whole script from line 1."
  )
}

warning_log <- data.frame(Context = character(), Warning = character())
with_warnlog <- function(context, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_log <<- rbind(
        warning_log,
        data.frame(Context = context, Warning = conditionMessage(w))
      )
      invokeRestart("muffleWarning")
    }
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

sample_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub("\\.cel$", "", x, ignore.case = TRUE)
  gsub("[^a-z0-9]", "", x)
}

read_required_csv <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ":\n", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("Cannot find ", label, ". Checked:\n", paste(paths, collapse = "\n"))
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

extract_item_value <- function(tab, item_pattern) {
  if (!all(c("Item", "Value") %in% names(tab))) return(NA_character_)
  idx <- grep(item_pattern, tab$Item, ignore.case = TRUE)
  if (!length(idx)) return(NA_character_)
  as.character(tab$Value[idx[1]])
}

assert_unique_key <- function(df, key_col, label) {
  k <- df[[key_col]]
  if (any(!nzchar(k) | is.na(k))) stop(label, " contains empty sample keys")
  if (anyDuplicated(k)) {
    dup <- unique(k[duplicated(k)])
    stop(label, " contains duplicated sample keys, e.g. ", paste(head(dup, 5), collapse = ", "))
  }
  invisible(TRUE)
}

cox_concordance <- function(time, event, marker) {
  d <- data.frame(time = safe_num(time), event = safe_num(event), marker = safe_num(marker))
  d <- d[complete.cases(d) & d$time > 0 & d$event %in% c(0, 1), , drop = FALSE]
  if (nrow(d) < 10 || sum(d$event) < 3 || length(unique(d$marker)) < 2) return(NA_real_)
  as.numeric(
    survival::concordance(
      survival::Surv(time, event) ~ marker,
      data = d, reverse = TRUE
    )$concordance
  )
}

time_auc <- function(time, event, marker, dataset, model_name) {
  d <- data.frame(time = safe_num(time), event = safe_num(event), marker = safe_num(marker))
  d <- d[complete.cases(d) & d$time > 0 & d$event %in% c(0, 1), , drop = FALSE]
  tt <- c(365, 1095, 1825)
  obj <- tryCatch(
    with_warnlog(
      paste(dataset, model_name, "timeROC"),
      timeROC::timeROC(
        T = d$time, delta = d$event, marker = d$marker,
        cause = 1, weighting = "marginal", times = tt, iid = FALSE
      )
    ),
    error = function(e) NULL
  )
  auc <- if (is.null(obj)) rep(NA_real_, 3) else as.numeric(obj$AUC)
  data.frame(
    Dataset = dataset, Model = model_name,
    AUC1year = auc[1], AUC3year = auc[2], AUC5year = auc[3],
    MeanAUC = if (all(is.na(auc))) NA_real_ else mean(auc, na.rm = TRUE)
  )
}

cox_coef_table <- function(fit, dataset, model_name) {
  sm <- summary(fit)
  cf <- sm$coefficients
  ci <- sm$conf.int
  data.frame(
    Dataset = dataset,
    Model = model_name,
    Term = rownames(cf),
    HR = ci[, "exp(coef)"],
    Lower95 = ci[, "lower .95"],
    Upper95 = ci[, "upper .95"],
    P = cf[, "Pr(>|z|)"],
    row.names = NULL
  )
}

lrt_upper_tail <- function(reduced_fit, full_fit, dataset, comparison) {
  ll0 <- as.numeric(logLik(reduced_fit))
  ll1 <- as.numeric(logLik(full_fit))
  df0 <- attr(logLik(reduced_fit), "df")
  df1 <- attr(logLik(full_fit), "df")
  chisq <- 2 * (ll1 - ll0)
  ddf <- df1 - df0
  p <- if (is.finite(chisq) && ddf > 0) pchisq(chisq, df = ddf, lower.tail = FALSE) else NA_real_
  data.frame(Dataset = dataset, Comparison = comparison, ChiSquare = chisq, Df = ddf, P = p)
}

performance_locked <- function(df, score_col, group_col, dataset, evidence_role) {
  d <- df[
    complete.cases(df[, c("OS_time", "OS_event", score_col, group_col)]) &
      df$OS_time > 0 & df$OS_event %in% c(0, 1),
    , drop = FALSE
  ]
  score <- safe_num(d[[score_col]])
  group <- factor(as.character(d[[group_col]]), levels = c("Low", "High"))

  fit_cont <- with_warnlog(
    paste(dataset, "continuous Cox"),
    survival::coxph(survival::Surv(OS_time, OS_event) ~ score, data = d)
  )
  sm <- summary(fit_cont)
  fit_group <- with_warnlog(
    paste(dataset, "locked-cutoff Cox"),
    survival::coxph(survival::Surv(OS_time, OS_event) ~ group, data = d)
  )
  smg <- summary(fit_group)
  lr <- with_warnlog(
    paste(dataset, "locked-cutoff log-rank"),
    survival::survdiff(survival::Surv(OS_time, OS_event) ~ group, data = d)
  )
  lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  auc <- time_auc(d$OS_time, d$OS_event, score, dataset, "M070")

  data.frame(
    Dataset = dataset,
    EvidenceRole = evidence_role,
    N = nrow(d), Events = sum(d$OS_event == 1),
    HR_per_1_M070_unit = sm$conf.int[1, "exp(coef)"],
    ContinuousLower95 = sm$conf.int[1, "lower .95"],
    ContinuousUpper95 = sm$conf.int[1, "upper .95"],
    ContinuousP = sm$coefficients[1, "Pr(>|z|)"],
    Cindex = cox_concordance(d$OS_time, d$OS_event, score),
    HighN = sum(group == "High"), LowN = sum(group == "Low"),
    HighEvents = sum(d$OS_event[group == "High"] == 1),
    LowEvents = sum(d$OS_event[group == "Low"] == 1),
    LockedGroupHR_High_vs_Low = smg$conf.int[1, "exp(coef)"],
    LockedGroupLower95 = smg$conf.int[1, "lower .95"],
    LockedGroupUpper95 = smg$conf.int[1, "upper .95"],
    LogRankP = lr_p,
    AUC1year = auc$AUC1year, AUC3year = auc$AUC3year,
    AUC5year = auc$AUC5year, MeanAUC = auc$MeanAUC
  )
}

plot_locked_km <- function(df, group_col, dataset_title, output_png) {
  d <- df[
    complete.cases(df[, c("OS_time", "OS_event", group_col)]) &
      df$OS_time > 0 & df$OS_event %in% c(0, 1),
    , drop = FALSE
  ]
  d$RiskGroup <- factor(as.character(d[[group_col]]), levels = c("Low", "High"))
  fit <- survival::survfit(survival::Surv(OS_time, OS_event) ~ RiskGroup, data = d)
  lr <- survival::survdiff(survival::Surv(OS_time, OS_event) ~ RiskGroup, data = d)
  p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  p_lab <- if (is.finite(p) && p < 0.001) "log-rank P < 0.001" else paste0("log-rank P = ", formatC(p, digits = 3, format = "fg"))

  png(output_png, width = 7, height = 6, units = "in", res = 300, bg = "white")
  op <- par(mar = c(5, 5, 3.5, 1.5), las = 1)
  on.exit({par(op); dev.off()}, add = TRUE)
  plot(
    fit, col = c("#2C7BB6", "#D7191C"), lwd = 2.2, mark.time = TRUE,
    xlab = "Time (days)", ylab = "Overall survival probability",
    main = dataset_title, conf.int = FALSE
  )
  legend(
    "bottomleft",
    legend = c(
      paste0("Low (n = ", sum(d$RiskGroup == "Low"), ")"),
      paste0("High (n = ", sum(d$RiskGroup == "High"), ")"),
      p_lab
    ),
    col = c("#2C7BB6", "#D7191C", NA), lwd = c(2.2, 2.2, NA),
    bty = "n", text.col = "black"
  )
}

bootstrap_cindex <- function(d, clinical_vars, dataset, B = 1000L, seed = 20260808L) {
  set.seed(seed)
  f_m070 <- survival::Surv(OS_time, OS_event) ~ M070_z
  f_clin <- as.formula(paste("survival::Surv(OS_time, OS_event) ~", paste(clinical_vars, collapse = " + ")))
  f_comb <- as.formula(paste("survival::Surv(OS_time, OS_event) ~ M070_z +", paste(clinical_vars, collapse = " + ")))
  res <- matrix(NA_real_, nrow = B, ncol = 5)
  colnames(res) <- c("M070", "Clinical", "Combined", "Combined_minus_Clinical", "Combined_minus_M070")

  for (b in seq_len(B)) {
    idx <- sample.int(nrow(d), size = nrow(d), replace = TRUE)
    db <- d[idx, , drop = FALSE]
    if (sum(db$OS_event == 1) < 3) next
    fits <- tryCatch(
      list(
        M070 = survival::coxph(f_m070, data = db, x = TRUE),
        Clinical = survival::coxph(f_clin, data = db, x = TRUE),
        Combined = survival::coxph(f_comb, data = db, x = TRUE)
      ),
      error = function(e) NULL
    )
    if (is.null(fits)) next
    lp <- lapply(fits, function(f) as.numeric(predict(f, type = "lp")))
    ci <- vapply(lp, function(z) cox_concordance(db$OS_time, db$OS_event, z), numeric(1))
    res[b, 1:3] <- ci[c("M070", "Clinical", "Combined")]
    res[b, 4] <- res[b, 3] - res[b, 2]
    res[b, 5] <- res[b, 3] - res[b, 1]
  }
  raw <- as.data.frame(res)
  raw$Bootstrap <- seq_len(B)
  raw$Dataset <- dataset
  long_summary <- do.call(
    rbind,
    lapply(colnames(res), function(nm) {
      x <- res[, nm]
      x <- x[is.finite(x)]
      data.frame(
        Dataset = dataset, Metric = nm, ValidBootstrapN = length(x),
        Mean = if (length(x)) mean(x) else NA_real_,
        Lower95 = if (length(x)) unname(quantile(x, 0.025, na.rm = TRUE)) else NA_real_,
        Upper95 = if (length(x)) unname(quantile(x, 0.975, na.rm = TRUE)) else NA_real_
      )
    })
  )
  list(raw = raw, summary = long_summary)
}

clinical_analysis <- function(d, clinical_vars, dataset, bootstrap_seed) {
  keep <- c("Sample", "OS_time", "OS_event", "M070", clinical_vars)
  cc <- d[complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
  cc <- cc[cc$OS_time > 0 & cc$OS_event %in% c(0, 1), , drop = FALSE]
  cc$M070_z <- as.numeric(scale(cc$M070))

  f_m070 <- survival::Surv(OS_time, OS_event) ~ M070_z
  f_clin <- as.formula(paste("survival::Surv(OS_time, OS_event) ~", paste(clinical_vars, collapse = " + ")))
  f_comb <- as.formula(paste("survival::Surv(OS_time, OS_event) ~ M070_z +", paste(clinical_vars, collapse = " + ")))
  fit_m070 <- with_warnlog(paste(dataset, "M070 Cox"), survival::coxph(f_m070, data = cc, x = TRUE))
  fit_clin <- with_warnlog(paste(dataset, "Clinical Cox"), survival::coxph(f_clin, data = cc, x = TRUE))
  fit_comb <- with_warnlog(paste(dataset, "Combined Cox"), survival::coxph(f_comb, data = cc, x = TRUE))

  pred <- data.frame(
    M070 = as.numeric(predict(fit_m070, type = "lp")),
    Clinical = as.numeric(predict(fit_clin, type = "lp")),
    Combined = as.numeric(predict(fit_comb, type = "lp"))
  )
  cindex <- data.frame(
    Dataset = dataset,
    Model = c("M070", "Clinical", "Combined"),
    Cindex = c(
      cox_concordance(cc$OS_time, cc$OS_event, pred$M070),
      cox_concordance(cc$OS_time, cc$OS_event, pred$Clinical),
      cox_concordance(cc$OS_time, cc$OS_event, pred$Combined)
    )
  )
  auc <- rbind(
    time_auc(cc$OS_time, cc$OS_event, pred$M070, dataset, "M070"),
    time_auc(cc$OS_time, cc$OS_event, pred$Clinical, dataset, "Clinical"),
    time_auc(cc$OS_time, cc$OS_event, pred$Combined, dataset, "Combined")
  )
  lrt <- rbind(
    lrt_upper_tail(fit_clin, fit_comb, dataset, "Clinical vs Combined: add M070"),
    lrt_upper_tail(fit_m070, fit_comb, dataset, "M070 vs Combined: add clinical covariates")
  )
  coef_tab <- rbind(
    cox_coef_table(fit_m070, dataset, "M070"),
    cox_coef_table(fit_clin, dataset, "Clinical"),
    cox_coef_table(fit_comb, dataset, "Combined")
  )
  boot <- bootstrap_cindex(cc, clinical_vars, dataset, B = 1000L, seed = bootstrap_seed)
  list(
    data = cc, fit_m070 = fit_m070, fit_clin = fit_clin, fit_comb = fit_comb,
    predictions = pred, cindex = cindex, auc = auc, lrt = lrt,
    coef = coef_tab, bootstrap = boot
  )
}

predicted_risk_at_time <- function(fit, newdata, tt) {
  bh <- survival::basehaz(fit, centered = FALSE)
  eligible <- which(bh$time <= tt)
  h0 <- if (length(eligible)) max(bh$hazard[eligible], na.rm = TRUE) else 0
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  pmin(pmax(1 - exp(-h0 * exp(lp)), 0), 1)
}

make_dca <- function(analysis, dataset, file_prefix) {
  d <- analysis$data
  for (tt in c(365, 1095, 1825)) {
    dca_dat <- data.frame(
      OS_time = d$OS_time, OS_event = d$OS_event,
      Clinical = predicted_risk_at_time(analysis$fit_clin, d, tt),
      M070 = predicted_risk_at_time(analysis$fit_m070, d, tt),
      Combined = predicted_risk_at_time(analysis$fit_comb, d, tt)
    )
    obj <- tryCatch(
      with_warnlog(
        paste(dataset, "DCA", tt),
        dcurves::dca(
          survival::Surv(OS_time, OS_event) ~ Clinical + M070 + Combined,
          data = dca_dat, time = tt,
          thresholds = seq(0.01, 0.80, by = 0.01),
          as_probability = c("Clinical", "M070", "Combined")
        )
      ),
      error = function(e) {
        warning_log <<- rbind(
          warning_log,
          data.frame(Context = paste(dataset, "DCA", tt), Warning = conditionMessage(e))
        )
        NULL
      }
    )
    if (!is.null(obj)) {
      p <- plot(obj, smooth = TRUE) +
        ggplot2::labs(
          title = paste0(dataset, ": ", round(tt / 365), "-year decision curve"),
          x = "Threshold probability", y = "Net benefit"
        ) +
        ggplot2::theme_bw(base_size = 11)
      ggplot2::ggsave(
        file.path(out_dir, "Figures", paste0(file_prefix, "_DCA_", round(tt / 365), "year.png")),
        plot = p, width = 7, height = 5.5, dpi = 300, bg = "white"
      )
    }
  }
}

make_nomogram_and_calibration <- function(analysis) {
  d <- analysis$data
  dd <- rms::datadist(d)
  assign("dd_07A", dd, envir = .GlobalEnv)
  old_dd <- getOption("datadist")
  options(datadist = "dd_07A")
  on.exit({
    options(datadist = old_dd)
    if (exists("dd_07A", envir = .GlobalEnv, inherits = FALSE)) rm("dd_07A", envir = .GlobalEnv)
  }, add = TRUE)

  fit <- with_warnlog(
    "E-MTAB-179 rms combined Cox",
    rms::cph(
      survival::Surv(OS_time, OS_event) ~ M070_z + Age_months + Male + INSS_stage4,
      data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095
    )
  )
  sfun <- rms::Survival(fit)
  nom <- rms::nomogram(
    fit,
    fun = list(
      function(lp) sfun(365, lp),
      function(lp) sfun(1095, lp),
      function(lp) sfun(1825, lp)
    ),
    funlabel = c("1-year overall survival", "3-year overall survival", "5-year overall survival"),
    lp = FALSE
  )
  pdf(file.path(out_dir, "Figures", "07_EMTAB179_M070_combined_nomogram.pdf"), width = 11, height = 7)
  plot(nom, xfrac = 0.45)
  dev.off()

  for (tt in c(365, 1095, 1825)) {
    cal <- tryCatch(
      with_warnlog(
        paste("E-MTAB-179 calibration", tt),
        rms::calibrate(fit, cmethod = "KM", method = "boot", u = tt, B = 200)
      ),
      error = function(e) {
        warning_log <<- rbind(
          warning_log,
          data.frame(Context = paste("E-MTAB-179 calibration", tt), Warning = conditionMessage(e))
        )
        NULL
      }
    )
    if (!is.null(cal)) {
      png(
        file.path(out_dir, "Figures", paste0("08_EMTAB179_calibration_", round(tt / 365), "year.png")),
        width = 6, height = 6, units = "in", res = 300, bg = "white"
      )
      plot(
        cal, lwd = 2, subtitles = FALSE,
        xlab = paste0("Predicted ", round(tt / 365), "-year overall survival"),
        ylab = paste0("Observed ", round(tt / 365), "-year overall survival")
      )
      abline(0, 1, lty = 2, col = "grey50")
      dev.off()
    }
  }
  invisible(fit)
}

## ----------------------------- Inputs ---------------------------------------
path_05_summary <- file.path(
  root_dir, "SC28_RECALC_05_RepeatedCV_SelectionAudit_20260808", "Tables",
  "08_model_stability_summary.csv"
)
path_06_scores <- file.path(
  root_dir, "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Tables",
  "06_EMTAB_locked_M070_scores.csv"
)
path_06_lock_record <- file.path(
  root_dir, "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Tables",
  "05_LOCK_RECORD_BEFORE_SEQC.csv"
)
path_06_model <- file.path(
  root_dir, "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Models",
  "M070_locked_model_BEFORE_SEQC.rds"
)
path_06b_scores <- file.path(
  root_dir, "SC28_RECALC_06B_M070_SEQC_CohortZ_External_20260808", "Tables",
  "06_SEQC_EXTERNAL_M070_scores_with_outcomes.csv"
)
path_06b_model <- file.path(
  root_dir, "SC28_RECALC_06B_M070_SEQC_CohortZ_External_20260808", "Models",
  "M070_locked_model_USED_UNCHANGED.rds"
)

clinical_base <- file.path(root_dir, "SC26B_M050_ClinicalIncrementalValue_MergeClinical", "Tables")
path_emtab_clin <- first_existing(
  c(
    file.path(clinical_base, "02_EMTAB179_merged_risk_clinical.csv"),
    file.path(root_dir, "SC26B_M050_ClinicalIncrementalValue_MergeClinical", "02_EMTAB179_merged_risk_clinical.csv")
  ),
  "prior curated E-MTAB-179 clinical merge"
)
path_seqc_clin <- first_existing(
  c(
    file.path(clinical_base, "03_SEQC_merged_risk_clinical.csv"),
    file.path(root_dir, "SC26B_M050_ClinicalIncrementalValue_MergeClinical", "03_SEQC_merged_risk_clinical.csv")
  ),
  "prior curated SEQC clinical merge"
)

for (p in c(path_05_summary, path_06_scores, path_06_lock_record, path_06_model, path_06b_scores, path_06b_model)) {
  if (!file.exists(p)) stop("Missing required RECALC input:\n", p)
}

## ---------------------- Frozen-model integrity audit -------------------------
expected_model_md5 <- "e719d0e6b626e75cc134e83fd0929083"
md5_06 <- unname(tools::md5sum(path_06_model))
md5_06b <- unname(tools::md5sum(path_06b_model))
integrity <- data.frame(
  Check = c(
    "RECALC06 locked model MD5",
    "RECALC06B copied model MD5",
    "06 equals 06B",
    "06 equals expected audited MD5"
  ),
  Value = c(md5_06, md5_06b, as.character(identical(md5_06, md5_06b)), as.character(identical(md5_06, expected_model_md5)))
)
write.csv(integrity, file.path(out_dir, "Tables", "01_frozen_model_integrity.csv"), row.names = FALSE)
if (!identical(md5_06, md5_06b) || !identical(md5_06, expected_model_md5)) {
  stop("Frozen M070 model integrity check failed. Do not continue.")
}

stability <- read_required_csv(path_05_summary, "RECALC05 model stability summary")
model_ref <- stability[stability$Model %in% c("M070_SuperPC_RSF", "M050_SuperPC_Ridge"), , drop = FALSE]
if (nrow(model_ref) != 2) stop("Cannot recover both M070 and M050 reference rows from RECALC05")
write.csv(model_ref, file.path(out_dir, "Tables", "02_model_selection_reference_from_repeatedCV.csv"), row.names = FALSE)

lock_record <- read_required_csv(path_06_lock_record, "RECALC06 lock record")
locked_cutoff <- safe_num(extract_item_value(lock_record, "Locked risk cutoff"))
if (!is.finite(locked_cutoff)) stop("Cannot recover locked E-MTAB-179 risk cutoff")

## -------------------- Final fixed M070 score datasets ------------------------
emtab <- read_required_csv(path_06_scores, "RECALC06 E-MTAB-179 M070 scores")
seqc <- read_required_csv(path_06b_scores, "RECALC06B SEQC M070 scores")
required_emtab <- c("Sample", "OS_time", "OS_event", "M070_RiskScore", "RiskGroup_LockedCutoff")
required_seqc <- c("Sample", "OS_time", "OS_event", "Corrected06B_M070_RiskScore", "RiskGroup_Locked_EMTAB_Cutoff")
if (!all(required_emtab %in% names(emtab))) stop("Unexpected E-MTAB-179 score table columns")
if (!all(required_seqc %in% names(seqc))) stop("Unexpected SEQC 06B score table columns")

emtab_final <- data.frame(
  Sample = as.character(emtab$Sample),
  OS_time = safe_num(emtab$OS_time), OS_event = safe_num(emtab$OS_event),
  M070 = safe_num(emtab$M070_RiskScore),
  RiskGroup = as.character(emtab$RiskGroup_LockedCutoff)
)
seqc_final <- data.frame(
  Sample = as.character(seqc$Sample),
  OS_time = safe_num(seqc$OS_time), OS_event = safe_num(seqc$OS_event),
  M070 = safe_num(seqc$Corrected06B_M070_RiskScore),
  RiskGroup = as.character(seqc$RiskGroup_Locked_EMTAB_Cutoff)
)
emtab_final$SampleKey07A <- sample_key(emtab_final$Sample)
seqc_final$SampleKey07A <- sample_key(seqc_final$Sample)
assert_unique_key(emtab_final, "SampleKey07A", "E-MTAB-179 final M070")
assert_unique_key(seqc_final, "SampleKey07A", "SEQC final M070")

if (nrow(emtab_final) != 478 || sum(emtab_final$OS_event == 1, na.rm = TRUE) != 91) {
  stop("Unexpected E-MTAB-179 full score population; expected N=478/events=91")
}
if (nrow(seqc_final) != 498 || sum(seqc_final$OS_event == 1, na.rm = TRUE) != 105) {
  stop("Unexpected SEQC full score population; expected N=498/events=105")
}

## Recompute locked groups only as an integrity check, never as a new cutoff.
group_emtab_check <- ifelse(emtab_final$M070 > locked_cutoff, "High", "Low")
group_seqc_check <- ifelse(seqc_final$M070 > locked_cutoff, "High", "Low")
if (!all(group_emtab_check == emtab_final$RiskGroup) || !all(group_seqc_check == seqc_final$RiskGroup)) {
  stop("Stored risk groups do not match the frozen E-MTAB-179 cutoff")
}

## ------------------------ Clinical merge ------------------------------------
emtab_old <- read_required_csv(path_emtab_clin, "prior E-MTAB-179 clinical merge")
seqc_old <- read_required_csv(path_seqc_clin, "prior SEQC clinical merge")

find_sample_col <- function(df, label) {
  candidates <- c("Sample", "SampleID", "sample", "sample_id", "ID", "id")
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit) && "SampleKey" %in% names(df)) hit <- "SampleKey"
  if (!length(hit)) stop("Cannot identify sample column in ", label)
  hit[1]
}

emtab_sample_col <- find_sample_col(emtab_old, "prior E-MTAB-179 clinical merge")
seqc_sample_col <- find_sample_col(seqc_old, "prior SEQC clinical merge")

if (!all(c("Age_months", "Male", "INSS_stage4") %in% names(emtab_old))) {
  stop("Prior E-MTAB-179 clinical merge lacks Age_months/Male/INSS_stage4")
}
if (!all(c("Age_months", "Male") %in% names(seqc_old))) {
  stop("Prior SEQC clinical merge lacks Age_months/Male")
}

emtab_clin <- data.frame(
  SampleKey07A = sample_key(emtab_old[[emtab_sample_col]]),
  Age_months = safe_num(emtab_old$Age_months),
  Male = safe_num(emtab_old$Male),
  INSS_stage4 = safe_num(emtab_old$INSS_stage4)
)
seqc_clin <- data.frame(
  SampleKey07A = sample_key(seqc_old[[seqc_sample_col]]),
  Age_months = safe_num(seqc_old$Age_months),
  Male = safe_num(seqc_old$Male)
)
assert_unique_key(emtab_clin, "SampleKey07A", "prior E-MTAB-179 clinical merge")
assert_unique_key(seqc_clin, "SampleKey07A", "prior SEQC clinical merge")

emtab_merged <- merge(emtab_final, emtab_clin, by = "SampleKey07A", all.x = TRUE, sort = FALSE)
seqc_merged <- merge(seqc_final, seqc_clin, by = "SampleKey07A", all.x = TRUE, sort = FALSE)
if (nrow(emtab_merged) != nrow(emtab_final) || nrow(seqc_merged) != nrow(seqc_final)) {
  stop("Clinical merge changed patient counts")
}

## Restore the final-score order after merge.
emtab_merged <- emtab_merged[match(emtab_final$SampleKey07A, emtab_merged$SampleKey07A), , drop = FALSE]
seqc_merged <- seqc_merged[match(seqc_final$SampleKey07A, seqc_merged$SampleKey07A), , drop = FALSE]

emtab_cc <- complete.cases(emtab_merged[, c("OS_time", "OS_event", "M070", "Age_months", "Male", "INSS_stage4")])
seqc_cc <- complete.cases(seqc_merged[, c("OS_time", "OS_event", "M070", "Age_months", "Male")])
clinical_audit <- data.frame(
  Dataset = c("E-MTAB-179", "SEQC"),
  FullScoreN = c(nrow(emtab_merged), nrow(seqc_merged)),
  FullScoreEvents = c(sum(emtab_merged$OS_event == 1), sum(seqc_merged$OS_event == 1)),
  ClinicalCompleteCaseN = c(sum(emtab_cc), sum(seqc_cc)),
  ClinicalCompleteCaseEvents = c(sum(emtab_merged$OS_event[emtab_cc] == 1), sum(seqc_merged$OS_event[seqc_cc] == 1)),
  ClinicalCovariates = c("Age_months + Male + INSS_stage4", "Age_months + Male")
)
write.csv(clinical_audit, file.path(out_dir, "Tables", "03_clinical_merge_and_complete_case_audit.csv"), row.names = FALSE)
if (sum(emtab_cc) != 389 || sum(emtab_merged$OS_event[emtab_cc] == 1) != 57) {
  stop("E-MTAB-179 clinical complete-case population drifted from audited N=389/events=57")
}
if (sum(seqc_cc) != 498 || sum(seqc_merged$OS_event[seqc_cc] == 1) != 105) {
  stop("SEQC clinical complete-case population drifted from audited N=498/events=105")
}

write.csv(emtab_merged, file.path(out_dir, "Tables", "04_EMTAB179_final_M070_with_curated_clinical.csv"), row.names = FALSE)
write.csv(seqc_merged, file.path(out_dir, "Tables", "05_SEQC_final_M070_with_curated_clinical.csv"), row.names = FALSE)

## -------------------- Locked prognostic performance -------------------------
perf_emtab <- performance_locked(
  emtab_merged, "M070", "RiskGroup", "E-MTAB-179",
  "Apparent full-development-cohort description; model-selection evidence is RECALC05 repeated CV"
)
perf_seqc <- performance_locked(
  seqc_merged, "M070", "RiskGroup", "SEQC",
  "External validation of frozen M070 after expression-only cohort-wise z harmonization"
)
perf <- rbind(perf_emtab, perf_seqc)
write.csv(perf, file.path(out_dir, "Tables", "06_final_M070_prognostic_performance.csv"), row.names = FALSE)

plot_locked_km(
  emtab_merged, "RiskGroup", "E-MTAB-179: locked M070 risk stratification",
  file.path(out_dir, "Figures", "01_EMTAB179_locked_M070_KM.png")
)
plot_locked_km(
  seqc_merged, "RiskGroup", "SEQC: external validation of locked M070 risk stratification",
  file.path(out_dir, "Figures", "02_SEQC_external_locked_M070_KM.png")
)

## ---------------- Clinical incremental value -------------------------------
ana_emtab <- clinical_analysis(
  emtab_merged, c("Age_months", "Male", "INSS_stage4"), "E-MTAB-179", 20260881L
)
ana_seqc <- clinical_analysis(
  seqc_merged, c("Age_months", "Male"), "SEQC", 20260882L
)

cindex_tab <- rbind(ana_emtab$cindex, ana_seqc$cindex)
auc_tab <- rbind(ana_emtab$auc, ana_seqc$auc)
lrt_tab <- rbind(ana_emtab$lrt, ana_seqc$lrt)
coef_tab <- rbind(ana_emtab$coef, ana_seqc$coef)
boot_summary <- rbind(ana_emtab$bootstrap$summary, ana_seqc$bootstrap$summary)
write.csv(cindex_tab, file.path(out_dir, "Tables", "07_clinical_incremental_Cindex.csv"), row.names = FALSE)
write.csv(auc_tab, file.path(out_dir, "Tables", "08_clinical_incremental_timeAUC.csv"), row.names = FALSE)
write.csv(lrt_tab, file.path(out_dir, "Tables", "09_clinical_incremental_likelihood_ratio_tests.csv"), row.names = FALSE)
write.csv(coef_tab, file.path(out_dir, "Tables", "10_clinical_incremental_Cox_coefficients.csv"), row.names = FALSE)
write.csv(boot_summary, file.path(out_dir, "Tables", "11_bootstrap_Cindex_95CI_and_differences.csv"), row.names = FALSE)
write.csv(ana_emtab$bootstrap$raw, file.path(out_dir, "Tables", "12_EMTAB179_bootstrap_Cindex_raw.csv"), row.names = FALSE)
write.csv(ana_seqc$bootstrap$raw, file.path(out_dir, "Tables", "13_SEQC_bootstrap_Cindex_raw.csv"), row.names = FALSE)

## Publication-facing comparison figures.
p_ci <- ggplot2::ggplot(cindex_tab, ggplot2::aes(x = Model, y = Cindex, fill = Model)) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::facet_wrap(~Dataset) +
  ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
  ggplot2::scale_fill_manual(values = c(M070 = "#4C78A8", Clinical = "#9E9E9E", Combined = "#E45756")) +
  ggplot2::labs(x = NULL, y = "Harrell C-index", title = "Incremental prognostic discrimination of M070") +
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "none")
ggplot2::ggsave(file.path(out_dir, "Figures", "03_clinical_incremental_Cindex.png"), p_ci, width = 7.5, height = 5, dpi = 300, bg = "white")

auc_long <- rbind(
  data.frame(Dataset = auc_tab$Dataset, Model = auc_tab$Model, Time = "1 year", AUC = auc_tab$AUC1year),
  data.frame(Dataset = auc_tab$Dataset, Model = auc_tab$Model, Time = "3 years", AUC = auc_tab$AUC3year),
  data.frame(Dataset = auc_tab$Dataset, Model = auc_tab$Model, Time = "5 years", AUC = auc_tab$AUC5year)
)
auc_long$Time <- factor(auc_long$Time, levels = c("1 year", "3 years", "5 years"))
p_auc <- ggplot2::ggplot(auc_long, ggplot2::aes(x = Time, y = AUC, color = Model, group = Model)) +
  ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~Dataset) + ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
  ggplot2::scale_color_manual(values = c(M070 = "#4C78A8", Clinical = "#7F7F7F", Combined = "#E45756")) +
  ggplot2::labs(x = NULL, y = "Time-dependent AUC", title = "Time-dependent prognostic discrimination") +
  ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "top")
ggplot2::ggsave(file.path(out_dir, "Figures", "04_clinical_incremental_timeAUC.png"), p_auc, width = 7.5, height = 5, dpi = 300, bg = "white")

## Main development-cohort nomogram/calibration; DCA in both cohorts.
make_nomogram_and_calibration(ana_emtab)
make_dca(ana_emtab, "E-MTAB-179", "05_EMTAB179")
make_dca(ana_seqc, "SEQC", "06_SEQC")

## ---------------------------- Final audit -----------------------------------
final_audit <- data.frame(
  Item = c(
    "Final model", "Frozen model MD5", "Model selection source",
    "SEQC role", "SEQC expression harmonization", "Model refitted in SEQC",
    "Features reselected in SEQC", "Cutoff tuned in SEQC", "Locked cutoff source",
    "Locked cutoff", "Clinical score transform", "E-MTAB clinical covariates",
    "SEQC clinical covariates", "Bootstrap repetitions", "Calibration bootstrap repetitions",
    "Interpretation caution"
  ),
  Value = c(
    "M070 SuperPC -> RSF", md5_06,
    "E-MTAB-179 repeated 10 x 5-fold CV in RECALC05",
    "External validation",
    "Within-SEQC gene-wise z standardization using expression only (RECALC06B)",
    "FALSE", "FALSE", "FALSE", "E-MTAB-179 locked cutoff from RECALC06",
    as.character(locked_cutoff),
    "Within-cohort z standardization of the already-fixed M070 score for clinical Cox models only",
    "Age_months + Male + INSS_stage4", "Age_months + Male", "1000", "200",
    "SEQC contains limited clinical covariates; claim additional information beyond available covariates, not independence from all established risk factors"
  )
)
write.csv(final_audit, file.path(out_dir, "Tables", "14_FINAL_07A_AUDIT_SUMMARY.csv"), row.names = FALSE)

write.csv(warning_log, file.path(out_dir, "Logs", "01_warnings_captured.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Logs", "02_sessionInfo.txt"))
writeLines(
  c(
    "SC28 RECALC 07A completed.",
    paste0("Frozen M070 MD5: ", md5_06),
    paste0("Locked E-MTAB-179 cutoff: ", format(locked_cutoff, digits = 15)),
    "E-MTAB-179 full score population: N=478, events=91.",
    "SEQC full external population: N=498, events=105.",
    "Clinical complete cases: E-MTAB-179 N=389/events=57; SEQC N=498/events=105.",
    "SEQC uses the frozen M070 with expression-only cohort-wise gene z harmonization from RECALC06B.",
    "Primary risk-group analyses use the locked E-MTAB-179 cutoff; no SEQC outcome-guided cutoff is used.",
    "Clinical M070_z is only a within-cohort affine rescaling of the fixed score for Cox interpretability.",
    "Important wording: SEQC supports prognostic information beyond the clinical variables available in that cohort; do not claim independence from all established neuroblastoma risk factors."
  ),
  file.path(out_dir, "Logs", "03_run_summary.txt")
)

cat("\n============================================================\n")
cat("SC28 RECALC 07A completed successfully.\n")
cat("Output: ", out_dir, "\n", sep = "")
cat("Frozen M070 MD5 verified: ", md5_06, "\n", sep = "")
cat("Clinical complete cases verified: E-MTAB-179 389/57; SEQC 498/105.\n")
cat("Warnings captured: ", nrow(warning_log), "\n", sep = "")
cat("Upload the ENTIRE 07A output folder as ZIP for strict audit before manuscript replacement.\n")
cat("============================================================\n")
