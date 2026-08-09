rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================================
## PRDX4_RECALC_07A_FIX1_ClinicalDCA_AgeUnit.R
##
## Targeted post-audit correction for RECALC 07A.
##
## This script DOES NOT touch:
##   - M070 model fitting or feature selection
##   - the frozen M070 model / MD5
##   - E-MTAB-179 repeated-CV model-selection evidence
##   - SEQC expression harmonization or external M070 risk scores
##   - the locked E-MTAB-179 cutoff
##   - the completed 1000-bootstrap C-index audit
##
## It corrects only:
##   A. E-MTAB-179 age unit. The old curated column called Age_months actually
##      contains day-scale values (audited range 0--8983; median 401). It is
##      retained as Age_days_source and converted to months using 365.25/12.
##   B. DCA. RECALC07A passed already-computed Cox probabilities to dca() and
##      also supplied as_probability=, which asks dcurves to convert the marker
##      to a probability a second time. FIX1 supplies the probabilities
##      directly and does NOT use as_probability=.
##   C. Nomogram labels use the corrected age-in-months variable.
##   D. The unstable 1/3/5-year rms bootstrap calibration figures from 07A are
##      NOT carried forward. A transparent, apparent 3-year grouped calibration
##      is generated for exploratory display only.
## ============================================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
source_07a <- file.path(root_dir, "SC28_RECALC_07A_M070_Prognostic_ClinicalIncremental_20260808")
out_dir <- file.path(root_dir, "SC28_RECALC_07A_FIX1_ClinicalDCA_AgeUnit_20260808")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

required_packages <- c("survival", "timeROC", "ggplot2", "rms", "dcurves", "tibble")
pkg_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
write.csv(data.frame(Package = required_packages, Available = unname(pkg_ok)),
          file.path(out_dir, "Tables", "00_package_status.csv"), row.names = FALSE)
if (any(!pkg_ok)) {
  miss <- required_packages[!pkg_ok]
  stop("Missing packages: ", paste(miss, collapse = ", "))
}

warning_log <- data.frame(Context = character(), Warning = character())
with_warnlog <- function(context, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_log <<- rbind(warning_log,
                            data.frame(Context = context, Warning = conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
read_req <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ":\n", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

## ------------------------ Load and verify 07A -------------------------------
path_audit <- file.path(source_07a, "Tables", "14_FINAL_07A_AUDIT_SUMMARY.csv")
path_emtab <- file.path(source_07a, "Tables", "04_EMTAB179_final_M070_with_curated_clinical.csv")
path_seqc <- file.path(source_07a, "Tables", "05_SEQC_final_M070_with_curated_clinical.csv")
path_boot <- file.path(source_07a, "Tables", "11_bootstrap_Cindex_95CI_and_differences.csv")

audit07a <- read_req(path_audit, "07A final audit")
emtab0 <- read_req(path_emtab, "07A E-MTAB-179 merged table")
seqc0 <- read_req(path_seqc, "07A SEQC merged table")
boot07a <- read_req(path_boot, "07A bootstrap summary")

get_audit <- function(item) {
  z <- audit07a$Value[audit07a$Item == item]
  if (!length(z)) NA_character_ else as.character(z[1])
}
expected_md5 <- "e719d0e6b626e75cc134e83fd0929083"
if (!identical(get_audit("Frozen model MD5"), expected_md5)) {
  stop("07A frozen model MD5 is not the audited M070 MD5. STOP.")
}
if (!identical(get_audit("Model refitted in SEQC"), "FALSE") ||
    !identical(get_audit("Features reselected in SEQC"), "FALSE") ||
    !identical(get_audit("Cutoff tuned in SEQC"), "FALSE")) {
  stop("07A external-validation safeguards failed. STOP.")
}

if (nrow(emtab0) != 478 || sum(safe_num(emtab0$OS_event) == 1, na.rm = TRUE) != 91) {
  stop("Unexpected E-MTAB-179 population")
}
if (nrow(seqc0) != 498 || sum(safe_num(seqc0$OS_event) == 1, na.rm = TRUE) != 105) {
  stop("Unexpected SEQC population")
}

## ------------------------ Correct age units ---------------------------------
## Audit evidence: the E-MTAB column named Age_months has day-scale integer
## values. Keep the source values, convert explicitly, and record the check.
emtab0$Age_days_source <- safe_num(emtab0$Age_months)
age_audit_before <- c(
  min = min(emtab0$Age_days_source, na.rm = TRUE),
  median = median(emtab0$Age_days_source, na.rm = TRUE),
  max = max(emtab0$Age_days_source, na.rm = TRUE)
)
if (!(age_audit_before["median"] > 100 && age_audit_before["max"] > 1000)) {
  stop("E-MTAB age values no longer look day-scaled. Do not apply automatic conversion.")
}
days_per_month <- 365.25 / 12
emtab0$Age_months <- emtab0$Age_days_source / days_per_month
seqc0$Age_months <- safe_num(seqc0$Age_months)

age_unit_audit <- data.frame(
  Dataset = c("E-MTAB-179 source", "E-MTAB-179 corrected", "SEQC unchanged"),
  Variable = c("Age_days_source", "Age_months", "Age_months"),
  Unit = c("days", "months", "months"),
  N = c(sum(is.finite(emtab0$Age_days_source)), sum(is.finite(emtab0$Age_months)), sum(is.finite(seqc0$Age_months))),
  Min = c(min(emtab0$Age_days_source, na.rm = TRUE), min(emtab0$Age_months, na.rm = TRUE), min(seqc0$Age_months, na.rm = TRUE)),
  Median = c(median(emtab0$Age_days_source, na.rm = TRUE), median(emtab0$Age_months, na.rm = TRUE), median(seqc0$Age_months, na.rm = TRUE)),
  Max = c(max(emtab0$Age_days_source, na.rm = TRUE), max(emtab0$Age_months, na.rm = TRUE), max(seqc0$Age_months, na.rm = TRUE)),
  Conversion = c("none; retained for provenance", "Age_days_source / (365.25/12)", "none")
)
write.csv(age_unit_audit, file.path(out_dir, "Tables", "01_age_unit_audit_and_correction.csv"), row.names = FALSE)

## Keep corrected merged tables for provenance.
write.csv(emtab0, file.path(out_dir, "Tables", "02_EMTAB179_M070_clinical_age_corrected.csv"), row.names = FALSE)
write.csv(seqc0, file.path(out_dir, "Tables", "03_SEQC_M070_clinical_unchanged.csv"), row.names = FALSE)

## ---------------- Clinical analysis helpers ---------------------------------
cox_concordance <- function(time, event, marker) {
  d <- data.frame(time = safe_num(time), event = safe_num(event), marker = safe_num(marker))
  d <- d[complete.cases(d) & d$time > 0 & d$event %in% c(0, 1), , drop = FALSE]
  as.numeric(survival::concordance(survival::Surv(time, event) ~ marker,
                                   data = d, reverse = TRUE)$concordance)
}

time_auc <- function(time, event, marker, dataset, model) {
  tt <- c(365, 1095, 1825)
  obj <- timeROC::timeROC(T = time, delta = event, marker = marker, cause = 1,
                          weighting = "marginal", times = tt, iid = FALSE)
  data.frame(Dataset = dataset, Model = model,
             AUC1year = obj$AUC[1], AUC3year = obj$AUC[2], AUC5year = obj$AUC[3])
}

coef_table <- function(fit, dataset, model) {
  s <- summary(fit)
  data.frame(
    Dataset = dataset, Model = model, Term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"], Lower95 = s$conf.int[, "lower .95"],
    Upper95 = s$conf.int[, "upper .95"], P = s$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
}

lrt <- function(reduced, full, dataset, comparison) {
  ll0 <- as.numeric(logLik(reduced)); ll1 <- as.numeric(logLik(full))
  df0 <- attr(logLik(reduced), "df"); df1 <- attr(logLik(full), "df")
  chi <- 2 * (ll1 - ll0); ddf <- df1 - df0
  data.frame(Dataset = dataset, Comparison = comparison,
             ChiSquare = chi, Df = ddf,
             P = pchisq(chi, df = ddf, lower.tail = FALSE))
}

prepare_and_fit <- function(d, clinical_vars, dataset) {
  keep <- c("Sample", "OS_time", "OS_event", "M070", clinical_vars)
  x <- d[complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
  x$OS_time <- safe_num(x$OS_time); x$OS_event <- safe_num(x$OS_event); x$M070 <- safe_num(x$M070)
  x <- x[x$OS_time > 0 & x$OS_event %in% c(0, 1), , drop = FALSE]
  x$M070_z <- as.numeric(scale(x$M070))
  fm <- survival::coxph(survival::Surv(OS_time, OS_event) ~ M070_z, data = x, x = TRUE)
  fc <- survival::coxph(as.formula(paste("survival::Surv(OS_time, OS_event) ~", paste(clinical_vars, collapse = " + "))), data = x, x = TRUE)
  fb <- survival::coxph(as.formula(paste("survival::Surv(OS_time, OS_event) ~ M070_z +", paste(clinical_vars, collapse = " + "))), data = x, x = TRUE)
  pred <- data.frame(M070 = as.numeric(predict(fm, type = "lp")),
                     Clinical = as.numeric(predict(fc, type = "lp")),
                     Combined = as.numeric(predict(fb, type = "lp")))
  list(data = x, fit_m070 = fm, fit_clin = fc, fit_comb = fb, pred = pred, dataset = dataset)
}

emtab <- prepare_and_fit(emtab0, c("Age_months", "Male", "INSS_stage4"), "E-MTAB-179")
seqc <- prepare_and_fit(seqc0, c("Age_months", "Male"), "SEQC")
if (nrow(emtab$data) != 389 || sum(emtab$data$OS_event == 1) != 57) stop("E-MTAB complete-case drift")
if (nrow(seqc$data) != 498 || sum(seqc$data$OS_event == 1) != 105) stop("SEQC complete-case drift")

coef_fix <- rbind(
  coef_table(emtab$fit_m070, "E-MTAB-179", "M070"),
  coef_table(emtab$fit_clin, "E-MTAB-179", "Clinical"),
  coef_table(emtab$fit_comb, "E-MTAB-179", "Combined"),
  coef_table(seqc$fit_m070, "SEQC", "M070"),
  coef_table(seqc$fit_clin, "SEQC", "Clinical"),
  coef_table(seqc$fit_comb, "SEQC", "Combined")
)
write.csv(coef_fix, file.path(out_dir, "Tables", "04_clinical_Cox_coefficients_correct_age_unit.csv"), row.names = FALSE)

cindex_fix <- rbind(
  data.frame(Dataset = "E-MTAB-179", Model = names(emtab$pred),
             Cindex = vapply(emtab$pred, function(z) cox_concordance(emtab$data$OS_time, emtab$data$OS_event, z), numeric(1))),
  data.frame(Dataset = "SEQC", Model = names(seqc$pred),
             Cindex = vapply(seqc$pred, function(z) cox_concordance(seqc$data$OS_time, seqc$data$OS_event, z), numeric(1)))
)
write.csv(cindex_fix, file.path(out_dir, "Tables", "05_Cindex_after_age_unit_correction.csv"), row.names = FALSE)

auc_fix <- rbind(
  time_auc(emtab$data$OS_time, emtab$data$OS_event, emtab$pred$M070, "E-MTAB-179", "M070"),
  time_auc(emtab$data$OS_time, emtab$data$OS_event, emtab$pred$Clinical, "E-MTAB-179", "Clinical"),
  time_auc(emtab$data$OS_time, emtab$data$OS_event, emtab$pred$Combined, "E-MTAB-179", "Combined"),
  time_auc(seqc$data$OS_time, seqc$data$OS_event, seqc$pred$M070, "SEQC", "M070"),
  time_auc(seqc$data$OS_time, seqc$data$OS_event, seqc$pred$Clinical, "SEQC", "Clinical"),
  time_auc(seqc$data$OS_time, seqc$data$OS_event, seqc$pred$Combined, "SEQC", "Combined")
)
write.csv(auc_fix, file.path(out_dir, "Tables", "06_timeAUC_after_age_unit_correction.csv"), row.names = FALSE)

lrt_fix <- rbind(
  lrt(emtab$fit_clin, emtab$fit_comb, "E-MTAB-179", "Clinical vs Combined: add M070"),
  lrt(emtab$fit_m070, emtab$fit_comb, "E-MTAB-179", "M070 vs Combined: add clinical covariates"),
  lrt(seqc$fit_clin, seqc$fit_comb, "SEQC", "Clinical vs Combined: add M070"),
  lrt(seqc$fit_m070, seqc$fit_comb, "SEQC", "M070 vs Combined: add clinical covariates")
)
write.csv(lrt_fix, file.path(out_dir, "Tables", "07_LRT_after_age_unit_correction.csv"), row.names = FALSE)

## The 07A bootstrap C-index results remain valid because dividing age by a
## positive constant is an invertible linear change of units in the same Cox
## model. Carry the audited table forward with an explicit provenance note.
write.csv(boot07a, file.path(out_dir, "Tables", "08_bootstrap_Cindex_from_07A_UNCHANGED.csv"), row.names = FALSE)

## ---------------- Correct DCA: probability inputs stay probabilities --------
predicted_risk <- function(fit, newdata, tt) {
  bh <- survival::basehaz(fit, centered = FALSE)
  idx <- which(bh$time <= tt)
  h0 <- if (length(idx)) max(bh$hazard[idx], na.rm = TRUE) else 0
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  pmin(pmax(1 - exp(-h0 * exp(lp)), 0), 1)
}

run_dca <- function(a, prefix) {
  d <- a$data
  all_summary <- list()
  for (tt in c(365, 1095, 1825)) {
    dd <- data.frame(
      OS_time = d$OS_time, OS_event = d$OS_event,
      Clinical = predicted_risk(a$fit_clin, d, tt),
      M070 = predicted_risk(a$fit_m070, d, tt),
      Combined = predicted_risk(a$fit_comb, d, tt)
    )
    yr <- round(tt / 365)
    ps <- do.call(rbind, lapply(c("Clinical", "M070", "Combined"), function(nm) {
      z <- dd[[nm]]
      data.frame(Dataset = a$dataset, TimeYear = yr, Model = nm,
                 Min = min(z), Q25 = unname(quantile(z, .25)), Median = median(z),
                 Q75 = unname(quantile(z, .75)), Max = max(z), Mean = mean(z))
    }))
    all_summary[[as.character(tt)]] <- ps

    ## IMPORTANT: dd already contains predicted probabilities. Do NOT pass
    ## as_probability= here; that option is only for converting non-probability
    ## markers inside dca().
    obj <- with_warnlog(
      paste(a$dataset, "corrected DCA", yr, "year"),
      dcurves::dca(
        survival::Surv(OS_time, OS_event) ~ Clinical + M070 + Combined,
        data = dd, time = tt, thresholds = seq(0.01, 0.50, by = 0.01),
        label = list(Clinical = "Clinical", M070 = "M070", Combined = "Combined")
      )
    )
    dca_raw <- as.data.frame(tibble::as_tibble(obj))
    write.csv(dca_raw,
              file.path(out_dir, "Tables", paste0(prefix, "_DCA_", yr, "year_raw.csv")),
              row.names = FALSE)
    p <- plot(obj, smooth = FALSE) +
      ggplot2::labs(title = paste0(a$dataset, ": ", yr, "-year decision curve"),
                    x = "Threshold probability", y = "Net benefit") +
      ggplot2::coord_cartesian(xlim = c(0.01, 0.50)) +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(out_dir, "Figures", paste0(prefix, "_DCA_", yr, "year_corrected.png")),
                    p, width = 7, height = 5.5, dpi = 300, bg = "white")
  }
  do.call(rbind, all_summary)
}

dca_prob_summary <- rbind(run_dca(emtab, "01_EMTAB179"), run_dca(seqc, "02_SEQC"))
write.csv(dca_prob_summary, file.path(out_dir, "Tables", "09_DCA_input_probability_summary.csv"), row.names = FALSE)

## ---------------- Corrected E-MTAB nomogram ---------------------------------
dd_rms <- rms::datadist(emtab$data)
assign("dd_fix1", dd_rms, envir = .GlobalEnv)
old_dd <- getOption("datadist")
options(datadist = "dd_fix1")

fit_rms <- with_warnlog(
  "corrected nomogram fit",
  rms::cph(survival::Surv(OS_time, OS_event) ~ M070_z + Age_months + Male + INSS_stage4,
           data = emtab$data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1095,
           iter.max = 50)
)
sfun <- rms::Survival(fit_rms)
nom <- rms::nomogram(
  fit_rms,
  fun = list(function(lp) sfun(365, lp), function(lp) sfun(1095, lp), function(lp) sfun(1825, lp)),
  funlabel = c("1-year overall survival", "3-year overall survival", "5-year overall survival"),
  lp = FALSE
)
pdf(file.path(out_dir, "Figures", "03_EMTAB179_M070_combined_nomogram_age_months_corrected.pdf"),
    width = 11, height = 7)
plot(nom, xfrac = 0.40)
dev.off()
options(datadist = old_dd)
if (exists("dd_fix1", envir = .GlobalEnv, inherits = FALSE)) rm("dd_fix1", envir = .GlobalEnv)

## ---------------- Exploratory 3-year apparent grouped calibration ------------
## The 07A bootstrap calibration was unstable: 1-year had only 10 observed
## deaths by 1 year; 5-year resampling repeatedly produced sparse intervals.
## Rather than silently forcing those figures, provide one transparent 3-year
## apparent calibration diagnostic in four equal-sized predicted-risk groups.
grouped_calibration <- function(a, tt = 1095, groups = 4L) {
  d <- a$data
  risk <- predicted_risk(a$fit_comb, d, tt)
  cuts <- unique(quantile(risk, probs = seq(0, 1, length.out = groups + 1), na.rm = TRUE))
  if (length(cuts) < 3) stop("Insufficient unique predicted risks for grouped calibration")
  grp <- cut(risk, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  out <- do.call(rbind, lapply(sort(unique(grp)), function(g) {
    di <- d[grp == g, , drop = FALSE]
    sf <- survival::survfit(survival::Surv(OS_time, OS_event) ~ 1, data = di)
    sm <- summary(sf, times = tt, extend = TRUE)
    data.frame(
      Group = g, N = nrow(di), Events = sum(di$OS_event == 1),
      PredictedSurvival = mean(1 - risk[grp == g]),
      ObservedSurvival = as.numeric(sm$surv[1]),
      Lower95 = as.numeric(sm$lower[1]), Upper95 = as.numeric(sm$upper[1])
    )
  }))
  out
}

cal3 <- grouped_calibration(emtab, 1095, 4L)
write.csv(cal3, file.path(out_dir, "Tables", "10_EMTAB179_3year_apparent_grouped_calibration.csv"), row.names = FALSE)
pcal <- ggplot2::ggplot(cal3, ggplot2::aes(x = PredictedSurvival, y = ObservedSurvival)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower95, ymax = Upper95), width = 0.01, color = "#4C78A8") +
  ggplot2::geom_point(size = 2.8, color = "#4C78A8") +
  ggplot2::geom_line(color = "#4C78A8") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(title = "E-MTAB-179: 3-year apparent grouped calibration",
                subtitle = "Exploratory diagnostic; four predicted-risk groups",
                x = "Mean predicted 3-year overall survival",
                y = "Kaplan-Meier observed 3-year overall survival") +
  ggplot2::theme_bw(base_size = 11)
ggplot2::ggsave(file.path(out_dir, "Figures", "04_EMTAB179_3year_apparent_grouped_calibration.png"),
                pcal, width = 6, height = 6, dpi = 300, bg = "white")

calibration_status <- data.frame(
  TimeYear = c(1, 3, 5),
  ObservedDeathsByTime = c(
    sum(emtab$data$OS_event == 1 & emtab$data$OS_time <= 365),
    sum(emtab$data$OS_event == 1 & emtab$data$OS_time <= 1095),
    sum(emtab$data$OS_event == 1 & emtab$data$OS_time <= 1825)
  ),
  RECALC07AWarning = c(
    "Ran out of iterations in bootstrap calibration",
    "Ran out of iterations in bootstrap calibration",
    "153 sparse-interval warnings in bootstrap calibration"
  ),
  FIX1Recommendation = c(
    "Do not use 1-year bootstrap calibration as primary evidence",
    "Use only the transparent apparent grouped calibration as exploratory supplementary evidence",
    "Do not use the unstable 5-year bootstrap calibration figure"
  )
)
write.csv(calibration_status, file.path(out_dir, "Tables", "11_calibration_audit_and_recommendation.csv"), row.names = FALSE)

final_audit <- data.frame(
  Item = c(
    "Frozen M070 MD5", "M070 refitted", "Features reselected", "SEQC cutoff retuned",
    "E-MTAB age source unit", "E-MTAB age corrected unit", "DCA probability conversion inside dca",
    "07A bootstrap C-index result changed", "Old 07A calibration figures approved for manuscript",
    "FIX1 calibration role"
  ),
  Value = c(
    expected_md5, "FALSE", "FALSE", "FALSE", "days", "months", "FALSE",
    "FALSE", "FALSE", "Exploratory 3-year apparent grouped calibration only"
  )
)
write.csv(final_audit, file.path(out_dir, "Tables", "12_FINAL_FIX1_AUDIT_SUMMARY.csv"), row.names = FALSE)

write.csv(warning_log, file.path(out_dir, "Logs", "01_warnings_captured.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "Logs", "02_sessionInfo.txt"))
writeLines(c(
  "07A FIX1 completed.",
  paste0("Frozen M070 MD5 retained: ", expected_md5),
  "No M070 refitting, feature reselection, SEQC cutoff tuning, or new machine learning was performed.",
  "E-MTAB-179 age was corrected from day-scale source values to months using 365.25/12 days per month.",
  "DCA now receives already-computed event probabilities directly; as_probability is not used.",
  "The 1000-bootstrap C-index results from 07A are unchanged by the linear age-unit conversion.",
  "Unstable 07A 1/3/5-year rms bootstrap calibration figures are not approved for manuscript use.",
  "A 3-year apparent grouped calibration figure is provided only as an exploratory supplementary diagnostic."
), file.path(out_dir, "Logs", "03_run_summary.txt"))

cat("\n============================================================\n")
cat("07A FIX1 completed successfully.\n")
cat("Frozen M070 untouched: TRUE\n")
cat("E-MTAB age unit corrected: days -> months\n")
cat("DCA double probability conversion removed: TRUE\n")
cat("Upload the ENTIRE FIX1 output folder as ZIP for audit.\n")
cat("============================================================\n")
