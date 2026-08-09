rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================
## PRDX4_RECALC_06_LockM070_SEQC_External.R
##
## Purpose
##   1. Use the model-selection decision already made from E-MTAB-179 only:
##      M070 = SuperPC selector -> Random Survival Forest (RSF) learner.
##   2. Refit/freeze that single pipeline on all E-MTAB-179 samples.
##   3. Freeze the selected genes, scaling parameters, RSF model, risk-score
##      orientation/scale, and E-MTAB-derived risk-group cutoff BEFORE SEQC
##      is opened.
##   4. Only after the lock record has been written, open SEQC once and perform
##      external validation. SEQC does not choose the model or cutoff.
##
## Model-selection evidence from SC28_RECALC_05 (E-MTAB repeated CV):
##   - M061 All->RSF:      mean C-index 0.8545114
##   - M070 SuperPC->RSF:  mean C-index 0.8541959; within 1-SE top tier;
##                         PRDX4 50/50; mean feature Jaccard 0.9620744
##   - M050 SuperPC->Ridge mean C-index 0.8437471; outside 1-SE top tier.
##   M070 is therefore locked by a transparent stability-informed decision rule
##   applied after the repeated-CV audit:
##   among models within 1-SE of the best repeated-CV C-index, prefer the
##   pipeline with the strongest feature stability.
##
## IMPORTANT
##   - This script does NOT rerun the 101-model competition.
##   - This script does NOT use SEQC outcomes for model selection.
##   - PRDX4 is NOT appended or forced into the selected set.
##   - All warnings are captured in Logs/01_warnings_captured.csv.
## ============================================================

seed_main <- 20260808L
root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"

emtab_expr_file <- file.path(
  root_dir,
  "SC13A_MANUAL_V6_EMTAB179_ADF_TabFix", "Tables",
  "05_EMTAB179_expr_gene_symbol_matched.rds"
)
emtab_clin_file <- file.path(
  root_dir,
  "SC13A_MANUAL_V5_EMTAB179_SemicolonFix", "Tables",
  "04_EMTAB179_clin_matched_raw.csv"
)
seqc_ready_file <- file.path(
  root_dir,
  "01_Preprocess", "SEQC_GSE49711_ready_FIX2",
  "SEQC_external_validation_ready_input_FIX2.rds"
)

out_dir <- file.path(
  root_dir,
  "SC28_RECALC_06_LockM070_SEQC_External_20260808"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

## ---------------------------- package gate ----------------------------
required_packages <- c("survival", "randomForestSRC", "superpc", "timeROC")
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

## ---------------------------- warning capture ----------------------------
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

## ---------------------------- helpers ----------------------------
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

scale_from_training <- function(train_mat) {
  mu <- rowMeans(train_mat, na.rm = TRUE)
  sdv <- apply(train_mat, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  z <- sweep(sweep(train_mat, 1, mu, "-"), 1, sdv, "/")
  list(z = z, center = mu, scale = sdv)
}

apply_locked_scaling <- function(mat, center, scale) {
  genes <- names(center)
  if (!all(genes %in% rownames(mat))) {
    stop("Expression matrix lacks locked gene(s): ",
         paste(setdiff(genes, rownames(mat)), collapse = ", "))
  }
  z <- mat[genes, , drop = FALSE]
  z <- sweep(z, 1, center[genes], "-")
  z <- sweep(z, 1, scale[genes], "/")
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

select_superpc_top20 <- function(train_z, genes, time, event) {
  x <- t(train_z[genes, , drop = FALSE])
  dat <- list(
    x = t(x),
    y = time,
    censoring.status = event,
    featurenames = genes
  )
  fit <- with_warnlog(
    "SuperPC selector on full E-MTAB-179",
    superpc::superpc.train(dat, type = "survival")
  )
  fs <- fit$feature.scores
  if (is.matrix(fs)) fs <- fs[, 1]
  fs <- as.numeric(fs)
  if (length(fs) != length(genes)) stop("Unexpected SuperPC feature-score length")
  names(fs) <- genes
  selected <- head(names(sort(abs(fs), decreasing = TRUE)), min(20L, length(fs)))
  list(model = fit, scores = fs, selected = selected)
}

fit_rsf_locked <- function(train_z, genes, time, event, seed) {
  set.seed(seed)
  x <- t(train_z[genes, , drop = FALSE])
  d <- as.data.frame(x, check.names = FALSE)
  d$time <- time
  d$event <- event
  rsf_formula <- stats::as.formula(
    "Surv(time,event) ~ .",
    env = asNamespace("survival")
  )
  fit <- with_warnlog(
    "Final M070 RSF fit on full E-MTAB-179",
    randomForestSRC::rfsrc(
      rsf_formula,
      data = d,
      ntree = 1000,
      nodesize = 10,
      importance = TRUE
    )
  )
  list(model = fit, train_data = d)
}

predict_rsf_mortality <- function(fit, new_z, genes, context) {
  nd <- as.data.frame(t(new_z[genes, , drop = FALSE]), check.names = FALSE)
  pr <- with_warnlog(context, predict(fit, newdata = nd))
  score <- safe_num(pr$predicted)
  if (length(score) != nrow(nd)) stop("RSF prediction length mismatch in ", context)
  names(score) <- rownames(nd)
  score
}

determine_risk_orientation <- function(raw_score, time, event) {
  d <- data.frame(time = time, event = event, score = raw_score)
  fit <- with_warnlog(
    "Training-only risk orientation Cox model",
    survival::coxph(survival::Surv(time, event) ~ score, data = d)
  )
  b <- as.numeric(coef(fit)[1])
  if (!is.finite(b)) stop("Could not determine training risk-score orientation")
  if (b < 0) -1 else 1
}

calc_metrics_locked <- function(time, event, score, cutoff, dataset_label) {
  d <- data.frame(
    time = safe_num(time),
    event = safe_num(event),
    score = safe_num(score)
  )
  d <- d[complete.cases(d) & d$time > 0 & d$event %in% c(0, 1), , drop = FALSE]
  if (nrow(d) < 30 || sum(d$event == 1) < 5 || length(unique(d$score)) < 2) {
    return(data.frame(
      Dataset = dataset_label, N = nrow(d), Events = sum(d$event == 1),
      HR_per_trainingSD = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_,
      CoxP = NA_real_, Cindex = NA_real_, HighN = NA_integer_, LowN = NA_integer_,
      KM_HR_High_vs_Low = NA_real_, LogRankP = NA_real_,
      AUC1 = NA_real_, AUC3 = NA_real_, AUC5 = NA_real_, MeanAUC = NA_real_
    ))
  }

  ## score is already expressed in units of the E-MTAB training-score SD.
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

  d$group <- factor(ifelse(d$score > cutoff, "High", "Low"), levels = c("Low", "High"))
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
  logrank_p <- if (is.null(lr_fit)) NA_real_ else 1 - pchisq(lr_fit$chisq, length(lr_fit$n) - 1)

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
    HighN = sum(d$group == "High"),
    LowN = sum(d$group == "Low"),
    KM_HR_High_vs_Low = km_hr,
    LogRankP = logrank_p,
    AUC1 = auc[1],
    AUC3 = auc[2],
    AUC5 = auc[3],
    MeanAUC = if (all(is.na(auc))) NA_real_ else mean(auc, na.rm = TRUE)
  )
}

## ---------------------------- input checks ----------------------------
for (f in c(emtab_expr_file, emtab_clin_file, seqc_ready_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

## ---------------------------- E-MTAB only ----------------------------
## No SEQC object is opened above or in this section.
emtab_expr <- collapse_duplicate_genes(as_expr_matrix(readRDS(emtab_expr_file)))
emtab_clin <- read.csv(emtab_clin_file, check.names = FALSE)
if (!"Sample" %in% colnames(emtab_clin)) stop("E-MTAB clinical lacks Sample")

common_samples <- intersect(colnames(emtab_expr), as.character(emtab_clin$Sample))
if (length(common_samples) < 30) stop("Too few matched E-MTAB samples")
emtab_expr <- emtab_expr[, common_samples, drop = FALSE]
emtab_clin <- emtab_clin[match(common_samples, emtab_clin$Sample), , drop = FALSE]
emtab_time <- safe_num(emtab_clin$OVERALLSURVIVAL)
emtab_event <- clean_event(emtab_clin$DEATHOFDISEASE)
keep <- is.finite(emtab_time) & emtab_time > 0 & !is.na(emtab_event) & emtab_event %in% c(0, 1)
emtab_expr <- emtab_expr[, keep, drop = FALSE]
emtab_clin <- emtab_clin[keep, , drop = FALSE]
emtab_time <- emtab_time[keep]
emtab_event <- emtab_event[keep]
emtab_samples <- colnames(emtab_expr)

## Frozen before this step from the audited cross-platform panel.
candidate_genes <- c(
  "PRDX4", "SLC16A1", "LDHA", "LDHB", "PGK1", "ENO1", "PDK1", "HK2",
  "GPI", "ALDOA", "PRDX6", "PRDX2", "TXN", "PRDX5", "TXNRD1", "NQO1",
  "NUDT5", "ELK1", "SLC25A5", "SNAPC1", "NUP37", "RUVBL1", "SLC1A5", "AHCY"
)
if (length(candidate_genes) != 24L || anyDuplicated(candidate_genes)) {
  stop("Frozen 24-gene candidate panel integrity check failed")
}
missing_emtab <- setdiff(candidate_genes, rownames(emtab_expr))
if (length(missing_emtab)) {
  stop("Frozen candidate gene(s) missing from E-MTAB: ", paste(missing_emtab, collapse = ", "))
}
write.csv(
  data.frame(Gene = candidate_genes, FrozenBeforeFinalRefit = TRUE),
  file.path(out_dir, "Tables", "01_frozen_24gene_panel.csv"),
  row.names = FALSE
)

## Fit E-MTAB expression scaling on ALL training patients.
emtab_scale <- scale_from_training(emtab_expr[candidate_genes, , drop = FALSE])
emtab_z24 <- emtab_scale$z
write.csv(
  data.frame(
    Gene = candidate_genes,
    Center = emtab_scale$center[candidate_genes],
    Scale = emtab_scale$scale[candidate_genes]
  ),
  file.path(out_dir, "Tables", "02_locked_expression_scaling.csv"),
  row.names = FALSE
)

## M070 selector: SuperPC, identical top-20 rule to the audited repeated-CV run.
set.seed(seed_main + 70001L)
spc <- select_superpc_top20(
  train_z = emtab_z24,
  genes = candidate_genes,
  time = emtab_time,
  event = emtab_event
)
selected_genes <- spc$selected
if (length(selected_genes) != 20L) stop("M070 SuperPC did not return exactly 20 genes")

spc_table <- data.frame(
  Gene = candidate_genes,
  SuperPC_feature_score = as.numeric(spc$scores[candidate_genes]),
  Abs_feature_score = abs(as.numeric(spc$scores[candidate_genes])),
  SelectedTop20 = candidate_genes %in% selected_genes,
  SelectionRank = match(candidate_genes, names(sort(abs(spc$scores), decreasing = TRUE)))
)
spc_table <- spc_table[order(spc_table$SelectionRank), , drop = FALSE]
write.csv(
  spc_table,
  file.path(out_dir, "Tables", "03_SuperPC_feature_scores_and_locked20.csv"),
  row.names = FALSE
)

## M070 learner: RSF with the same hyperparameters as the repeated-CV audit.
rsf_seed <- seed_main + 70002L
rsf_obj <- fit_rsf_locked(
  train_z = emtab_z24,
  genes = selected_genes,
  time = emtab_time,
  event = emtab_event,
  seed = rsf_seed
)
rsf_fit <- rsf_obj$model

train_raw <- predict_rsf_mortality(
  rsf_fit,
  emtab_z24,
  selected_genes,
  "Final M070 RSF prediction on E-MTAB-179"
)
train_raw <- train_raw[emtab_samples]
orientation_sign <- determine_risk_orientation(train_raw, emtab_time, emtab_event)
train_oriented <- orientation_sign * train_raw
train_risk_mean <- mean(train_oriented, na.rm = TRUE)
train_risk_sd <- sd(train_oriented, na.rm = TRUE)
if (!is.finite(train_risk_mean) || !is.finite(train_risk_sd) || train_risk_sd <= 0) {
  stop("Invalid E-MTAB training risk-score distribution")
}
train_score <- (train_oriented - train_risk_mean) / train_risk_sd
locked_cutoff <- median(train_score, na.rm = TRUE)

## Save RSF variable importance. The 20 selected genes are not filtered again.
rsf_importance <- rsf_fit$importance
if (is.null(rsf_importance)) {
  rsf_imp_table <- data.frame(Gene = selected_genes, RSF_importance = NA_real_)
} else {
  rsf_importance <- as.numeric(rsf_importance)
  names(rsf_importance) <- names(rsf_fit$importance)
  rsf_imp_table <- data.frame(
    Gene = selected_genes,
    RSF_importance = as.numeric(rsf_importance[selected_genes])
  )
  rsf_imp_table <- rsf_imp_table[order(-rsf_imp_table$RSF_importance), , drop = FALSE]
}
write.csv(
  rsf_imp_table,
  file.path(out_dir, "Tables", "04_locked_M070_RSF_importance.csv"),
  row.names = FALSE
)

## ---------------------------- HARD LOCK BEFORE SEQC ----------------------------
## The following record and model objects are written before readRDS(seqc_ready_file).
lock_record <- data.frame(
  Item = c(
    "Lock timestamp", "Locked model ID", "Selector", "Learner",
    "Model-selection dataset", "Repeated-CV rule", "Candidate gene N",
    "Selected gene N", "Selected genes", "PRDX4 selected", "RSF ntree",
    "RSF nodesize", "Main seed", "RSF seed", "Risk orientation sign",
    "Training risk mean", "Training risk SD", "Locked risk cutoff",
    "SEQC opened before this lock", "SEQC outcomes used for model selection",
    "SEQC outcomes used for cutoff"
  ),
  Value = c(
    as.character(Sys.time()), "M070_SuperPC_RSF", "SuperPC", "RSF",
    "E-MTAB-179 only",
    "Among models within 1-SE of best repeated-CV mean C-index, prefer strongest feature stability",
    length(candidate_genes), length(selected_genes), paste(selected_genes, collapse = ";"),
    "PRDX4" %in% selected_genes, 1000, 10, seed_main, rsf_seed,
    orientation_sign, train_risk_mean, train_risk_sd, locked_cutoff,
    "FALSE", "FALSE", "FALSE"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  lock_record,
  file.path(out_dir, "Tables", "05_LOCK_RECORD_BEFORE_SEQC.csv"),
  row.names = FALSE
)

locked_model <- list(
  ModelID = "M070_SuperPC_RSF",
  Selector = "SuperPC",
  Learner = "RSF",
  CandidateGenes = candidate_genes,
  SelectedGenes = selected_genes,
  PRDX4Selected = "PRDX4" %in% selected_genes,
  SuperPCFeatureScores = spc$scores,
  ExpressionCenter = emtab_scale$center[candidate_genes],
  ExpressionScale = emtab_scale$scale[candidate_genes],
  RSFModel = rsf_fit,
  RSFNtree = 1000L,
  RSFNodesize = 10L,
  RiskOrientationSign = orientation_sign,
  TrainingRiskMean = train_risk_mean,
  TrainingRiskSD = train_risk_sd,
  LockedRiskCutoff = locked_cutoff,
  MainSeed = seed_main,
  RSFSeed = rsf_seed,
  SEQCUsedForModelSelection = FALSE,
  SEQCUsedForCutoff = FALSE
)
saveRDS(
  locked_model,
  file.path(out_dir, "Models", "M070_locked_model_BEFORE_SEQC.rds")
)
saveRDS(
  rsf_fit,
  file.path(out_dir, "Models", "M070_locked_RSF_model_BEFORE_SEQC.rds")
)

train_scores <- data.frame(
  Sample = emtab_samples,
  OS_time = emtab_time,
  OS_event = emtab_event,
  M070_RiskScore = as.numeric(train_score),
  RiskGroup_LockedCutoff = ifelse(train_score > locked_cutoff, "High", "Low")
)
write.csv(
  train_scores,
  file.path(out_dir, "Tables", "06_EMTAB_locked_M070_scores.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("M070 IS NOW LOCKED BEFORE SEQC.\n")
cat("Selected genes (", length(selected_genes), "): ", paste(selected_genes, collapse = ", "), "\n", sep = "")
cat("PRDX4 selected: ", "PRDX4" %in% selected_genes, "\n", sep = "")
cat("Locked training cutoff: ", signif(locked_cutoff, 7), "\n", sep = "")
cat("============================================================\n\n")

## ---------------------------- NOW OPEN SEQC ----------------------------
## From this line onward the complete M070 pipeline and cutoff are immutable.
seqc_obj <- readRDS(seqc_ready_file)
if (is.null(seqc_obj$expr) || is.null(seqc_obj$clin)) {
  stop("SEQC ready object must contain both $expr and $clin")
}
seqc_expr <- collapse_duplicate_genes(as_expr_matrix(seqc_obj$expr))
missing_seqc <- setdiff(candidate_genes, rownames(seqc_expr))
write.csv(
  data.frame(
    Gene = candidate_genes,
    In_SEQC = candidate_genes %in% rownames(seqc_expr),
    SelectedByLockedSuperPC = candidate_genes %in% selected_genes
  ),
  file.path(out_dir, "Tables", "07_SEQC_locked_gene_coverage.csv"),
  row.names = FALSE
)
if (length(missing_seqc)) {
  stop("Frozen candidate gene(s) unexpectedly missing from SEQC: ", paste(missing_seqc, collapse = ", "))
}

seqc_z24 <- apply_locked_scaling(
  seqc_expr,
  center = emtab_scale$center[candidate_genes],
  scale = emtab_scale$scale[candidate_genes]
)
seqc_raw <- predict_rsf_mortality(
  rsf_fit,
  seqc_z24,
  selected_genes,
  "Locked M070 RSF prediction in SEQC"
)
seqc_oriented <- orientation_sign * seqc_raw
seqc_score <- (seqc_oriented - train_risk_mean) / train_risk_sd

## Only now access SEQC survival outcomes.
seqc_clin <- as.data.frame(seqc_obj$clin)
if (!"Sample" %in% colnames(seqc_clin)) seqc_clin$Sample <- rownames(seqc_clin)
seqc_clin$Sample <- as.character(seqc_clin$Sample)
te <- find_time_event_cols(seqc_clin)
if (is.na(te$time) || is.na(te$event)) {
  stop("Cannot identify SEQC survival time/event columns")
}
seqc_time_all <- safe_num(seqc_clin[[te$time]])
seqc_event_all <- clean_event(seqc_clin[[te$event]])

seqc_pred <- data.frame(
  Sample = names(seqc_score),
  M070_RiskScore = as.numeric(seqc_score),
  stringsAsFactors = FALSE
)
idx <- match(seqc_pred$Sample, seqc_clin$Sample)
seqc_pred$OS_time <- seqc_time_all[idx]
seqc_pred$OS_event <- seqc_event_all[idx]
seqc_pred <- seqc_pred[
  complete.cases(seqc_pred[, c("OS_time", "OS_event", "M070_RiskScore")]) &
    seqc_pred$OS_time > 0 & seqc_pred$OS_event %in% c(0, 1),
  , drop = FALSE
]
if (nrow(seqc_pred) < 30) stop("Too few SEQC samples after matching survival data")
seqc_pred$RiskGroup_LockedCutoff <- ifelse(
  seqc_pred$M070_RiskScore > locked_cutoff,
  "High", "Low"
)
seqc_pred <- seqc_pred[, c(
  "Sample", "OS_time", "OS_event", "M070_RiskScore", "RiskGroup_LockedCutoff"
)]
write.csv(
  seqc_pred,
  file.path(out_dir, "Tables", "08_SEQC_EXTERNAL_locked_M070_scores.csv"),
  row.names = FALSE
)

## ---------------------------- performance ----------------------------
perf <- rbind(
  calc_metrics_locked(
    train_scores$OS_time,
    train_scores$OS_event,
    train_scores$M070_RiskScore,
    cutoff = locked_cutoff,
    dataset_label = "E-MTAB-179 full refit (apparent)"
  ),
  calc_metrics_locked(
    seqc_pred$OS_time,
    seqc_pred$OS_event,
    seqc_pred$M070_RiskScore,
    cutoff = locked_cutoff,
    dataset_label = "SEQC external validation (locked model and cutoff)"
  )
)
write.csv(
  perf,
  file.path(out_dir, "Tables", "09_M070_performance_EMTAB_and_SEQC_external.csv"),
  row.names = FALSE
)

final_summary <- data.frame(
  Item = c(
    "Locked model", "Lock based on dataset", "SEQC role", "Candidate genes",
    "Selected genes", "PRDX4 selected", "E-MTAB N", "E-MTAB events",
    "SEQC external N", "SEQC external events", "Locked cutoff source",
    "SEQC survival time column", "SEQC survival event column",
    "Warnings captured"
  ),
  Value = c(
    "M070_SuperPC_RSF", "E-MTAB-179 repeated CV only",
    "External validation only; not used for model selection or cutoff",
    length(candidate_genes), length(selected_genes), "PRDX4" %in% selected_genes,
    nrow(train_scores), sum(train_scores$OS_event == 1),
    nrow(seqc_pred), sum(seqc_pred$OS_event == 1),
    "Median of full-refit E-MTAB M070 training score",
    te$time, te$event, nrow(warning_log)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  final_summary,
  file.path(out_dir, "Tables", "10_FINAL_M070_EXTERNAL_VALIDATION_SUMMARY.csv"),
  row.names = FALSE
)

write.csv(
  warning_log,
  file.path(out_dir, "Logs", "01_warnings_captured.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "Logs", "02_sessionInfo.txt")
)

summary_lines <- c(
  "SC28 RECALC 06 completed.",
  paste0("End: ", as.character(Sys.time())),
  "Locked model: M070 = SuperPC -> RSF",
  paste0("Selected genes: ", paste(selected_genes, collapse = ", ")),
  paste0("PRDX4 selected: ", "PRDX4" %in% selected_genes),
  paste0("E-MTAB N/events: ", nrow(train_scores), "/", sum(train_scores$OS_event == 1)),
  paste0("SEQC external N/events: ", nrow(seqc_pred), "/", sum(seqc_pred$OS_event == 1)),
  paste0("Locked cutoff: ", signif(locked_cutoff, 8)),
  paste0("Warnings captured: ", nrow(warning_log)),
  "SEQC was not used for model selection or cutoff.",
  "Primary results: Tables/09_M070_performance_EMTAB_and_SEQC_external.csv",
  "Lock evidence: Tables/05_LOCK_RECORD_BEFORE_SEQC.csv"
)
writeLines(summary_lines, file.path(out_dir, "Logs", "03_run_summary.txt"))

cat("\n============================================================\n")
cat("SC28 RECALC 06 completed successfully.\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("============================================================\n")
