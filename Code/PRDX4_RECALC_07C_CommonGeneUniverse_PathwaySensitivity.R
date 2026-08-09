rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260809)

## ============================================================================
## PRDX4_RECALC_07C_CommonGeneUniverse_PathwaySensitivity.R
##
## Purpose
##   Strict cross-platform sensitivity analysis for the pathway layer of the
##   locked M070 = SuperPC -> RSF model.
##
## Why this analysis is run
##   In 07B FIX3, Hallmark/KEGG gene sets were intersected with each platform
##   separately. That is acceptable for within-cohort ssGSEA, but the same
##   nominal pathway can contain different measurable genes on E-MTAB-179 and
##   SEQC. 07C removes this platform-composition difference.
##
## Frozen analysis boundary
##   1. M070 is NOT refitted.
##   2. No model genes/features are reselected.
##   3. The locked E-MTAB-179 cutoff is NOT retuned.
##   4. SEQC survival outcomes are NOT used.
##   5. Exact MSigDB membership is read from the completed 07B FIX3 output;
##      no new MSigDB download or gene-set update occurs here.
##   6. Both cohorts are restricted to the identical common expression-gene
##      universe before ssGSEA.
##   7. Each Hallmark/KEGG pathway uses the identical shared gene membership in
##      both cohorts (minimum 5 genes, maximum 500 genes).
##   8. Results are sensitivity/association evidence, not causal proof.
## ============================================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
fix3_dir <- file.path(
  root_dir,
  "SC28_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism_20260809"
)
out_dir <- file.path(
  root_dir,
  "SC28_RECALC_07C_CommonGeneUniverse_PathwaySensitivity_20260809"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

## ----------------------------- fixed inputs ----------------------------------
emtab_expr_file <- file.path(
  root_dir,
  "SC13A_MANUAL_V6_EMTAB179_ADF_TabFix", "Tables",
  "05_EMTAB179_expr_gene_symbol_matched.rds"
)
seqc_ready_file <- file.path(
  root_dir,
  "01_Preprocess", "SEQC_GSE49711_ready_FIX2",
  "SEQC_external_validation_ready_input_FIX2.rds"
)
emtab_score_file <- file.path(
  root_dir,
  "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Tables",
  "06_EMTAB_locked_M070_scores.csv"
)
seqc_score_file <- file.path(
  root_dir,
  "SC28_RECALC_06B_M070_SEQC_CohortZ_External_20260808", "Tables",
  "06_SEQC_EXTERNAL_M070_scores_with_outcomes.csv"
)
locked_model_file <- file.path(
  root_dir,
  "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Models",
  "M070_locked_model_BEFORE_SEQC.rds"
)

fix3_membership_file <- file.path(
  fix3_dir, "Tables", "04_MSigDB_gene_set_membership_USED.csv"
)
fix3_hall_group_file <- file.path(
  fix3_dir, "Tables", "10_Hallmark_cross_cohort_locked_group_concordance.csv"
)
fix3_hall_cont_file <- file.path(
  fix3_dir, "Tables", "11_Hallmark_cross_cohort_continuous_concordance.csv"
)
fix3_kegg_group_file <- file.path(
  fix3_dir, "Tables", "17_KEGG_cross_cohort_locked_group_concordance.csv"
)
fix3_kegg_cont_file <- file.path(
  fix3_dir, "Tables", "18_KEGG_cross_cohort_continuous_concordance.csv"
)

required_files <- c(
  emtab_expr_file, seqc_ready_file, emtab_score_file, seqc_score_file,
  locked_model_file, fix3_membership_file, fix3_hall_group_file,
  fix3_hall_cont_file, fix3_kegg_group_file, fix3_kegg_cont_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\n\n07C requires the completed 07B FIX3 output folder at the exact path shown above."
  )
}

## ----------------------------- packages --------------------------------------
required_packages <- c("GSVA", "limma", "ggplot2")
pkg_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
write.csv(
  data.frame(Package = required_packages, Available = unname(pkg_ok)),
  file.path(out_dir, "Tables", "00_package_status.csv"),
  row.names = FALSE
)
if (any(!pkg_ok)) {
  stop(
    "Missing R package(s): ",
    paste(required_packages[!pkg_ok], collapse = ", "),
    ". These packages were already used by 07B FIX3; use the same R library and rerun."
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

as_expr_matrix <- function(x) {
  if (is.matrix(x)) {
    m <- x
  } else if (inherits(x, "Matrix")) {
    m <- as.matrix(x)
  } else if (is.data.frame(x)) {
    m <- as.matrix(x)
  } else {
    stop("Unsupported expression object class: ", paste(class(x), collapse = "/"))
  }
  storage.mode(m) <- "double"
  if (is.null(rownames(m)) || is.null(colnames(m))) {
    stop("Expression matrix must have gene row names and sample column names")
  }
  if (any(!is.finite(m))) stop("Expression matrix contains non-finite values")
  m
}

collapse_duplicate_genes <- function(mat) {
  rn <- toupper(trimws(rownames(mat)))
  keep <- nzchar(rn) & !is.na(rn)
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  if (!anyDuplicated(rn)) {
    rownames(mat) <- rn
    return(mat)
  }
  summed <- rowsum(mat, group = rn, reorder = FALSE)
  counts <- table(rn)[rownames(summed)]
  out <- summed / as.numeric(counts)
  colnames(out) <- colnames(mat)
  out
}

align_expression_to_scores <- function(
    expr, scores, score_col, group_col, dataset_label, expected_n) {
  needed <- c("Sample", score_col, group_col)
  if (!all(needed %in% names(scores))) {
    stop(dataset_label, " score table lacks required columns")
  }
  scores <- scores[, needed, drop = FALSE]
  scores$Sample <- as.character(scores$Sample)
  scores$M070 <- safe_num(scores[[score_col]])
  scores$RiskGroup <- factor(as.character(scores[[group_col]]), levels = c("Low", "High"))
  if (any(!is.finite(scores$M070)) || any(is.na(scores$RiskGroup))) {
    stop(dataset_label, " contains invalid M070 score or locked risk group")
  }

  expr_names <- colnames(expr)
  idx <- match(scores$Sample, expr_names)
  if (anyNA(idx)) {
    ek <- sample_key(expr_names)
    sk <- sample_key(scores$Sample)
    if (anyDuplicated(ek)) stop(dataset_label, " expression sample keys are duplicated")
    if (anyDuplicated(sk)) stop(dataset_label, " score sample keys are duplicated")
    idx <- match(sk, ek)
  }
  if (anyNA(idx)) {
    stop(dataset_label, " score-to-expression sample matching failed")
  }
  expr2 <- expr[, idx, drop = FALSE]
  colnames(expr2) <- scores$Sample
  if (ncol(expr2) != expected_n || nrow(scores) != expected_n) {
    stop(dataset_label, " unexpected sample count after alignment")
  }
  list(expr = expr2, scores = scores)
}

## -------------------------- frozen M070 checks -------------------------------
frozen_m070 <- readRDS(locked_model_file)
required_frozen_fields <- c(
  "ModelID", "Selector", "Learner", "SelectedGenes", "LockedRiskCutoff",
  "SEQCUsedForModelSelection", "SEQCUsedForCutoff"
)
missing_fields <- setdiff(required_frozen_fields, names(frozen_m070))
if (length(missing_fields)) {
  stop("Frozen M070 RDS lacks required field(s): ", paste(missing_fields, collapse = ", "))
}

model_id <- as.character(frozen_m070$ModelID)[1]
selector <- as.character(frozen_m070$Selector)[1]
learner <- as.character(frozen_m070$Learner)[1]
selected_genes <- as.character(frozen_m070$SelectedGenes)
lock_cutoff <- safe_num(frozen_m070$LockedRiskCutoff)[1]
model_md5 <- unname(tools::md5sum(locked_model_file))
expected_md5 <- "e719d0e6b626e75cc134e83fd0929083"

if (!identical(model_id, "M070_SuperPC_RSF")) stop("Unexpected frozen model ID")
if (!identical(selector, "SuperPC") || !identical(learner, "RSF")) {
  stop("Frozen M070 is not SuperPC -> RSF")
}
if (!"PRDX4" %in% selected_genes) stop("PRDX4 is absent from frozen M070 selected genes")
if (!is.finite(lock_cutoff)) stop("Frozen M070 cutoff is not finite")
if (!identical(frozen_m070$SEQCUsedForModelSelection, FALSE)) {
  stop("Frozen provenance indicates SEQC was used for model selection")
}
if (!identical(frozen_m070$SEQCUsedForCutoff, FALSE)) {
  stop("Frozen provenance indicates SEQC was used for cutoff selection")
}
if (!identical(model_md5, expected_md5)) {
  stop("Frozen M070 MD5 mismatch. Expected ", expected_md5, "; observed ", model_md5)
}

## -------------------------- load/align expression ----------------------------
emtab_expr <- collapse_duplicate_genes(as_expr_matrix(readRDS(emtab_expr_file)))
seqc_obj <- readRDS(seqc_ready_file)
if (is.null(seqc_obj$expr)) stop("SEQC ready object lacks $expr")
seqc_expr <- collapse_duplicate_genes(as_expr_matrix(seqc_obj$expr))

emtab <- align_expression_to_scores(
  emtab_expr, read.csv(emtab_score_file, check.names = FALSE),
  score_col = "M070_RiskScore",
  group_col = "RiskGroup_LockedCutoff",
  dataset_label = "E-MTAB-179", expected_n = 478L
)
seqc <- align_expression_to_scores(
  seqc_expr, read.csv(seqc_score_file, check.names = FALSE),
  score_col = "Corrected06B_M070_RiskScore",
  group_col = "RiskGroup_Locked_EMTAB_Cutoff",
  dataset_label = "SEQC", expected_n = 498L
)

emtab_expected_group <- ifelse(emtab$scores$M070 >= lock_cutoff, "High", "Low")
seqc_expected_group <- ifelse(seqc$scores$M070 >= lock_cutoff, "High", "Low")
if (!all(emtab_expected_group == as.character(emtab$scores$RiskGroup))) {
  stop("E-MTAB-179 risk groups do not match frozen cutoff")
}
if (!all(seqc_expected_group == as.character(seqc$scores$RiskGroup))) {
  stop("SEQC risk groups do not match frozen E-MTAB-179 cutoff")
}

## Exact common expression universe. Both matrices are ordered identically.
common_genes <- sort(intersect(rownames(emtab$expr), rownames(seqc$expr)))
if (length(common_genes) < 1000L) {
  stop("Unexpectedly small common expression-gene universe: ", length(common_genes))
}
emtab_common <- emtab$expr[common_genes, , drop = FALSE]
seqc_common <- seqc$expr[common_genes, , drop = FALSE]
if (!identical(rownames(emtab_common), rownames(seqc_common))) {
  stop("Common expression-gene universe ordering mismatch")
}

input_audit <- data.frame(
  Item = c(
    "Frozen model ID", "Selector", "Learner", "Frozen M070 MD5",
    "Expected frozen MD5", "Locked cutoff", "M070 refitted in 07C",
    "Model features reselected in 07C", "Cutoff retuned in 07C",
    "SEQC survival outcomes used in 07C", "E-MTAB-179 original genes",
    "SEQC original genes", "Identical common expression-gene universe",
    "E-MTAB-179 N", "SEQC N", "FIX3 MSigDB membership reused"
  ),
  Value = c(
    model_id, selector, learner, model_md5, expected_md5, lock_cutoff,
    "FALSE", "FALSE", "FALSE", "FALSE", nrow(emtab$expr), nrow(seqc$expr),
    length(common_genes), ncol(emtab_common), ncol(seqc_common), "TRUE"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  input_audit,
  file.path(out_dir, "Tables", "01_design_and_common_universe_audit.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Gene = common_genes),
  file.path(out_dir, "Tables", "02_identical_common_expression_gene_universe.csv"),
  row.names = FALSE
)

## ---------------------- reuse exact FIX3 MSigDB sets -------------------------
membership <- read.csv(fix3_membership_file, check.names = FALSE)
needed_membership_cols <- c("Collection", "GeneSet", "Gene")
if (!all(needed_membership_cols %in% names(membership))) {
  stop("FIX3 MSigDB membership table lacks required columns")
}
membership$Gene <- toupper(trimws(as.character(membership$Gene)))
membership$GeneSet <- as.character(membership$GeneSet)
membership$Collection <- as.character(membership$Collection)

hall_rows <- membership$Collection == "Hallmark"
kegg_rows <- grepl("^KEGG_", membership$Collection)
if (!any(hall_rows) || !any(kegg_rows)) {
  stop("Could not reconstruct Hallmark/KEGG sets from FIX3 membership table")
}

hallmark_sets_original <- split(membership$Gene[hall_rows], membership$GeneSet[hall_rows])
kegg_sets_original <- split(membership$Gene[kegg_rows], membership$GeneSet[kegg_rows])
hallmark_sets_original <- lapply(hallmark_sets_original, unique)
kegg_sets_original <- lapply(kegg_sets_original, unique)

build_shared_sets <- function(sets, common_universe, family_label, min_size = 5L, max_size = 500L) {
  rows <- lapply(names(sets), function(nm) {
    planned <- unique(toupper(sets[[nm]]))
    shared <- intersect(planned, common_universe)
    data.frame(
      AnalysisFamily = family_label,
      GeneSet = nm,
      FIX3MembershipN = length(planned),
      SharedCommonUniverseN = length(shared),
      CoverageFraction = ifelse(length(planned) > 0, length(shared) / length(planned), NA_real_),
      Evaluable = length(shared) >= min_size && length(shared) <= max_size,
      SharedGenes = paste(shared, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  coverage <- do.call(rbind, rows)
  shared_sets <- setNames(
    lapply(seq_len(nrow(coverage)), function(i) {
      z <- coverage$SharedGenes[i]
      if (!nzchar(z)) character() else strsplit(z, ";", fixed = TRUE)[[1]]
    }),
    coverage$GeneSet
  )
  shared_sets <- shared_sets[coverage$Evaluable]
  list(sets = shared_sets, coverage = coverage)
}

hall_shared <- build_shared_sets(hallmark_sets_original, common_genes, "Hallmark")
kegg_shared <- build_shared_sets(kegg_sets_original, common_genes, "KEGG")
if (!length(hall_shared$sets)) stop("No Hallmark set survived common-universe filtering")
if (!length(kegg_shared$sets)) stop("No KEGG set survived common-universe filtering")

coverage_all <- rbind(hall_shared$coverage, kegg_shared$coverage)
write.csv(
  coverage_all,
  file.path(out_dir, "Tables", "03_pathway_common_gene_coverage.csv"),
  row.names = FALSE
)

shared_gene_long <- function(sets, family_label) {
  do.call(rbind, lapply(names(sets), function(nm) {
    data.frame(
      AnalysisFamily = family_label,
      GeneSet = nm,
      Gene = sets[[nm]],
      stringsAsFactors = FALSE
    )
  }))
}
write.csv(
  rbind(
    shared_gene_long(hall_shared$sets, "Hallmark"),
    shared_gene_long(kegg_shared$sets, "KEGG")
  ),
  file.path(out_dir, "Tables", "04_IDENTICAL_pathway_gene_membership_used_in_both_cohorts.csv"),
  row.names = FALSE
)

## ----------------------------- ssGSEA ----------------------------------------
run_ssgsea_fixed <- function(expr_common, fixed_sets, context) {
  if (!all(unlist(fixed_sets, use.names = FALSE) %in% rownames(expr_common))) {
    stop(context, ": fixed gene-set member absent from common expression universe")
  }
  ex <- getNamespaceExports("GSVA")
  if ("ssgseaParam" %in% ex) {
    par <- with_warnlog(
      paste(context, "ssgseaParam"),
      GSVA::ssgseaParam(
        expr_common, fixed_sets,
        minSize = 5, maxSize = 500,
        alpha = 0.25, normalize = TRUE,
        checkNA = "yes", use = "all.obs",
        verbose = FALSE
      )
    )
    es <- with_warnlog(
      paste(context, "GSVA::gsva"),
      GSVA::gsva(par, verbose = FALSE)
    )
  } else {
    es <- with_warnlog(
      paste(context, "legacy GSVA::gsva"),
      GSVA::gsva(
        expr_common, fixed_sets, method = "ssgsea",
        min.sz = 5, max.sz = 500,
        ssgsea.norm = TRUE,
        parallel.sz = 1, verbose = FALSE
      )
    )
  }
  es <- as.matrix(es)
  if (any(!is.finite(es))) stop(context, ": non-finite ssGSEA score produced")
  es
}

limma_group_test <- function(score_mat, score_df, dataset, family_label) {
  rg <- factor(score_df$RiskGroup, levels = c("Low", "High"))
  design <- model.matrix(~ 0 + rg)
  colnames(design) <- c("Low", "High")
  cont <- limma::makeContrasts(HighMinusLow = High - Low, levels = design)
  fit <- limma::lmFit(score_mat, design)
  fit <- limma::contrasts.fit(fit, cont)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, number = Inf, sort.by = "none", adjust.method = "BH")
  tt$GeneSet <- rownames(tt)
  rownames(tt) <- NULL
  data.frame(
    Dataset = dataset, AnalysisFamily = family_label, GeneSet = tt$GeneSet,
    MeanDifference_HighMinusLow = tt$logFC, P = tt$P.Value, FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

continuous_spearman <- function(score_mat, score_df, dataset, family_label) {
  out <- lapply(seq_len(nrow(score_mat)), function(i) {
    ct <- suppressWarnings(cor.test(
      safe_num(score_mat[i, ]), score_df$M070,
      method = "spearman", exact = FALSE
    ))
    data.frame(
      Dataset = dataset, AnalysisFamily = family_label,
      GeneSet = rownames(score_mat)[i], Rho = unname(ct$estimate), P = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$FDR <- p.adjust(out$P, method = "BH")
  out
}

cross_group <- function(a, b, family_label) {
  x <- merge(
    a[, c("GeneSet", "MeanDifference_HighMinusLow", "P", "FDR")],
    b[, c("GeneSet", "MeanDifference_HighMinusLow", "P", "FDR")],
    by = "GeneSet", suffixes = c("_EMTAB", "_SEQC")
  )
  x$DirectionConcordant <- sign(x$MeanDifference_HighMinusLow_EMTAB) == sign(x$MeanDifference_HighMinusLow_SEQC)
  x$FDRlt005_Both <- x$FDR_EMTAB < 0.05 & x$FDR_SEQC < 0.05
  x$RobustBoth <- x$DirectionConcordant & x$FDRlt005_Both
  x$AnalysisFamily <- family_label
  x
}

cross_cont <- function(a, b, family_label) {
  x <- merge(
    a[, c("GeneSet", "Rho", "P", "FDR")],
    b[, c("GeneSet", "Rho", "P", "FDR")],
    by = "GeneSet", suffixes = c("_EMTAB", "_SEQC")
  )
  x$DirectionConcordant <- sign(x$Rho_EMTAB) == sign(x$Rho_SEQC)
  x$FDRlt005_Both <- x$FDR_EMTAB < 0.05 & x$FDR_SEQC < 0.05
  x$RobustBoth <- x$DirectionConcordant & x$FDRlt005_Both
  x$AnalysisFamily <- family_label
  x
}

write_score_matrix <- function(mat, path) {
  write.csv(
    data.frame(GeneSet = rownames(mat), mat, check.names = FALSE),
    path, row.names = FALSE
  )
}

hall_emtab <- run_ssgsea_fixed(emtab_common, hall_shared$sets, "E-MTAB-179 Hallmark common-universe")
hall_seqc <- run_ssgsea_fixed(seqc_common, hall_shared$sets, "SEQC Hallmark common-universe")
kegg_emtab <- run_ssgsea_fixed(emtab_common, kegg_shared$sets, "E-MTAB-179 KEGG common-universe")
kegg_seqc <- run_ssgsea_fixed(seqc_common, kegg_shared$sets, "SEQC KEGG common-universe")

write_score_matrix(hall_emtab, file.path(out_dir, "Tables", "05_Hallmark_commonUniverse_EMTAB_scores.csv"))
write_score_matrix(hall_seqc, file.path(out_dir, "Tables", "06_Hallmark_commonUniverse_SEQC_scores.csv"))
write_score_matrix(kegg_emtab, file.path(out_dir, "Tables", "07_KEGG_commonUniverse_EMTAB_scores.csv"))
write_score_matrix(kegg_seqc, file.path(out_dir, "Tables", "08_KEGG_commonUniverse_SEQC_scores.csv"))

hall_grp_em <- limma_group_test(hall_emtab, emtab$scores, "E-MTAB-179", "Hallmark")
hall_grp_sq <- limma_group_test(hall_seqc, seqc$scores, "SEQC", "Hallmark")
hall_cor_em <- continuous_spearman(hall_emtab, emtab$scores, "E-MTAB-179", "Hallmark")
hall_cor_sq <- continuous_spearman(hall_seqc, seqc$scores, "SEQC", "Hallmark")
kegg_grp_em <- limma_group_test(kegg_emtab, emtab$scores, "E-MTAB-179", "KEGG")
kegg_grp_sq <- limma_group_test(kegg_seqc, seqc$scores, "SEQC", "KEGG")
kegg_cor_em <- continuous_spearman(kegg_emtab, emtab$scores, "E-MTAB-179", "KEGG")
kegg_cor_sq <- continuous_spearman(kegg_seqc, seqc$scores, "SEQC", "KEGG")

hall_grp_cross <- cross_group(hall_grp_em, hall_grp_sq, "Hallmark")
hall_cor_cross <- cross_cont(hall_cor_em, hall_cor_sq, "Hallmark")
kegg_grp_cross <- cross_group(kegg_grp_em, kegg_grp_sq, "KEGG")
kegg_cor_cross <- cross_cont(kegg_cor_em, kegg_cor_sq, "KEGG")

write.csv(hall_grp_cross, file.path(out_dir, "Tables", "09_Hallmark_commonUniverse_group_concordance.csv"), row.names = FALSE)
write.csv(hall_cor_cross, file.path(out_dir, "Tables", "10_Hallmark_commonUniverse_continuous_concordance.csv"), row.names = FALSE)
write.csv(kegg_grp_cross, file.path(out_dir, "Tables", "11_KEGG_commonUniverse_group_concordance.csv"), row.names = FALSE)
write.csv(kegg_cor_cross, file.path(out_dir, "Tables", "12_KEGG_commonUniverse_continuous_concordance.csv"), row.names = FALSE)

## ---------------------- compare with 07B FIX3 --------------------------------
old_hall_group <- read.csv(fix3_hall_group_file, check.names = FALSE)
old_hall_cont <- read.csv(fix3_hall_cont_file, check.names = FALSE)
old_kegg_group <- read.csv(fix3_kegg_group_file, check.names = FALSE)
old_kegg_cont <- read.csv(fix3_kegg_cont_file, check.names = FALSE)

compare_group <- function(old, new, family_label) {
  old2 <- old[, c(
    "GeneSet", "MeanDifference_HighMinusLow_EMTAB", "FDR_EMTAB",
    "MeanDifference_HighMinusLow_SEQC", "FDR_SEQC",
    "DirectionConcordant", "FDRlt005_Both"
  )]
  new2 <- new[, c(
    "GeneSet", "MeanDifference_HighMinusLow_EMTAB", "FDR_EMTAB",
    "MeanDifference_HighMinusLow_SEQC", "FDR_SEQC",
    "DirectionConcordant", "FDRlt005_Both"
  )]
  x <- merge(old2, new2, by = "GeneSet", suffixes = c("_FIX3", "_07C"))
  x$AnalysisFamily <- family_label
  x$DirectionStable_EMTAB <- sign(x$MeanDifference_HighMinusLow_EMTAB_FIX3) == sign(x$MeanDifference_HighMinusLow_EMTAB_07C)
  x$DirectionStable_SEQC <- sign(x$MeanDifference_HighMinusLow_SEQC_FIX3) == sign(x$MeanDifference_HighMinusLow_SEQC_07C)
  x$RobustBoth_FIX3 <- x$DirectionConcordant_FIX3 & x$FDRlt005_Both_FIX3
  x$RobustBoth_07C <- x$DirectionConcordant_07C & x$FDRlt005_Both_07C
  x$RobustStatusPreserved <- x$RobustBoth_FIX3 == x$RobustBoth_07C
  x
}

compare_cont <- function(old, new, family_label) {
  old2 <- old[, c(
    "GeneSet", "Rho_EMTAB", "FDR_EMTAB", "Rho_SEQC", "FDR_SEQC",
    "DirectionConcordant", "FDRlt005_Both"
  )]
  new2 <- new[, c(
    "GeneSet", "Rho_EMTAB", "FDR_EMTAB", "Rho_SEQC", "FDR_SEQC",
    "DirectionConcordant", "FDRlt005_Both"
  )]
  x <- merge(old2, new2, by = "GeneSet", suffixes = c("_FIX3", "_07C"))
  x$AnalysisFamily <- family_label
  x$DirectionStable_EMTAB <- sign(x$Rho_EMTAB_FIX3) == sign(x$Rho_EMTAB_07C)
  x$DirectionStable_SEQC <- sign(x$Rho_SEQC_FIX3) == sign(x$Rho_SEQC_07C)
  x$RobustBoth_FIX3 <- x$DirectionConcordant_FIX3 & x$FDRlt005_Both_FIX3
  x$RobustBoth_07C <- x$DirectionConcordant_07C & x$FDRlt005_Both_07C
  x$RobustStatusPreserved <- x$RobustBoth_FIX3 == x$RobustBoth_07C
  x
}

hall_group_compare <- compare_group(old_hall_group, hall_grp_cross, "Hallmark")
hall_cont_compare <- compare_cont(old_hall_cont, hall_cor_cross, "Hallmark")
kegg_group_compare <- compare_group(old_kegg_group, kegg_grp_cross, "KEGG")
kegg_cont_compare <- compare_cont(old_kegg_cont, kegg_cor_cross, "KEGG")

write.csv(hall_group_compare, file.path(out_dir, "Tables", "13_Hallmark_FIX3_vs_07C_group_sensitivity.csv"), row.names = FALSE)
write.csv(hall_cont_compare, file.path(out_dir, "Tables", "14_Hallmark_FIX3_vs_07C_continuous_sensitivity.csv"), row.names = FALSE)
write.csv(kegg_group_compare, file.path(out_dir, "Tables", "15_KEGG_FIX3_vs_07C_group_sensitivity.csv"), row.names = FALSE)
write.csv(kegg_cont_compare, file.path(out_dir, "Tables", "16_KEGG_FIX3_vs_07C_continuous_sensitivity.csv"), row.names = FALSE)

## ----------------------------- core pathways ---------------------------------
core_hallmark <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)
core_kegg <- c(
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_CITRATE_CYCLE_TCA_CYCLE",
  "KEGG_PYRUVATE_METABOLISM",
  "KEGG_CELL_CYCLE",
  "KEGG_DNA_REPLICATION",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
  "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"
)

core_hall_result <- hall_cont_compare[hall_cont_compare$GeneSet %in% core_hallmark, , drop = FALSE]
core_kegg_result <- kegg_cont_compare[kegg_cont_compare$GeneSet %in% core_kegg, , drop = FALSE]
core_hall_result$CoreFamily <- "Hallmark"
core_kegg_result$CoreFamily <- "KEGG"
if (!identical(names(core_hall_result), names(core_kegg_result))) {
  stop("Unexpected Hallmark/KEGG comparison schema mismatch")
}
core_all <- rbind(core_hall_result, core_kegg_result)
write.csv(
  core_all,
  file.path(out_dir, "Tables", "17_CORE_PATHWAY_SENSITIVITY_AUDIT.csv"),
  row.names = FALSE
)

## ------------------------------- figures -------------------------------------
theme_pub <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    strip.text = ggplot2::element_text(face = "bold")
  )

make_rho_stability_plot <- function(comp, title, output_file) {
  d1 <- data.frame(
    Dataset = "E-MTAB-179", GeneSet = comp$GeneSet,
    FIX3 = comp$Rho_EMTAB_FIX3, CommonUniverse = comp$Rho_EMTAB_07C
  )
  d2 <- data.frame(
    Dataset = "SEQC", GeneSet = comp$GeneSet,
    FIX3 = comp$Rho_SEQC_FIX3, CommonUniverse = comp$Rho_SEQC_07C
  )
  d <- rbind(d1, d2)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = FIX3, y = CommonUniverse)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey75") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey75") +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#2166AC", linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8, alpha = 0.7, color = "#333333") +
    ggplot2::facet_wrap(~Dataset, nrow = 1) +
    ggplot2::labs(
      title = title,
      subtitle = "Continuous M070 Spearman correlations",
      x = "Original FIX3 rho", y = "Common-gene-universe sensitivity rho"
    ) + theme_pub
  ggplot2::ggsave(output_file, p, width = 9.5, height = 4.7, dpi = 320, bg = "white")
}

make_rho_stability_plot(
  hall_cont_compare,
  "Hallmark pathway stability after common-gene-universe harmonization",
  file.path(out_dir, "Figures", "01_Hallmark_commonUniverse_stability.png")
)
make_rho_stability_plot(
  kegg_cont_compare,
  "KEGG pathway stability after common-gene-universe harmonization",
  file.path(out_dir, "Figures", "02_KEGG_commonUniverse_stability.png")
)

pretty_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

core_plot_df <- rbind(
  data.frame(
    Family = "Hallmark", Program = core_hall_result$GeneSet,
    EMTAB = core_hall_result$Rho_EMTAB_07C,
    SEQC = core_hall_result$Rho_SEQC_07C,
    RobustBoth = core_hall_result$RobustBoth_07C
  ),
  data.frame(
    Family = "KEGG", Program = core_kegg_result$GeneSet,
    EMTAB = core_kegg_result$Rho_EMTAB_07C,
    SEQC = core_kegg_result$Rho_SEQC_07C,
    RobustBoth = core_kegg_result$RobustBoth_07C
  )
)
core_long <- rbind(
  data.frame(Family = core_plot_df$Family, Program = core_plot_df$Program, Dataset = "E-MTAB-179", Rho = core_plot_df$EMTAB, RobustBoth = core_plot_df$RobustBoth),
  data.frame(Family = core_plot_df$Family, Program = core_plot_df$Program, Dataset = "SEQC", Rho = core_plot_df$SEQC, RobustBoth = core_plot_df$RobustBoth)
)
core_long$ProgramLabel <- pretty_name(core_long$Program)
core_long$Dataset <- factor(core_long$Dataset, levels = c("E-MTAB-179", "SEQC"))

p_core <- ggplot2::ggplot(
  core_long,
  ggplot2::aes(x = Dataset, y = ProgramLabel, fill = Rho)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = ifelse(RobustBoth, sprintf("%.2f*", Rho), sprintf("%.2f", Rho))),
    size = 3
  ) +
  ggplot2::facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  ggplot2::labs(
    title = "Core M070-associated programs under strict common-gene harmonization",
    subtitle = "* Direction-concordant and FDR < 0.05 in both cohorts",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) + theme_pub
ggplot2::ggsave(
  file.path(out_dir, "Figures", "03_CorePathways_commonUniverse_heatmap.png"),
  p_core, width = 9, height = 9.5, dpi = 320, bg = "white"
)

## ----------------------------- final summary ---------------------------------
hall_old_robust_group <- sum(hall_group_compare$RobustBoth_FIX3, na.rm = TRUE)
hall_new_robust_group <- sum(hall_group_compare$RobustBoth_07C, na.rm = TRUE)
hall_old_robust_cont <- sum(hall_cont_compare$RobustBoth_FIX3, na.rm = TRUE)
hall_new_robust_cont <- sum(hall_cont_compare$RobustBoth_07C, na.rm = TRUE)
kegg_old_robust_group <- sum(kegg_group_compare$RobustBoth_FIX3, na.rm = TRUE)
kegg_new_robust_group <- sum(kegg_group_compare$RobustBoth_07C, na.rm = TRUE)
kegg_old_robust_cont <- sum(kegg_cont_compare$RobustBoth_FIX3, na.rm = TRUE)
kegg_new_robust_cont <- sum(kegg_cont_compare$RobustBoth_07C, na.rm = TRUE)

core_hall_robust <- sum(core_hall_result$RobustBoth_07C, na.rm = TRUE)
core_kegg_robust <- sum(core_kegg_result$RobustBoth_07C, na.rm = TRUE)
core_hall_dir_stable <- sum(core_hall_result$DirectionStable_EMTAB & core_hall_result$DirectionStable_SEQC, na.rm = TRUE)
core_kegg_dir_stable <- sum(core_kegg_result$DirectionStable_EMTAB & core_kegg_result$DirectionStable_SEQC, na.rm = TRUE)

hypoxia_row <- core_hall_result[core_hall_result$GeneSet == "HALLMARK_HYPOXIA", , drop = FALSE]
hypoxia_robust <- nrow(hypoxia_row) == 1L && isTRUE(hypoxia_row$RobustBoth_07C)

final_summary <- data.frame(
  Item = c(
    "Frozen M070 MD5",
    "M070 refitted in 07C",
    "Features reselected in 07C",
    "Cutoff retuned in 07C",
    "SEQC survival outcomes used in 07C",
    "Common expression-gene universe N",
    "Hallmark sets evaluable with identical membership",
    "KEGG sets evaluable with identical membership",
    "Hallmark robust group results in FIX3",
    "Hallmark robust group results in 07C",
    "Hallmark robust continuous results in FIX3",
    "Hallmark robust continuous results in 07C",
    "KEGG robust group results in FIX3",
    "KEGG robust group results in 07C",
    "KEGG robust continuous results in FIX3",
    "KEGG robust continuous results in 07C",
    "Core Hallmark pathways direction-stable in both cohorts",
    "Core Hallmark pathways robust in 07C",
    "Core KEGG pathways direction-stable in both cohorts",
    "Core KEGG pathways robust in 07C",
    "Hallmark hypoxia robust in both cohorts after harmonization",
    "Interpretation"
  ),
  Value = c(
    model_md5, "FALSE", "FALSE", "FALSE", "FALSE",
    length(common_genes), length(hall_shared$sets), length(kegg_shared$sets),
    hall_old_robust_group, hall_new_robust_group,
    hall_old_robust_cont, hall_new_robust_cont,
    kegg_old_robust_group, kegg_new_robust_group,
    kegg_old_robust_cont, kegg_new_robust_cont,
    paste0(core_hall_dir_stable, "/", nrow(core_hall_result)),
    paste0(core_hall_robust, "/", nrow(core_hall_result)),
    paste0(core_kegg_dir_stable, "/", nrow(core_kegg_result)),
    paste0(core_kegg_robust, "/", nrow(core_kegg_result)),
    as.character(hypoxia_robust),
    "Strict cross-platform sensitivity analysis using identical expression background and identical pathway membership; association evidence only"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  final_summary,
  file.path(out_dir, "Tables", "18_FINAL_07C_SENSITIVITY_SUMMARY.csv"),
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
writeLines(
  c(
    "SC28 RECALC 07C completed successfully.",
    paste0("Frozen M070 MD5: ", model_md5),
    "M070 was not refitted; features were not reselected; cutoff was not retuned.",
    "SEQC survival outcomes were not used.",
    paste0("Identical common expression-gene universe: ", length(common_genes), " genes"),
    paste0("Hallmark evaluable sets: ", length(hall_shared$sets)),
    paste0("KEGG evaluable sets: ", length(kegg_shared$sets)),
    paste0("Core Hallmark robust in 07C: ", core_hall_robust, "/", nrow(core_hall_result)),
    paste0("Core KEGG robust in 07C: ", core_kegg_robust, "/", nrow(core_kegg_result)),
    paste0("Hallmark hypoxia robust in both cohorts: ", hypoxia_robust),
    paste0("Warnings captured: ", nrow(warning_log)),
    "Upload the ENTIRE SC28_RECALC_07C_CommonGeneUniverse_PathwaySensitivity_20260809 folder as ZIP for final audit."
  ),
  file.path(out_dir, "Logs", "03_run_summary.txt")
)

cat("\n============================================================\n")
cat("SC28 RECALC 07C completed successfully.\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("Frozen M070 MD5: ", model_md5, "\n", sep = "")
cat("Common expression-gene universe: ", length(common_genes), " genes\n", sep = "")
cat("Hallmark evaluable sets: ", length(hall_shared$sets), "\n", sep = "")
cat("KEGG evaluable sets: ", length(kegg_shared$sets), "\n", sep = "")
cat("Warnings captured: ", nrow(warning_log), "\n", sep = "")
cat("Upload the ENTIRE 07C output folder as ZIP for final audit.\n")
cat("============================================================\n")
