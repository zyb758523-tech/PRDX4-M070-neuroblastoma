rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260808)

## ============================================================================
## PRDX4_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism.R
##
## Purpose
##   Refresh the bulk-transcriptomic mechanism layer after the final prognostic
##   model was locked as M070 = SuperPC -> RSF.
##
## Frozen analysis boundary
##   1. M070 is NOT refitted here.
##   2. No genes/features are reselected here.
##   3. The E-MTAB-179-derived locked risk cutoff is NOT retuned here.
##   4. SEQC outcomes are NOT used for pathway or immune-signature selection.
##   5. E-MTAB-179 and SEQC are analysed independently, then direction/FDR
##      concordance is summarized across cohorts.
##   6. All pathway/immune results are association evidence, not causal proof.
##
## Analysis blocks
##   A. Hallmark ssGSEA + locked-group limma + continuous M070 correlation
##   B. KEGG ssGSEA + locked-group limma + continuous M070 correlation
##   C. Curated immune transcriptomic signatures (mean gene-wise z score)
##   D. Focus-gene bridge: PRDX4, MIF, CD74, CXCR4, CD44, LDHA vs M070
##   E. PRDX4 vs selected pathway/immune programs (descriptive; PRDX4 is a
##      component of M070 and is therefore NOT treated as independent evidence)
##
## FIX3 immune-signature safeguard
##   Immune signatures are evaluated only when at least 3 identical genes are
##   measurable and nonconstant in BOTH E-MTAB-179 and SEQC. The exact same
##   shared gene list is used in both cohorts. Signatures below this threshold
##   are recorded as not evaluable and are not forced into the analysis.
##
## Software note
##   Current GSVA releases use ssgseaParam() -> gsva(). A compatibility branch
##   for older GSVA releases is retained so the same script can run if the
##   existing R library is older.
## ============================================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
out_dir <- file.path(
  root_dir,
  "SC28_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism_20260809"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

## ----------------------------- fixed input files -----------------------------
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
m070_lock_record_file <- file.path(
  root_dir,
  "SC28_RECALC_06_LockM070_SEQC_External_20260808", "Tables",
  "05_LOCK_RECORD_BEFORE_SEQC.csv"
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

required_files <- c(
  emtab_expr_file,
  seqc_ready_file,
  m070_lock_record_file,
  emtab_score_file,
  seqc_score_file,
  locked_model_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

## ----------------------------- required packages -----------------------------
required_packages <- c("GSVA", "limma", "msigdbr", "ggplot2")
pkg_ok <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
write.csv(
  data.frame(Package = required_packages, Available = unname(pkg_ok)),
  file.path(out_dir, "Tables", "00_package_status.csv"),
  row.names = FALSE
)
if (any(!pkg_ok)) {
  missing_pkgs <- required_packages[!pkg_ok]
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "), "\n\n",
    "Install CRAN packages with:\n",
    "install.packages(c(\"msigdbr\", \"ggplot2\"), repos=\"https://cloud.r-project.org\")\n",
    "Install Bioconductor packages with:\n",
    "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"GSVA\", \"limma\"))\n",
    "Then restart R and rerun this WHOLE script from line 1."
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
  if (!all(c("Sample", score_col, group_col) %in% names(scores))) {
    stop(dataset_label, " score table lacks required columns")
  }
  scores$Sample <- as.character(scores$Sample)
  scores$M070 <- safe_num(scores[[score_col]])
  scores$RiskGroup <- factor(as.character(scores[[group_col]]), levels = c("Low", "High"))
  if (any(!is.finite(scores$M070)) || any(is.na(scores$RiskGroup))) {
    stop(dataset_label, " contains invalid M070 score or locked risk group")
  }

  expr_names <- colnames(expr)
  exact <- match(scores$Sample, expr_names)
  if (anyNA(exact)) {
    ek <- sample_key(expr_names)
    sk <- sample_key(scores$Sample)
    if (anyDuplicated(ek)) stop(dataset_label, " expression sample keys are duplicated")
    if (anyDuplicated(sk)) stop(dataset_label, " score sample keys are duplicated")
    exact <- match(sk, ek)
  }
  if (anyNA(exact)) {
    bad <- scores$Sample[is.na(exact)]
    stop(
      dataset_label, " score-to-expression matching failed for ", length(bad),
      " sample(s), e.g. ", paste(head(bad, 8), collapse = ", ")
    )
  }
  expr2 <- expr[, exact, drop = FALSE]
  colnames(expr2) <- scores$Sample
  if (ncol(expr2) != expected_n || nrow(scores) != expected_n) {
    stop(
      dataset_label, " unexpected sample count after alignment. Expected ",
      expected_n, "; expression=", ncol(expr2), "; scores=", nrow(scores)
    )
  }
  list(expr = expr2, scores = scores)
}

extract_lock_value <- function(lock_tab, item_pattern) {
  if (!all(c("Item", "Value") %in% names(lock_tab))) return(NA_character_)
  ix <- grep(item_pattern, lock_tab$Item, ignore.case = TRUE)
  if (!length(ix)) return(NA_character_)
  as.character(lock_tab$Value[ix[1]])
}

pretty_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

## ------------------------------ freeze checks --------------------------------
## FIX2: use the authoritative frozen RDS fields directly. The human-readable
## lock-record CSV remains an audit artifact only and is not used to recover the
## numerical cutoff.
lock_tab <- read.csv(m070_lock_record_file, check.names = FALSE)
frozen_m070 <- readRDS(locked_model_file)
required_frozen_fields <- c(
  "ModelID", "Selector", "Learner", "SelectedGenes", "LockedRiskCutoff",
  "SEQCUsedForModelSelection", "SEQCUsedForCutoff"
)
missing_frozen_fields <- setdiff(required_frozen_fields, names(frozen_m070))
if (length(missing_frozen_fields)) {
  stop(
    "Frozen M070 RDS lacks required field(s): ",
    paste(missing_frozen_fields, collapse = ", ")
  )
}

lock_model_id <- as.character(frozen_m070$ModelID)[1]
lock_selector <- as.character(frozen_m070$Selector)[1]
lock_learner <- as.character(frozen_m070$Learner)[1]
lock_selected_vec <- as.character(frozen_m070$SelectedGenes)
lock_selected <- paste(lock_selected_vec, collapse = ";")
lock_prdx4 <- ifelse("PRDX4" %in% lock_selected_vec, "TRUE", "FALSE")
lock_cutoff <- safe_num(frozen_m070$LockedRiskCutoff)[1]

if (!identical(lock_model_id, "M070_SuperPC_RSF")) stop("Unexpected locked model ID: ", lock_model_id)
if (!identical(lock_selector, "SuperPC") || !identical(lock_learner, "RSF")) {
  stop("Frozen M070 RDS is not SuperPC -> RSF")
}
if (!identical(lock_prdx4, "TRUE")) stop("PRDX4 is not recorded as selected in frozen M070")
if (!is.finite(lock_cutoff)) {
  stop("LockedRiskCutoff inside frozen M070 RDS is not finite")
}
if (!identical(frozen_m070$SEQCUsedForModelSelection, FALSE)) {
  stop("Frozen M070 provenance indicates SEQC was used for model selection")
}
if (!identical(frozen_m070$SEQCUsedForCutoff, FALSE)) {
  stop("Frozen M070 provenance indicates SEQC was used for cutoff selection")
}

model_md5 <- unname(tools::md5sum(locked_model_file))
expected_md5 <- "e719d0e6b626e75cc134e83fd0929083"
if (!identical(model_md5, expected_md5)) {
  stop("Frozen M070 MD5 mismatch. Expected ", expected_md5, "; observed ", model_md5)
}

input_audit <- data.frame(
  Item = c(
    "Locked model ID", "Selector", "Learner", "Frozen M070 MD5",
    "Expected frozen MD5", "PRDX4 selected", "Locked risk cutoff",
    "M070 refitted in 07B", "Features reselected in 07B",
    "Risk cutoff retuned in 07B", "SEQC outcomes used for mechanism selection"
  ),
  Value = c(
    lock_model_id, lock_selector, lock_learner, model_md5, expected_md5,
    lock_prdx4, lock_cutoff, "FALSE", "FALSE", "FALSE", "FALSE"
  )
)
write.csv(
  input_audit,
  file.path(out_dir, "Tables", "01_frozen_M070_and_design_audit.csv"),
  row.names = FALSE
)

## --------------------------- load/align expression ----------------------------
emtab_expr <- collapse_duplicate_genes(as_expr_matrix(readRDS(emtab_expr_file)))
seqc_obj <- readRDS(seqc_ready_file)
if (is.null(seqc_obj$expr)) stop("SEQC ready object lacks $expr")
seqc_expr <- collapse_duplicate_genes(as_expr_matrix(seqc_obj$expr))

emtab_scores <- read.csv(emtab_score_file, check.names = FALSE)
seqc_scores <- read.csv(seqc_score_file, check.names = FALSE)

emtab <- align_expression_to_scores(
  emtab_expr, emtab_scores,
  score_col = "M070_RiskScore",
  group_col = "RiskGroup_LockedCutoff",
  dataset_label = "E-MTAB-179", expected_n = 478L
)
seqc <- align_expression_to_scores(
  seqc_expr, seqc_scores,
  score_col = "Corrected06B_M070_RiskScore",
  group_col = "RiskGroup_Locked_EMTAB_Cutoff",
  dataset_label = "SEQC", expected_n = 498L
)

## Verify that the saved groups are exactly those produced by the frozen
## E-MTAB-179 cutoff. This catches accidental median re-splitting in 07B.
emtab_group_from_cutoff <- ifelse(emtab$scores$M070 >= lock_cutoff, "High", "Low")
seqc_group_from_cutoff <- ifelse(seqc$scores$M070 >= lock_cutoff, "High", "Low")
if (!all(emtab_group_from_cutoff == as.character(emtab$scores$RiskGroup))) {
  stop("E-MTAB-179 saved risk groups do not match the frozen M070 cutoff")
}
if (!all(seqc_group_from_cutoff == as.character(seqc$scores$RiskGroup))) {
  stop("SEQC saved risk groups do not match the frozen E-MTAB-179 M070 cutoff")
}

alignment_audit <- rbind(
  data.frame(
    Dataset = "E-MTAB-179", N = ncol(emtab$expr),
    LowN = sum(emtab$scores$RiskGroup == "Low"),
    HighN = sum(emtab$scores$RiskGroup == "High"),
    GeneN = nrow(emtab$expr), PRDX4 = "PRDX4" %in% rownames(emtab$expr),
    MIF = "MIF" %in% rownames(emtab$expr)
  ),
  data.frame(
    Dataset = "SEQC", N = ncol(seqc$expr),
    LowN = sum(seqc$scores$RiskGroup == "Low"),
    HighN = sum(seqc$scores$RiskGroup == "High"),
    GeneN = nrow(seqc$expr), PRDX4 = "PRDX4" %in% rownames(seqc$expr),
    MIF = "MIF" %in% rownames(seqc$expr)
  )
)
write.csv(
  alignment_audit,
  file.path(out_dir, "Tables", "02_expression_M070_alignment_audit.csv"),
  row.names = FALSE
)

## ------------------------------- MSigDB sets ----------------------------------
get_msigdb_sets <- function(collection, subcollection = NULL) {
  fml <- names(formals(msigdbr::msigdbr))
  args <- list(species = "Homo sapiens")
  if ("collection" %in% fml) {
    args$collection <- collection
    if (!is.null(subcollection)) args$subcollection <- subcollection
  } else {
    args$category <- collection
    if (!is.null(subcollection)) args$subcategory <- subcollection
  }
  tab <- do.call(msigdbr::msigdbr, args)
  gene_col <- intersect(c("gene_symbol", "human_gene_symbol"), names(tab))[1]
  set_col <- intersect(c("gs_name", "gs_name"), names(tab))[1]
  if (is.na(gene_col) || is.na(set_col)) stop("Cannot identify MSigDB gene-set columns")
  sets <- split(toupper(as.character(tab[[gene_col]])), as.character(tab[[set_col]]))
  sets <- lapply(sets, function(g) unique(g[nzchar(g) & !is.na(g)]))
  list(table = tab, sets = sets)
}

hallmark_obj <- with_warnlog("msigdbr Hallmark retrieval", get_msigdb_sets("H"))
hallmark_sets <- hallmark_obj$sets

## Prefer KEGG Legacy to preserve continuity with the older manuscript. If the
## installed MSigDB snapshot lacks it, fall back to KEGG Medicus and record it.
kegg_subcollection <- "CP:KEGG_LEGACY"
kegg_obj <- tryCatch(
  with_warnlog(
    "msigdbr KEGG Legacy retrieval",
    get_msigdb_sets("C2", kegg_subcollection)
  ),
  error = function(e) NULL
)
if (is.null(kegg_obj) || !length(kegg_obj$sets)) {
  kegg_subcollection <- "CP:KEGG_MEDICUS"
  kegg_obj <- with_warnlog(
    "msigdbr KEGG Medicus retrieval",
    get_msigdb_sets("C2", kegg_subcollection)
  )
}
kegg_sets <- kegg_obj$sets

msigdb_audit <- data.frame(
  Collection = c("Hallmark", "KEGG"),
  MSigDB_collection = c("H", "C2"),
  MSigDB_subcollection = c(NA, kegg_subcollection),
  GeneSetN = c(length(hallmark_sets), length(kegg_sets)),
  stringsAsFactors = FALSE
)
write.csv(
  msigdb_audit,
  file.path(out_dir, "Tables", "03_MSigDB_collection_audit.csv"),
  row.names = FALSE
)

## Save exact gene-set membership used for reproducibility.
gene_set_long <- function(sets, collection_label) {
  do.call(
    rbind,
    lapply(names(sets), function(nm) {
      data.frame(Collection = collection_label, GeneSet = nm, Gene = sets[[nm]])
    })
  )
}
write.csv(
  rbind(
    gene_set_long(hallmark_sets, "Hallmark"),
    gene_set_long(kegg_sets, paste0("KEGG_", kegg_subcollection))
  ),
  file.path(out_dir, "Tables", "04_MSigDB_gene_set_membership_USED.csv"),
  row.names = FALSE
)

## ----------------------------- ssGSEA functions -------------------------------
filter_sets_for_expr <- function(sets, expr, min_size = 5L, max_size = 500L) {
  ans <- lapply(sets, function(g) intersect(unique(toupper(g)), rownames(expr)))
  keep <- lengths(ans) >= min_size & lengths(ans) <= max_size
  ans[keep]
}

run_ssgsea <- function(expr, sets, context) {
  sets2 <- filter_sets_for_expr(sets, expr, min_size = 5L, max_size = 500L)
  if (!length(sets2)) stop(context, ": no usable gene sets after expression matching")
  ex <- getNamespaceExports("GSVA")
  if ("ssgseaParam" %in% ex) {
    par <- with_warnlog(
      paste(context, "ssgseaParam"),
      GSVA::ssgseaParam(
        expr, sets2,
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
        expr, sets2, method = "ssgsea",
        min.sz = 5, max.sz = 500,
        ssgsea.norm = TRUE,
        parallel.sz = 1, verbose = FALSE
      )
    )
  }
  es <- as.matrix(es)
  if (any(!is.finite(es))) stop(context, ": non-finite ssGSEA scores produced")
  es
}

limma_group_test <- function(score_mat, score_df, dataset, analysis_family) {
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
    Dataset = dataset,
    AnalysisFamily = analysis_family,
    GeneSet = tt$GeneSet,
    MeanDifference_HighMinusLow = tt$logFC,
    T = tt$t,
    P = tt$P.Value,
    FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}

continuous_spearman <- function(score_mat, score_df, dataset, analysis_family) {
  out <- lapply(seq_len(nrow(score_mat)), function(i) {
    ct <- suppressWarnings(cor.test(
      safe_num(score_mat[i, ]), score_df$M070,
      method = "spearman", exact = FALSE
    ))
    data.frame(
      Dataset = dataset,
      AnalysisFamily = analysis_family,
      GeneSet = rownames(score_mat)[i],
      Rho = unname(ct$estimate),
      P = ct$p.value
    )
  })
  out <- do.call(rbind, out)
  out$FDR <- p.adjust(out$P, method = "BH")
  out
}

cross_cohort_group <- function(a, b, family_label) {
  x <- merge(
    a[, c("GeneSet", "MeanDifference_HighMinusLow", "P", "FDR")],
    b[, c("GeneSet", "MeanDifference_HighMinusLow", "P", "FDR")],
    by = "GeneSet", suffixes = c("_EMTAB", "_SEQC")
  )
  x$DirectionConcordant <- sign(x$MeanDifference_HighMinusLow_EMTAB) == sign(x$MeanDifference_HighMinusLow_SEQC)
  x$FDRlt005_Both <- x$FDR_EMTAB < 0.05 & x$FDR_SEQC < 0.05
  x$ConcordanceClass <- ifelse(
    x$DirectionConcordant & x$FDRlt005_Both & x$MeanDifference_HighMinusLow_EMTAB > 0,
    "Concordant high-risk enriched; FDR < 0.05 in both",
    ifelse(
      x$DirectionConcordant & x$FDRlt005_Both & x$MeanDifference_HighMinusLow_EMTAB < 0,
      "Concordant low-risk enriched; FDR < 0.05 in both",
      ifelse(x$DirectionConcordant, "Direction-concordant only", "Direction-discordant")
    )
  )
  x$AnalysisFamily <- family_label
  x
}

cross_cohort_cor <- function(a, b, family_label) {
  x <- merge(
    a[, c("GeneSet", "Rho", "P", "FDR")],
    b[, c("GeneSet", "Rho", "P", "FDR")],
    by = "GeneSet", suffixes = c("_EMTAB", "_SEQC")
  )
  x$DirectionConcordant <- sign(x$Rho_EMTAB) == sign(x$Rho_SEQC)
  x$FDRlt005_Both <- x$FDR_EMTAB < 0.05 & x$FDR_SEQC < 0.05
  x$AnalysisFamily <- family_label
  x
}

write_score_matrix <- function(mat, path) {
  write.csv(
    data.frame(GeneSet = rownames(mat), mat, check.names = FALSE),
    path, row.names = FALSE
  )
}

## ------------------------------ Hallmark ssGSEA -------------------------------
hall_emtab <- run_ssgsea(emtab$expr, hallmark_sets, "E-MTAB-179 Hallmark")
hall_seqc <- run_ssgsea(seqc$expr, hallmark_sets, "SEQC Hallmark")
write_score_matrix(hall_emtab, file.path(out_dir, "Tables", "05_Hallmark_ssGSEA_EMTAB_scores.csv"))
write_score_matrix(hall_seqc, file.path(out_dir, "Tables", "06_Hallmark_ssGSEA_SEQC_scores.csv"))

hall_grp_emtab <- limma_group_test(hall_emtab, emtab$scores, "E-MTAB-179", "Hallmark")
hall_grp_seqc <- limma_group_test(hall_seqc, seqc$scores, "SEQC", "Hallmark")
hall_cor_emtab <- continuous_spearman(hall_emtab, emtab$scores, "E-MTAB-179", "Hallmark")
hall_cor_seqc <- continuous_spearman(hall_seqc, seqc$scores, "SEQC", "Hallmark")
hall_cross_grp <- cross_cohort_group(hall_grp_emtab, hall_grp_seqc, "Hallmark")
hall_cross_cor <- cross_cohort_cor(hall_cor_emtab, hall_cor_seqc, "Hallmark")

write.csv(hall_grp_emtab, file.path(out_dir, "Tables", "07_Hallmark_EMTAB_locked_group_limma.csv"), row.names = FALSE)
write.csv(hall_grp_seqc, file.path(out_dir, "Tables", "08_Hallmark_SEQC_locked_group_limma.csv"), row.names = FALSE)
write.csv(rbind(hall_cor_emtab, hall_cor_seqc), file.path(out_dir, "Tables", "09_Hallmark_continuous_M070_Spearman.csv"), row.names = FALSE)
write.csv(hall_cross_grp, file.path(out_dir, "Tables", "10_Hallmark_cross_cohort_locked_group_concordance.csv"), row.names = FALSE)
write.csv(hall_cross_cor, file.path(out_dir, "Tables", "11_Hallmark_cross_cohort_continuous_concordance.csv"), row.names = FALSE)

## -------------------------------- KEGG ssGSEA --------------------------------
kegg_emtab <- run_ssgsea(emtab$expr, kegg_sets, "E-MTAB-179 KEGG")
kegg_seqc <- run_ssgsea(seqc$expr, kegg_sets, "SEQC KEGG")
write_score_matrix(kegg_emtab, file.path(out_dir, "Tables", "12_KEGG_ssGSEA_EMTAB_scores.csv"))
write_score_matrix(kegg_seqc, file.path(out_dir, "Tables", "13_KEGG_ssGSEA_SEQC_scores.csv"))

kegg_grp_emtab <- limma_group_test(kegg_emtab, emtab$scores, "E-MTAB-179", "KEGG")
kegg_grp_seqc <- limma_group_test(kegg_seqc, seqc$scores, "SEQC", "KEGG")
kegg_cor_emtab <- continuous_spearman(kegg_emtab, emtab$scores, "E-MTAB-179", "KEGG")
kegg_cor_seqc <- continuous_spearman(kegg_seqc, seqc$scores, "SEQC", "KEGG")
kegg_cross_grp <- cross_cohort_group(kegg_grp_emtab, kegg_grp_seqc, "KEGG")
kegg_cross_cor <- cross_cohort_cor(kegg_cor_emtab, kegg_cor_seqc, "KEGG")

write.csv(kegg_grp_emtab, file.path(out_dir, "Tables", "14_KEGG_EMTAB_locked_group_limma.csv"), row.names = FALSE)
write.csv(kegg_grp_seqc, file.path(out_dir, "Tables", "15_KEGG_SEQC_locked_group_limma.csv"), row.names = FALSE)
write.csv(rbind(kegg_cor_emtab, kegg_cor_seqc), file.path(out_dir, "Tables", "16_KEGG_continuous_M070_Spearman.csv"), row.names = FALSE)
write.csv(kegg_cross_grp, file.path(out_dir, "Tables", "17_KEGG_cross_cohort_locked_group_concordance.csv"), row.names = FALSE)
write.csv(kegg_cross_cor, file.path(out_dir, "Tables", "18_KEGG_cross_cohort_continuous_concordance.csv"), row.names = FALSE)

## ------------------------- curated immune signatures --------------------------
## These are transcriptomic marker signatures, not direct cell-count estimates.
immune_sets <- list(
  "Cytotoxic T/NK" = c("CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "CTSW"),
  "CD8 T cells" = c("CD8A", "CD8B", "CD3D", "CD3E", "TRAC", "CCL5", "LCK"),
  "CD4 T cells" = c("CD4", "CD3D", "CD3E", "IL7R", "CCR7", "LTB", "MAL"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CCR8"),
  "NK cells" = c("NKG7", "GNLY", "KLRD1", "KLRK1", "PRF1", "FCGR3A", "TRDC"),
  "B cells" = c("CD19", "MS4A1", "CD79A", "CD37", "CD74", "CD22", "CD83"),
  "Plasma cells" = c("MZB1", "JCHAIN", "SDC1", "XBP1", "IGKC", "DERL3"),
  "Macrophage/monocyte" = c("LST1", "LILRB1", "CSF1R", "CD14", "FCGR3A", "CTSS", "TYROBP"),
  "M2 macrophage" = c("MRC1", "CD163", "MSR1", "C1QA", "C1QB", "C1QC", "TGFB1"),
  "Dendritic cells" = c("FCER1A", "CD1C", "CLEC10A", "CST3", "HLA-DRA", "ITGAX"),
  "Neutrophil" = c("FCGR3B", "CSF3R", "S100A8", "S100A9", "CXCR2", "FPR1"),
  "Antigen presentation" = c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "TAP1", "TAP2", "PSMB8", "PSMB9"),
  "IFN-gamma response" = c("STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "IDO1", "IFIT1", "IFIT3"),
  "Immune checkpoint" = c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "PDCD1LG2"),
  "Stromal/fibroblast" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL5A1", "FAP", "PDGFRA"),
  "Endothelial" = c("PECAM1", "VWF", "EMCN", "KDR", "CDH5", "RAMP2")
)
immune_sets <- lapply(immune_sets, unique)

## FIX3: build one common, nonconstant gene definition for every immune
## signature before scoring. This avoids two problems at once:
##   (i) aborting the entire mechanism analysis when a sparse microarray
##       signature has fewer than 3 measurable genes; and
##   (ii) comparing nominally identical signatures that were actually computed
##        from different genes in E-MTAB-179 and SEQC.
## SEQC survival outcomes are not used here; this is expression-coverage QC.
immune_signature_coverage <- lapply(names(immune_sets), function(nm) {
  planned <- unique(toupper(immune_sets[[nm]]))
  em_avail <- intersect(planned, rownames(emtab$expr))
  sq_avail <- intersect(planned, rownames(seqc$expr))
  common <- intersect(em_avail, sq_avail)

  if (length(common)) {
    em_sd <- apply(emtab$expr[common, , drop = FALSE], 1, sd)
    sq_sd <- apply(seqc$expr[common, , drop = FALSE], 1, sd)
    keep <- is.finite(em_sd) & em_sd > 0 & is.finite(sq_sd) & sq_sd > 0
    common_variable <- common[keep]
  } else {
    common_variable <- character()
  }

  data.frame(
    Signature = nm,
    PlannedGeneN = length(planned),
    EMTABAvailableGeneN = length(em_avail),
    SEQCAvailableGeneN = length(sq_avail),
    SharedAvailableGeneN = length(common),
    SharedNonconstantGeneN = length(common_variable),
    EMTABAvailableGenes = paste(em_avail, collapse = ";"),
    SEQCAvailableGenes = paste(sq_avail, collapse = ";"),
    SharedNonconstantGenes = paste(common_variable, collapse = ";"),
    EvaluableInBoth = length(common_variable) >= 3L,
    Status = ifelse(
      length(common_variable) >= 3L,
      "evaluable_same_genes_both_cohorts",
      "not_evaluable_lt3_shared_nonconstant_genes"
    ),
    stringsAsFactors = FALSE
  )
})
immune_signature_coverage <- do.call(rbind, immune_signature_coverage)

shared_immune_sets <- setNames(
  lapply(seq_len(nrow(immune_signature_coverage)), function(i) {
    x <- immune_signature_coverage$SharedNonconstantGenes[i]
    if (!nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1]]
  }),
  immune_signature_coverage$Signature
)
shared_immune_sets <- shared_immune_sets[immune_signature_coverage$EvaluableInBoth]

if (!length(shared_immune_sets)) {
  stop("No immune signature has >=3 shared nonconstant genes in both cohorts")
}

write.csv(
  immune_signature_coverage,
  file.path(out_dir, "Tables", "19_immune_signature_gene_coverage.csv"),
  row.names = FALSE
)

not_evaluable <- immune_signature_coverage$Signature[!immune_signature_coverage$EvaluableInBoth]
if (length(not_evaluable)) {
  message(
    "Immune signatures recorded as not evaluable (<3 shared nonconstant genes): ",
    paste(not_evaluable, collapse = ", ")
  )
}

immune_signature_scores <- function(expr, sets, dataset) {
  score_mat <- do.call(
    rbind,
    lapply(names(sets), function(nm) {
      genes <- sets[[nm]]
      if (length(genes) < 3L || !all(genes %in% rownames(expr))) {
        stop(dataset, " shared-gene QC failed unexpectedly for ", nm)
      }
      x <- expr[genes, , drop = FALSE]
      sdv <- apply(x, 1, sd)
      if (any(!is.finite(sdv) | sdv <= 0)) {
        stop(dataset, " shared-gene variance QC failed unexpectedly for ", nm)
      }
      z <- t(scale(t(x)))
      colMeans(z)
    })
  )
  rownames(score_mat) <- names(sets)
  colnames(score_mat) <- colnames(expr)
  list(score = score_mat)
}

imm_emtab <- immune_signature_scores(emtab$expr, shared_immune_sets, "E-MTAB-179")
imm_seqc <- immune_signature_scores(seqc$expr, shared_immune_sets, "SEQC")
write_score_matrix(imm_emtab$score, file.path(out_dir, "Tables", "20_immune_signature_EMTAB_scores.csv"))
write_score_matrix(imm_seqc$score, file.path(out_dir, "Tables", "21_immune_signature_SEQC_scores.csv"))

immune_assoc <- function(score_mat, score_df, dataset) {
  out <- lapply(seq_len(nrow(score_mat)), function(i) {
    y <- safe_num(score_mat[i, ])
    ct <- suppressWarnings(cor.test(y, score_df$M070, method = "spearman", exact = FALSE))
    lo <- y[score_df$RiskGroup == "Low"]
    hi <- y[score_df$RiskGroup == "High"]
    wt <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
    data.frame(
      Dataset = dataset,
      Signature = rownames(score_mat)[i],
      Rho_M070 = unname(ct$estimate),
      SpearmanP = ct$p.value,
      MeanLow = mean(lo),
      MeanHigh = mean(hi),
      MeanDifference_HighMinusLow = mean(hi) - mean(lo),
      WilcoxonP = wt$p.value
    )
  })
  out <- do.call(rbind, out)
  out$SpearmanFDR <- p.adjust(out$SpearmanP, method = "BH")
  out$WilcoxonFDR <- p.adjust(out$WilcoxonP, method = "BH")
  out
}

imm_assoc_emtab <- immune_assoc(imm_emtab$score, emtab$scores, "E-MTAB-179")
imm_assoc_seqc <- immune_assoc(imm_seqc$score, seqc$scores, "SEQC")
immune_cross <- merge(
  imm_assoc_emtab,
  imm_assoc_seqc,
  by = "Signature", suffixes = c("_EMTAB", "_SEQC")
)
immune_cross$RhoDirectionConcordant <- sign(immune_cross$Rho_M070_EMTAB) == sign(immune_cross$Rho_M070_SEQC)
immune_cross$SpearmanFDRlt005_Both <- immune_cross$SpearmanFDR_EMTAB < 0.05 & immune_cross$SpearmanFDR_SEQC < 0.05
immune_cross$GroupDirectionConcordant <- sign(immune_cross$MeanDifference_HighMinusLow_EMTAB) == sign(immune_cross$MeanDifference_HighMinusLow_SEQC)
immune_cross$WilcoxonFDRlt005_Both <- immune_cross$WilcoxonFDR_EMTAB < 0.05 & immune_cross$WilcoxonFDR_SEQC < 0.05

write.csv(
  rbind(imm_assoc_emtab, imm_assoc_seqc),
  file.path(out_dir, "Tables", "22_immune_signature_M070_associations.csv"),
  row.names = FALSE
)
write.csv(
  immune_cross,
  file.path(out_dir, "Tables", "23_immune_signature_cross_cohort_concordance.csv"),
  row.names = FALSE
)

## ----------------------------- focus-gene bridge ------------------------------
focus_genes <- c("PRDX4", "MIF", "CD74", "CXCR4", "CD44", "LDHA")
focus_gene_assoc <- function(expr, score_df, dataset) {
  out <- lapply(focus_genes, function(g) {
    if (!g %in% rownames(expr)) {
      return(data.frame(
        Dataset = dataset, Gene = g, Available = FALSE,
        Rho_M070 = NA, SpearmanP = NA,
        MeanLow = NA, MeanHigh = NA,
        MeanDifference_HighMinusLow = NA, WilcoxonP = NA
      ))
    }
    y <- safe_num(expr[g, ])
    ct <- suppressWarnings(cor.test(y, score_df$M070, method = "spearman", exact = FALSE))
    lo <- y[score_df$RiskGroup == "Low"]
    hi <- y[score_df$RiskGroup == "High"]
    wt <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
    data.frame(
      Dataset = dataset, Gene = g, Available = TRUE,
      Rho_M070 = unname(ct$estimate), SpearmanP = ct$p.value,
      MeanLow = mean(lo), MeanHigh = mean(hi),
      MeanDifference_HighMinusLow = mean(hi) - mean(lo),
      WilcoxonP = wt$p.value
    )
  })
  out <- do.call(rbind, out)
  ok <- which(out$Available)
  out$SpearmanFDR <- NA_real_
  out$WilcoxonFDR <- NA_real_
  out$SpearmanFDR[ok] <- p.adjust(out$SpearmanP[ok], method = "BH")
  out$WilcoxonFDR[ok] <- p.adjust(out$WilcoxonP[ok], method = "BH")
  out
}

focus_emtab <- focus_gene_assoc(emtab$expr, emtab$scores, "E-MTAB-179")
focus_seqc <- focus_gene_assoc(seqc$expr, seqc$scores, "SEQC")
focus_all <- rbind(focus_emtab, focus_seqc)
write.csv(
  focus_all,
  file.path(out_dir, "Tables", "24_focus_genes_PRDX4_MIF_axis_vs_M070.csv"),
  row.names = FALSE
)

## -------------------- PRDX4 vs mechanism programs (descriptive) ---------------
focus_hallmark_patterns <- c(
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
focus_immune_names <- c(
  "Cytotoxic T/NK", "CD8 T cells", "NK cells", "Antigen presentation",
  "IFN-gamma response", "Immune checkpoint", "M2 macrophage",
  "Stromal/fibroblast"
)

prdx4_program_cor <- function(expr, hall_mat, immune_mat, dataset) {
  if (!"PRDX4" %in% rownames(expr)) stop(dataset, " lacks PRDX4")
  px <- safe_num(expr["PRDX4", ])
  hm <- intersect(focus_hallmark_patterns, rownames(hall_mat))
  im <- intersect(focus_immune_names, rownames(immune_mat))
  program_list <- c(
    setNames(lapply(hm, function(nm) hall_mat[nm, ]), hm),
    setNames(lapply(im, function(nm) immune_mat[nm, ]), im)
  )
  out <- lapply(names(program_list), function(nm) {
    y <- safe_num(program_list[[nm]])
    ct <- suppressWarnings(cor.test(px, y, method = "spearman", exact = FALSE))
    data.frame(
      Dataset = dataset,
      Program = nm,
      ProgramType = ifelse(nm %in% hm, "Hallmark", "Immune signature"),
      Rho_PRDX4 = unname(ct$estimate),
      P = ct$p.value
    )
  })
  out <- do.call(rbind, out)
  out$FDR <- p.adjust(out$P, method = "BH")
  out
}

prdx4_prog_emtab <- prdx4_program_cor(emtab$expr, hall_emtab, imm_emtab$score, "E-MTAB-179")
prdx4_prog_seqc <- prdx4_program_cor(seqc$expr, hall_seqc, imm_seqc$score, "SEQC")
write.csv(
  rbind(prdx4_prog_emtab, prdx4_prog_seqc),
  file.path(out_dir, "Tables", "25_PRDX4_selected_program_correlations.csv"),
  row.names = FALSE
)

## --------------------------------- figures ------------------------------------
theme_pub <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    legend.title = ggplot2::element_text(face = "bold")
  )

plot_cross_cohort <- function(tab, title, output_file, label_n = 14L) {
  tab$Label <- pretty_name(tab$GeneSet)
  ord <- order(
    tab$FDRlt005_Both,
    abs(tab$MeanDifference_HighMinusLow_EMTAB) + abs(tab$MeanDifference_HighMinusLow_SEQC),
    decreasing = TRUE
  )
  lab_idx <- head(ord, min(label_n, nrow(tab)))
  tab$PointLabel <- ""
  tab$PointLabel[lab_idx] <- tab$Label[lab_idx]
  p <- ggplot2::ggplot(
    tab,
    ggplot2::aes(
      x = MeanDifference_HighMinusLow_EMTAB,
      y = MeanDifference_HighMinusLow_SEQC,
      color = ConcordanceClass
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
    ggplot2::geom_point(size = 2.2, alpha = 0.85) +
    ggplot2::geom_text(
      data = tab[tab$PointLabel != "", , drop = FALSE],
      ggplot2::aes(label = PointLabel),
      color = "black", size = 3, check_overlap = TRUE, vjust = -0.55
    ) +
    ggplot2::scale_color_manual(values = c(
      "Concordant high-risk enriched; FDR < 0.05 in both" = "#B2182B",
      "Concordant low-risk enriched; FDR < 0.05 in both" = "#2166AC",
      "Direction-concordant only" = "#999999",
      "Direction-discordant" = "#D9D9D9"
    )) +
    ggplot2::labs(
      title = title,
      subtitle = "Locked M070 risk groups; high minus low",
      x = "E-MTAB-179 mean ssGSEA difference",
      y = "SEQC mean ssGSEA difference",
      color = "Cross-cohort status"
    ) + theme_pub
  ggplot2::ggsave(output_file, p, width = 9.5, height = 7, dpi = 320, bg = "white")
}

plot_cross_cohort(
  hall_cross_grp,
  "Hallmark programs associated with the locked M070 high-risk state",
  file.path(out_dir, "Figures", "01_Hallmark_M070_cross_cohort.png"),
  label_n = 16L
)
plot_cross_cohort(
  kegg_cross_grp,
  paste0("KEGG programs associated with the locked M070 high-risk state (", kegg_subcollection, ")"),
  file.path(out_dir, "Figures", "02_KEGG_M070_cross_cohort.png"),
  label_n = 18L
)

## Hallmark focus heatmap: continuous M070 correlation in both cohorts.
hall_focus <- rbind(hall_cor_emtab, hall_cor_seqc)
hall_focus <- hall_focus[hall_focus$GeneSet %in% focus_hallmark_patterns, , drop = FALSE]
hall_focus$Program <- pretty_name(hall_focus$GeneSet)
hall_focus$Dataset <- factor(hall_focus$Dataset, levels = c("E-MTAB-179", "SEQC"))
p_hf <- ggplot2::ggplot(
  hall_focus,
  ggplot2::aes(x = Dataset, y = Program, fill = Rho)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Rho)), size = 3.4) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  ggplot2::labs(
    title = "Selected biological programs associated with continuous M070 risk",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) + theme_pub
ggplot2::ggsave(
  file.path(out_dir, "Figures", "03_Hallmark_focus_M070_heatmap.png"),
  p_hf, width = 7.5, height = 5.7, dpi = 320, bg = "white"
)

## Immune lollipop plot, continuous score correlation.
imm_long <- rbind(imm_assoc_emtab, imm_assoc_seqc)
imm_long$Dataset <- factor(imm_long$Dataset, levels = c("E-MTAB-179", "SEQC"))
imm_long$Direction <- ifelse(imm_long$Rho_M070 >= 0, "Positive", "Negative")
p_imm <- ggplot2::ggplot(
  imm_long,
  ggplot2::aes(x = Rho_M070, y = Signature, color = Direction)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Rho_M070, yend = Signature),
    linewidth = 0.6, color = "grey70"
  ) +
  ggplot2::geom_point(ggplot2::aes(size = -log10(pmax(SpearmanFDR, 1e-300))), alpha = 0.9) +
  ggplot2::facet_wrap(~Dataset, nrow = 1, scales = "free_x") +
  ggplot2::scale_color_manual(values = c("Negative" = "#2166AC", "Positive" = "#B2182B")) +
  ggplot2::labs(
    title = "Immune transcriptomic signatures associated with continuous M070 risk",
    subtitle = "Signature score = mean gene-wise z-scored expression; not a direct cell-count estimate",
    x = "Spearman rho", y = NULL, color = "Direction", size = "-log10(FDR)"
  ) + theme_pub
ggplot2::ggsave(
  file.path(out_dir, "Figures", "04_Immune_M070_risk_correlations.png"),
  p_imm, width = 10.5, height = 6.2, dpi = 320, bg = "white"
)

## Focus-gene correlation plot (descriptive bridge only).
fg <- focus_all[focus_all$Available, , drop = FALSE]
fg$Dataset <- factor(fg$Dataset, levels = c("E-MTAB-179", "SEQC"))
fg$Direction <- ifelse(fg$Rho_M070 >= 0, "Positive", "Negative")
p_fg <- ggplot2::ggplot(
  fg,
  ggplot2::aes(x = Rho_M070, y = Gene, color = Direction)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = Rho_M070, yend = Gene),
    color = "grey70", linewidth = 0.7
  ) +
  ggplot2::geom_point(ggplot2::aes(size = -log10(pmax(SpearmanFDR, 1e-300)))) +
  ggplot2::facet_wrap(~Dataset, nrow = 1) +
  ggplot2::scale_color_manual(values = c("Negative" = "#2166AC", "Positive" = "#B2182B")) +
  ggplot2::labs(
    title = "PRDX4/MIF-axis genes and continuous M070 risk",
    subtitle = "Descriptive bulk-transcriptomic bridge; PRDX4 is a component of M070",
    x = "Spearman rho", y = NULL, color = "Direction", size = "-log10(FDR)"
  ) + theme_pub
ggplot2::ggsave(
  file.path(out_dir, "Figures", "05_PRDX4_MIF_axis_genes_vs_M070.png"),
  p_fg, width = 9, height = 4.7, dpi = 320, bg = "white"
)

## PRDX4-to-program heatmap.
pp <- rbind(prdx4_prog_emtab, prdx4_prog_seqc)
pp$Dataset <- factor(pp$Dataset, levels = c("E-MTAB-179", "SEQC"))
pp$ProgramLabel <- ifelse(pp$ProgramType == "Hallmark", pretty_name(pp$Program), pp$Program)
p_pp <- ggplot2::ggplot(
  pp,
  ggplot2::aes(x = Dataset, y = ProgramLabel, fill = Rho_PRDX4)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Rho_PRDX4)), size = 3.2) +
  ggplot2::facet_grid(ProgramType ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  ggplot2::labs(
    title = "PRDX4 expression and selected metabolic-redox/immune programs",
    subtitle = "Descriptive correlation analysis in each bulk cohort",
    x = NULL, y = NULL, fill = "Spearman rho"
  ) + theme_pub
ggplot2::ggsave(
  file.path(out_dir, "Figures", "06_PRDX4_selected_program_correlations.png"),
  p_pp, width = 8.5, height = 8.3, dpi = 320, bg = "white"
)

## ----------------------------- final audit summary ----------------------------
core_hall <- hall_cross_grp[hall_cross_grp$GeneSet %in% focus_hallmark_patterns, , drop = FALSE]
core_immune <- immune_cross[immune_cross$Signature %in% focus_immune_names, , drop = FALSE]
focus_mif_em <- focus_emtab[focus_emtab$Gene == "MIF", , drop = FALSE]
focus_mif_sq <- focus_seqc[focus_seqc$Gene == "MIF", , drop = FALSE]

final_summary <- data.frame(
  Item = c(
    "Frozen M070 MD5",
    "M070 refitted in 07B",
    "Features reselected in 07B",
    "Cutoff retuned in 07B",
    "E-MTAB-179 N",
    "SEQC N",
    "Hallmark gene sets tested in E-MTAB-179",
    "Hallmark gene sets tested in SEQC",
    "Hallmark FDR<0.05 and direction-concordant in both cohorts",
    "Core Hallmark programs direction-concordant",
    "KEGG collection used",
    "KEGG FDR<0.05 and direction-concordant in both cohorts",
    "Immune signatures planned",
    "Immune signatures evaluable with identical genes in both cohorts",
    "Immune signatures recorded as not evaluable",
    "Immune signatures with concordant continuous-M070 direction",
    "Core immune signatures direction-concordant",
    "MIF available in both cohorts",
    "Mechanistic interpretation level",
    "Recommended manuscript role"
  ),
  Value = c(
    model_md5,
    "FALSE", "FALSE", "FALSE",
    ncol(emtab$expr), ncol(seqc$expr),
    nrow(hall_emtab), nrow(hall_seqc),
    sum(hall_cross_grp$DirectionConcordant & hall_cross_grp$FDRlt005_Both, na.rm = TRUE),
    sum(core_hall$DirectionConcordant, na.rm = TRUE),
    kegg_subcollection,
    sum(kegg_cross_grp$DirectionConcordant & kegg_cross_grp$FDRlt005_Both, na.rm = TRUE),
    length(immune_sets),
    length(shared_immune_sets),
    length(not_evaluable),
    sum(immune_cross$RhoDirectionConcordant, na.rm = TRUE),
    sum(core_immune$RhoDirectionConcordant, na.rm = TRUE),
    as.character(nrow(focus_mif_em) == 1 && focus_mif_em$Available && nrow(focus_mif_sq) == 1 && focus_mif_sq$Available),
    "Association only; not causal proof",
    "Refresh M050-dependent pathway/immune results with locked M070; retain single-cell CellChat and virtual-knockout results only if they were not defined by M050"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  final_summary,
  file.path(out_dir, "Tables", "26_FINAL_07B_AUDIT_SUMMARY.csv"),
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
    "SC28 RECALC 07B FIX3 completed.",
    paste0("Frozen M070 MD5: ", model_md5),
    "M070 was not refitted; genes were not reselected; cutoff was not retuned.",
    paste0("E-MTAB-179 N: ", ncol(emtab$expr), "; SEQC N: ", ncol(seqc$expr)),
    paste0("Hallmark sets analysed: E-MTAB-179 ", nrow(hall_emtab), "; SEQC ", nrow(hall_seqc)),
    paste0("KEGG collection: ", kegg_subcollection),
    paste0("Immune signatures planned: ", length(immune_sets)),
    paste0("Immune signatures analysed with identical shared genes: ", length(shared_immune_sets)),
    paste0("Immune signatures not evaluable (<3 shared nonconstant genes): ", length(not_evaluable)),
    paste0("Warnings captured: ", nrow(warning_log)),
    "Interpret pathway/immune/PRDX4-MIF results as transcriptomic associations, not causal mechanisms.",
    "Upload the ENTIRE SC28_RECALC_07B_FIX3_M070_Pathway_ImmuneMechanism_20260809 folder as ZIP for audit."
  ),
  file.path(out_dir, "Logs", "03_run_summary.txt")
)

cat("\n============================================================\n")
cat("SC28 RECALC 07B FIX3 completed.\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("Frozen M070 MD5: ", model_md5, "\n", sep = "")
cat("Warnings captured: ", nrow(warning_log), "\n", sep = "")
cat("Upload the ENTIRE output folder as ZIP for audit.\n")
cat("============================================================\n")
