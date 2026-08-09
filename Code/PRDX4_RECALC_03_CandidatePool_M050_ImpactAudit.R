rm(list = ls())
options(stringsAsFactors = FALSE)

## ============================================================
## PRDX4_RECALC_03_CandidatePool_M050_ImpactAudit.R
##
## Purpose
## 1. Audit the historical SC12A candidate panel actually supplied to 101ML.
## 2. Determine whether PRDX4/LDHA were already present BEFORE the historical
##    code appended them to candidate_genes_query.
## 3. Compare the historical candidate panel and locked M050 genes against the
##    corrected proliferation-balanced single-cell pseudobulk results.
## 4. Produce the evidence needed to decide how the 101ML step should be rerun.
##
## This script does NOT refit or select a prognostic model.
## It is an audit step before the clean 101ML rerun.
## ============================================================

root_dir <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"

phase2_dir <- file.path(
  root_dir,
  "SC28_RECALC_02_ProliferationBalanced_Pseudobulk_20260806"
)

sc12a_panel_rds <- file.path(
  root_dir,
  "SC12A_PRDX4_CandidatePool_Refinement_ModelInput",
  "Tables",
  "07_SC12A_SEQC_model_input_PRDX4_panels.rds"
)

emtab_expr_file <- file.path(
  root_dir,
  "SC13A_MANUAL_V6_EMTAB179_ADF_TabFix",
  "Tables",
  "05_EMTAB179_expr_gene_symbol_matched.rds"
)

seqc_ready_file <- file.path(
  root_dir,
  "01_Preprocess",
  "SEQC_GSE49711_ready_FIX2",
  "SEQC_external_validation_ready_input_FIX2.rds"
)

m050_coef_file <- file.path(
  root_dir,
  "SC18A_RiskModel_StandardPrognosticFigures",
  "Tables",
  "00_M050_locked_coefficients.csv"
)

historical_sc13e_coverage <- file.path(
  root_dir,
  "SC13E_V2_Standard101ML_EMTAB179_Train_SEQC_Validation",
  "Tables",
  "02_candidate_gene_coverage.csv"
)

historical_sc13e_mapping <- file.path(
  root_dir,
  "SC13E_V2_Standard101ML_EMTAB179_Train_SEQC_Validation",
  "Tables",
  "02B_candidate_gene_train_valid_mapping_detail.csv"
)

historical_sc13e_ranking <- file.path(
  root_dir,
  "SC13E_V2_Standard101ML_EMTAB179_Train_SEQC_Validation",
  "Tables",
  "12_SEQC_validation_model_ranking.csv"
)

out_dir <- file.path(
  root_dir,
  "SC28_RECALC_03_CandidatePool_M050_ImpactAudit_20260806"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "Logs", "SC28_RECALC_03_log.txt")
sink(log_file, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("SC28 RECALC 03: Candidate-pool and M050 impact audit\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("============================================================\n\n")

## ============================================================
## Helpers
## ============================================================

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

norm_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

detect_gene_col <- function(df) {
  hit <- intersect(
    c("Gene", "gene", "Symbol", "symbol", "GeneSymbol", "gene_symbol"),
    colnames(df)
  )
  if (length(hit) > 0) return(hit[1])
  colnames(df)[1]
}

get_expr_gene_names <- function(x) {
  if (is.list(x) && !is.null(x$expr)) x <- x$expr
  if (is.data.frame(x) || is.matrix(x)) {
    rn <- rownames(x)
    if (!is.null(rn) && any(nzchar(rn))) return(norm_gene(rn))
    first <- x[[1]]
    if (!is.numeric(first)) return(norm_gene(first))
  }
  stop("Unable to determine expression gene names from object.")
}

add_corrected_single_cell_status <- function(gene_df, deg_all, ranked_all, ranked_sig) {
  gene_df$GeneKey <- norm_gene(gene_df$Gene)

  deg2 <- deg_all
  deg2$GeneKey <- norm_gene(deg2$Gene)
  deg2 <- deg2[!duplicated(deg2$GeneKey), , drop = FALSE]

  cols_need <- intersect(
    c(
      "GeneKey", "logFC", "logCPM", "PValue_edgeR", "FDR_edgeR",
      "Is_CellCycle_Proliferation", "Direction"
    ),
    colnames(deg2)
  )
  out <- merge(
    gene_df,
    deg2[, cols_need, drop = FALSE],
    by = "GeneKey",
    all.x = TRUE,
    sort = FALSE
  )

  r1 <- ranked_all
  r1$GeneKey <- norm_gene(r1$Gene)
  rank1_col <- "Rank_NonProliferative_HighUp"
  if (rank1_col %in% colnames(r1)) {
    rank_map1 <- setNames(r1[[rank1_col]], r1$GeneKey)
    out[[rank1_col]] <- unname(rank_map1[out$GeneKey])
  }

  r2 <- ranked_sig
  r2$GeneKey <- norm_gene(r2$Gene)
  rank2_col <- "Rank_Significant_NonProliferative_HighUp"
  if (rank2_col %in% colnames(r2)) {
    rank_map2 <- setNames(r2[[rank2_col]], r2$GeneKey)
    out[[rank2_col]] <- unname(rank_map2[out$GeneKey])
  }

  out$Corrected_Significant_NonProliferative_HighUp <-
    !is.na(out$FDR_edgeR) &
    out$FDR_edgeR < 0.05 &
    !is.na(out$logFC) &
    out$logFC > 0 &
    !is.na(out$Is_CellCycle_Proliferation) &
    !as.logical(out$Is_CellCycle_Proliferation)

  out
}

## ============================================================
## 1. Check required corrected Phase-2 outputs
## ============================================================

deg_all_file <- file.path(phase2_dir, "Tables", "10_paired_edgeR_DEG_all.csv")
rank_all_file <- file.path(phase2_dir, "Tables", "11_nonproliferative_HighUp_all_ranked.csv")
rank_sig_file <- file.path(phase2_dir, "Tables", "12_nonproliferative_HighUp_significant_ranked.csv")
prdx4_file <- file.path(phase2_dir, "Tables", "14_PRDX4_rank_and_inference_summary.csv")

required_files <- c(
  deg_all_file,
  rank_all_file,
  rank_sig_file,
  prdx4_file,
  sc12a_panel_rds,
  emtab_expr_file,
  seqc_ready_file,
  m050_coef_file
)

missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required) > 0) {
  cat("Missing required files:\n")
  cat(paste0("  - ", missing_required, collapse = "\n"), "\n")
  stop("Required inputs are missing. See log for exact paths.")
}

deg_all <- read.csv(deg_all_file, check.names = FALSE, stringsAsFactors = FALSE)
rank_all <- read.csv(rank_all_file, check.names = FALSE, stringsAsFactors = FALSE)
rank_sig <- read.csv(rank_sig_file, check.names = FALSE, stringsAsFactors = FALSE)
prdx4_summary <- read.csv(prdx4_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Corrected Phase-2 candidates:\n")
cat("  All non-proliferative High-up:", nrow(rank_all), "\n")
cat("  Significant non-proliferative High-up:", nrow(rank_sig), "\n\n")

## ============================================================
## 2. Recover the historical SC12A panel BEFORE PRDX4/LDHA append
## ============================================================

sc12a_obj <- readRDS(sc12a_panel_rds)
if (is.null(sc12a_obj$panels) || !is.list(sc12a_obj$panels)) {
  stop("SC12A RDS does not contain a $panels list.")
}

panels <- sc12a_obj$panels
panel_names <- names(panels)
if (is.null(panel_names) || length(panel_names) == 0) {
  stop("SC12A $panels is unnamed; cannot reproduce historical panel choice.")
}

if ("BALANCED_40" %in% panel_names) {
  historical_panel_name <- "BALANCED_40"
} else {
  hit <- grep("BALANCED", panel_names, value = TRUE, ignore.case = TRUE)
  if (length(hit) > 0) {
    historical_panel_name <- hit[1]
  } else {
    historical_panel_name <- panel_names[1]
  }
}

panel_before_append <- unique(norm_gene(panels[[historical_panel_name]]))
panel_before_append <- panel_before_append[nzchar(panel_before_append) & !is.na(panel_before_append)]

prdx4_before_append <- "PRDX4" %in% panel_before_append
ldha_before_append <- "LDHA" %in% panel_before_append

## Reproduce the exact historical SC13E line:
## candidate_genes_query <- unique(c(candidate_genes_query, "PRDX4", "LDHA"))
panel_after_append <- unique(c(panel_before_append, "PRDX4", "LDHA"))

panel_inventory <- data.frame(
  Panel = panel_names,
  GeneN = vapply(panels, function(x) length(unique(norm_gene(x))), integer(1)),
  ContainsPRDX4 = vapply(panels, function(x) "PRDX4" %in% norm_gene(x), logical(1)),
  ContainsLDHA = vapply(panels, function(x) "LDHA" %in% norm_gene(x), logical(1)),
  stringsAsFactors = FALSE
)
write.csv(
  panel_inventory,
  file.path(out_dir, "Tables", "01_SC12A_panel_inventory.csv"),
  row.names = FALSE
)

before_after <- data.frame(
  Item = c(
    "Historical panel name",
    "Panel genes before explicit append",
    "PRDX4 present before explicit append",
    "LDHA present before explicit append",
    "Panel genes after unique(c(panel, PRDX4, LDHA))",
    "PRDX4 was newly added by explicit append",
    "LDHA was newly added by explicit append"
  ),
  Value = c(
    historical_panel_name,
    length(panel_before_append),
    prdx4_before_append,
    ldha_before_append,
    length(panel_after_append),
    !prdx4_before_append,
    !ldha_before_append
  ),
  stringsAsFactors = FALSE
)
write.csv(
  before_after,
  file.path(out_dir, "Tables", "02_historical_candidate_append_audit.csv"),
  row.names = FALSE
)

cat("Historical candidate-panel audit:\n")
print(before_after)
cat("\n")

## ============================================================
## 3. Reconstruct bulk cross-platform availability
## ============================================================

emtab_expr <- readRDS(emtab_expr_file)
seqc_obj <- readRDS(seqc_ready_file)

emtab_genes <- unique(get_expr_gene_names(emtab_expr))
seqc_genes <- unique(get_expr_gene_names(seqc_obj))

panel_map <- data.frame(
  Gene = panel_after_append,
  Present_EMTAB179 = panel_after_append %in% emtab_genes,
  Present_SEQC = panel_after_append %in% seqc_genes,
  stringsAsFactors = FALSE
)
panel_map$Present_Both <- panel_map$Present_EMTAB179 & panel_map$Present_SEQC

write.csv(
  panel_map,
  file.path(out_dir, "Tables", "03_historical_panel_bulk_availability.csv"),
  row.names = FALSE
)

cat("Historical panel after explicit append:\n")
cat("  Genes:", nrow(panel_map), "\n")
cat("  Available in both E-MTAB-179 and SEQC:", sum(panel_map$Present_Both), "\n\n")

## ============================================================
## 4. Compare historical panel with corrected single-cell results
## ============================================================

panel_status <- add_corrected_single_cell_status(
  data.frame(Gene = panel_after_append, stringsAsFactors = FALSE),
  deg_all,
  rank_all,
  rank_sig
)
panel_status <- merge(
  panel_status,
  panel_map,
  by.x = "GeneKey",
  by.y = "Gene",
  all.x = TRUE,
  sort = FALSE
)

## Restore a clean display gene column.
if ("Gene.x" %in% colnames(panel_status)) {
  panel_status$Gene <- panel_status$Gene.x
} else if (!"Gene" %in% colnames(panel_status)) {
  panel_status$Gene <- panel_status$GeneKey
}

front_cols <- c(
  "Gene", "Corrected_Significant_NonProliferative_HighUp", "logFC",
  "PValue_edgeR", "FDR_edgeR", "Direction",
  "Rank_NonProliferative_HighUp",
  "Rank_Significant_NonProliferative_HighUp",
  "Present_EMTAB179", "Present_SEQC", "Present_Both"
)
front_cols <- intersect(front_cols, colnames(panel_status))
panel_status <- panel_status[, c(front_cols, setdiff(colnames(panel_status), front_cols)), drop = FALSE]

write.csv(
  panel_status,
  file.path(out_dir, "Tables", "04_historical_panel_vs_corrected_singlecell.csv"),
  row.names = FALSE
)

## ============================================================
## 5. Audit locked M050 genes against corrected single-cell results
## ============================================================

m050 <- read.csv(m050_coef_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- detect_gene_col(m050)
m050$Gene <- norm_gene(m050[[gene_col]])
m050 <- m050[!duplicated(m050$Gene), , drop = FALSE]

m050_status <- add_corrected_single_cell_status(
  data.frame(Gene = m050$Gene, stringsAsFactors = FALSE),
  deg_all,
  rank_all,
  rank_sig
)

coef_cols <- setdiff(colnames(m050), "Gene")
coef_keep <- m050[, c("Gene", coef_cols), drop = FALSE]
coef_keep$GeneKey <- norm_gene(coef_keep$Gene)

m050_status <- merge(
  m050_status,
  coef_keep[, setdiff(colnames(coef_keep), "Gene"), drop = FALSE],
  by = "GeneKey",
  all.x = TRUE,
  sort = FALSE
)

write.csv(
  m050_status,
  file.path(out_dir, "Tables", "05_M050_20genes_vs_corrected_singlecell.csv"),
  row.names = FALSE
)

m050_n <- nrow(m050_status)
m050_sig_n <- sum(m050_status$Corrected_Significant_NonProliferative_HighUp, na.rm = TRUE)

## ============================================================
## 6. Historical SC13E files, if available
## ============================================================

historical_files <- data.frame(
  Item = c("SC13E candidate coverage", "SC13E candidate mapping", "SC13E SEQC ranking"),
  Path = c(historical_sc13e_coverage, historical_sc13e_mapping, historical_sc13e_ranking),
  Exists = file.exists(c(historical_sc13e_coverage, historical_sc13e_mapping, historical_sc13e_ranking)),
  stringsAsFactors = FALSE
)
write.csv(
  historical_files,
  file.path(out_dir, "Tables", "06_historical_SC13E_file_availability.csv"),
  row.names = FALSE
)

if (file.exists(historical_sc13e_coverage)) {
  file.copy(
    historical_sc13e_coverage,
    file.path(out_dir, "Tables", "07_original_SC13E_candidate_gene_coverage.csv"),
    overwrite = TRUE
  )
}
if (file.exists(historical_sc13e_mapping)) {
  file.copy(
    historical_sc13e_mapping,
    file.path(out_dir, "Tables", "08_original_SC13E_candidate_mapping.csv"),
    overwrite = TRUE
  )
}
if (file.exists(historical_sc13e_ranking)) {
  file.copy(
    historical_sc13e_ranking,
    file.path(out_dir, "Tables", "09_original_SC13E_SEQC_model_ranking.csv"),
    overwrite = TRUE
  )
}

## ============================================================
## 7. Decision summary for the next clean 101ML rerun
## ============================================================

prdx4_row <- deg_all[norm_gene(deg_all$Gene) == "PRDX4", , drop = FALSE]
if (nrow(prdx4_row) == 0) stop("PRDX4 missing from corrected Phase-2 DEG table.")

prdx4_rank_all <- NA_real_
hit <- rank_all[norm_gene(rank_all$Gene) == "PRDX4", , drop = FALSE]
if (nrow(hit) > 0 && "Rank_NonProliferative_HighUp" %in% colnames(hit)) {
  prdx4_rank_all <- hit$Rank_NonProliferative_HighUp[1]
}

prdx4_rank_sig <- NA_real_
hit2 <- rank_sig[norm_gene(rank_sig$Gene) == "PRDX4", , drop = FALSE]
if (nrow(hit2) > 0 && "Rank_Significant_NonProliferative_HighUp" %in% colnames(hit2)) {
  prdx4_rank_sig <- hit2$Rank_Significant_NonProliferative_HighUp[1]
}

decision <- data.frame(
  Item = c(
    "Corrected significant non-proliferative High-up genes",
    "PRDX4 corrected logFC",
    "PRDX4 corrected edgeR P",
    "PRDX4 corrected edgeR FDR",
    "PRDX4 corrected rank among all non-proliferative High-up",
    "PRDX4 corrected rank among significant non-proliferative High-up",
    "Historical selected SC12A panel",
    "Historical panel size before explicit append",
    "PRDX4 present in historical panel before explicit append",
    "LDHA present in historical panel before explicit append",
    "PRDX4 newly added by explicit append",
    "LDHA newly added by explicit append",
    "Historical panel genes available in both bulk cohorts",
    "Locked M050 gene count",
    "Locked M050 genes supported by corrected significant non-proliferative High-up criterion",
    "SEQC was historically used for final 101ML ranking",
    "Clean 101ML rerun needed before calling SEQC an external validation cohort"
  ),
  Value = c(
    nrow(rank_sig),
    prdx4_row$logFC[1],
    prdx4_row$PValue_edgeR[1],
    prdx4_row$FDR_edgeR[1],
    prdx4_rank_all,
    prdx4_rank_sig,
    historical_panel_name,
    length(panel_before_append),
    prdx4_before_append,
    ldha_before_append,
    !prdx4_before_append,
    !ldha_before_append,
    sum(panel_map$Present_Both),
    m050_n,
    m050_sig_n,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  decision,
  file.path(out_dir, "Tables", "10_DECISION_SUMMARY.csv"),
  row.names = FALSE
)

cat("Decision summary:\n")
print(decision)
cat("\n")

## Plain-language guardrails for manuscript revision.
guardrails <- c(
  "1. The corrected proliferation-balanced pseudobulk result is the authoritative single-cell candidate result.",
  "2. Do not retain the historical statement that PRDX4 ranked 14th after proliferation matching if the corrected rank differs.",
  "3. If PRDX4 was absent from the SC12A panel before the explicit append, do not state that PRDX4 entered the 101ML candidate pool without explicit inclusion.",
  "4. Even if PRDX4 was already present in SC12A, the historical SC13E code ranked final models using SEQC; therefore SEQC cannot simultaneously be described as an untouched independent external validation cohort for that historical model-selection procedure.",
  "5. To use SEQC as external validation, the next 101ML rerun must develop/select/lock the model using E-MTAB-179 only, and evaluate SEQC only after locking the model.",
  "6. Do not force PRDX4 or LDHA into the corrected 101ML candidate set; use an explicit objective candidate rule and report it transparently."
)
writeLines(guardrails, file.path(out_dir, "Tables", "11_MANUSCRIPT_GUARDRAILS.txt"))

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "Logs", "sessionInfo.txt")
)

cat("============================================================\n")
cat("SC28 RECALC 03 completed successfully.\n")
cat("Please upload the ENTIRE output folder as a ZIP.\n")
cat("Most important files:\n")
cat("  Tables/02_historical_candidate_append_audit.csv\n")
cat("  Tables/04_historical_panel_vs_corrected_singlecell.csv\n")
cat("  Tables/05_M050_20genes_vs_corrected_singlecell.csv\n")
cat("  Tables/10_DECISION_SUMMARY.csv\n")
cat("  Tables/11_MANUSCRIPT_GUARDRAILS.txt\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("============================================================\n")

sink()

