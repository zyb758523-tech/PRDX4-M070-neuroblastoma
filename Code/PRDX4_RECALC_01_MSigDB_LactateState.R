#!/usr/bin/env Rscript

# PRDX4 / neuroblastoma single-cell repair analysis
# Phase 1-2: rebuild the MSigDB lactate gene set, fix feature-ID handling,
#            calculate LactateScore, and define exact within-sample states.
#
# This script deliberately does NOT use SC28A/FIX2 LactateState labels.
# Historical outputs are left untouched.

options(stringsAsFactors = FALSE)

project_root <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
seurat_path <- "C:/Users/Administrator/Desktop/NB_Lactate_Bulk_ServerInput_20260522_112238/GSE137804_seurat_final.rds"
annotation_path <- file.path(
  project_root,
  "SC27C_FIXED_RevisedManualAnnotation_AfterOriginalComparison",
  "Tables",
  "02_SC27C_cell_level_revised_annotation.csv"
)
out_dir <- file.path(project_root, "SC28_RECALC_01B_MSigDB_LactateState_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Rdata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "Logs", "SC28_RECALC_01_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

stop_clean <- function(...) {
  msg("ERROR: ", paste0(..., collapse = ""))
  stop(paste0(..., collapse = ""), call. = FALSE)
}

install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg("Installing CRAN package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in c("data.table", "Matrix", "msigdbr")) install_cran_if_missing(pkg)

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop_clean("Package 'SeuratObject' is required. It should already be present in the R environment used for the original Seurat analysis.")
}

# PRDX4 needs a complete human Ensembl <-> symbol map because it is not part of
# the lactate signature itself.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg("Installing Bioconductor package: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(msigdbr)
  library(SeuratObject)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

msg("============================================================")
msg("SC28_RECALC_01: MSigDB lactate state reconstruction")
msg("R version: ", R.version.string)
msg("msigdbr version: ", as.character(utils::packageVersion("msigdbr")))
msg("Output: ", out_dir)
msg("============================================================")

if (!file.exists(seurat_path)) stop_clean("Missing Seurat object: ", seurat_path)
if (!file.exists(annotation_path)) stop_clean("Missing SC27C annotation: ", annotation_path)

# -----------------------------------------------------------------------------
# 1. Rebuild the five MSigDB lactate-related gene sets from MSigDB 2026.1.Hs.
#    Stable GO/HPO exact-source IDs are used to survive gene-set renaming.
# -----------------------------------------------------------------------------

collection_info <- as.data.frame(msigdbr::msigdbr_collections(db_species = "HS"))
db_versions <- unique(as.character(collection_info$db_version))
msg("MSigDB version(s) exposed by msigdbr: ", paste(db_versions, collapse = ", "))

if (!any(db_versions == "2026.1.Hs")) {
  stop_clean(
    "This analysis is locked to MSigDB 2026.1.Hs, but msigdbr exposes: ",
    paste(db_versions, collapse = ", "),
    ". Install msigdbr 26.1.0 before running the analysis."
  )
}

target_sets <- data.table(
  ManuscriptName = c(
    "GOBP_LACTATE_METABOLIC_PROCESS",
    "HP_INCREASED_SERUM_LACTATE",
    "HP_LACTIC_ACIDOSIS",
    "HP_LACTIC_ACIDURIA",
    "HP_SEVERE_LACTIC_ACIDOSIS"
  ),
  ExactSource = c("GO:0006089", "HP:0002151", "HP:0003128", "HP:0003648", "HP:0004900")
)

msg("Loading MSigDB C5 collection...")
c5 <- as.data.table(msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C5"
))

required_cols <- c("gs_name", "gs_exact_source", "gene_symbol", "ensembl_gene", "ncbi_gene", "db_version")
missing_cols <- setdiff(required_cols, names(c5))
if (length(missing_cols)) stop_clean("msigdbr output is missing columns: ", paste(missing_cols, collapse = ", "))

sel <- c5[gs_exact_source %in% target_sets$ExactSource]
if (!nrow(sel)) stop_clean("None of the five stable GO/HPO source IDs were found in MSigDB C5.")

found_sources <- unique(sel$gs_exact_source)
missing_sources <- setdiff(target_sets$ExactSource, found_sources)
if (length(missing_sources)) {
  stop_clean("Missing MSigDB exact-source IDs: ", paste(missing_sources, collapse = ", "))
}

set_summary <- sel[, .(
  CurrentMSigDBName = paste(sort(unique(gs_name)), collapse = " | "),
  GeneN = uniqueN(gene_symbol),
  MSigDBVersion = paste(sort(unique(db_version)), collapse = " | ")
), by = gs_exact_source]
set_summary <- merge(target_sets, set_summary, by.x = "ExactSource", by.y = "gs_exact_source", all.x = TRUE)
set_summary[, SetOrder := match(ExactSource, target_sets$ExactSource)]
setorder(set_summary, SetOrder)
set_summary[, SetOrder := NULL]

union_symbols <- sort(unique(sel$gene_symbol[!is.na(sel$gene_symbol) & nzchar(sel$gene_symbol)]))
union_n <- length(union_symbols)

membership <- unique(sel[, .(
  MSigDBVersion = db_version,
  CurrentMSigDBName = gs_name,
  ExactSource = gs_exact_source,
  GeneSymbol = gene_symbol,
  EnsemblGene = ensembl_gene,
  EntrezGene = ncbi_gene
)])
membership <- merge(membership, target_sets, by = "ExactSource", all.x = TRUE)
setcolorder(membership, c("ManuscriptName", "CurrentMSigDBName", "ExactSource", "GeneSymbol", "EnsemblGene", "EntrezGene", "MSigDBVersion"))
setorder(membership, ManuscriptName, GeneSymbol)

fwrite(set_summary, file.path(out_dir, "Tables", "01_MSIGDB_2026.1_five_set_summary.csv"))
fwrite(membership, file.path(out_dir, "Tables", "02_MSIGDB_2026.1_five_set_membership.csv"))
fwrite(data.table(GeneSymbol = union_symbols), file.path(out_dir, "Tables", "03_MSIGDB_2026.1_lactate_union.csv"))

union_check <- data.table(
  Item = c("MSigDB version", "Five exact-source sets found", "Unique union genes", "Historical manuscript value", "Exact match to historical 284"),
  Value = c("2026.1.Hs", "5", as.character(union_n), "284", as.character(union_n == 284L))
)
fwrite(union_check, file.path(out_dir, "Tables", "04_MSIGDB_union_check.csv"))
msg("Verified five-set union size = ", union_n, " unique genes. Historical manuscript value = 284.")

# Build target identifier table. Direct HGNC symbols are preferred; Ensembl IDs
# are a fallback for layers whose row names are ENSG identifiers.
target_id_map <- unique(membership[!is.na(GeneSymbol) & nzchar(GeneSymbol), .(
  GeneSymbol,
  EnsemblGene = sub("\\..*$", "", EnsemblGene),
  EntrezGene = as.character(EntrezGene)
)])

prdx4_ens <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "PRDX4",
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL", "ENTREZID")
)
prdx4_ens <- unique(sub("\\..*$", "", prdx4_ens$ENSEMBL[!is.na(prdx4_ens$ENSEMBL)]))
prdx4_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "PRDX4",
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
prdx4_entrez <- unique(as.character(prdx4_entrez[!is.na(prdx4_entrez)]))
msg("PRDX4 Ensembl ID(s): ", paste(prdx4_ens, collapse = ", "))
msg("PRDX4 Entrez ID(s): ", paste(prdx4_entrez, collapse = ", "))

# -----------------------------------------------------------------------------
# 2. Load SC27C final annotations and the original GSE137804 Seurat object.
# -----------------------------------------------------------------------------

msg("Loading SC27C final cell-level annotation...")
ann <- fread(annotation_path)
msg("Annotation dimensions: ", nrow(ann), " x ", ncol(ann))

msg("Loading GSE137804 Seurat object. This is the slowest step (~4.5 GB)...")
obj <- readRDS(seurat_path)
msg("Seurat cells: ", ncol(obj), "; assays: ", paste(names(obj@assays), collapse = ", "))

obj_cells <- colnames(obj)

# SC27C contains several historical/manual annotation columns. Do not infer the
# final annotation column from its values: the audited SC28A code explicitly
# used RevisedCellType_SC27C.
required_sc27c_cols <- c("CellID", "Sample_sc", "RevisedCellType_SC27C")
missing_sc27c_cols <- setdiff(required_sc27c_cols, names(ann))
if (length(missing_sc27c_cols)) {
  stop_clean(
    "SC27C final annotation table is missing audited required column(s): ",
    paste(missing_sc27c_cols, collapse = ", ")
  )
}

cell_col <- "CellID"
cell_overlap <- sum(as.character(ann[[cell_col]]) %in% obj_cells)
if (cell_overlap < 0.9 * ncol(obj)) {
  stop_clean("Could not identify SC27C cell-ID column reliably. Best column = ", cell_col, "; overlap = ", cell_overlap)
}

sample_col <- "Sample_sc"
annotation_col <- "RevisedCellType_SC27C"
if (!any(as.character(ann[[annotation_col]]) == "Tumor-like neural", na.rm = TRUE)) {
  stop_clean("The audited final SC27C column RevisedCellType_SC27C does not contain 'Tumor-like neural'.")
}

msg("Detected columns: CellID = ", cell_col, "; Sample = ", sample_col, "; Revised annotation = ", annotation_col)

ann[, CellID_recalc := as.character(get(cell_col))]
ann[, Sample_recalc := as.character(get(sample_col))]
ann[, Annotation_recalc := as.character(get(annotation_col))]

# Main analysis: core Tumor-like neural cells from primary-tumor samples only.
# GSE137804 T-prefix samples are primary neuroblastomas; F-prefix samples are not
# included in the main lactate-state discovery analysis.
core <- ann[
  Annotation_recalc == "Tumor-like neural" & grepl("_T[0-9]+$", Sample_recalc)
]
core <- core[CellID_recalc %in% obj_cells]

tumor_samples <- sort(unique(core$Sample_recalc))
if (length(tumor_samples) != 16L) {
  stop_clean("Expected 16 primary-tumor T samples, but detected ", length(tumor_samples), ": ", paste(tumor_samples, collapse = ", "))
}

sample_counts <- core[, .(CoreTumorLikeCellN = .N), by = Sample_recalc]
setnames(sample_counts, "Sample_recalc", "Sample_sc")
setorder(sample_counts, Sample_sc)
fwrite(sample_counts, file.path(out_dir, "Tables", "05_core_tumorlike_tumor_sample_counts.csv"))
msg("Core Tumor-like neural cells in 16 T samples: ", nrow(core))

# -----------------------------------------------------------------------------
# 3. Diagnose layer-specific feature IDs and target-gene coverage.
# -----------------------------------------------------------------------------

rna_assay <- obj[["RNA"]]
layer_names <- SeuratObject::Layers(rna_assay)
msg("RNA layers: ", length(layer_names))

find_counts_layer <- function(sample_name) {
  exact <- paste0("counts.", sample_name)
  if (exact %in% layer_names) return(exact)
  hits <- layer_names[grepl("^counts", layer_names) & endsWith(layer_names, sample_name)]
  if (length(hits) == 1L) return(hits)
  NA_character_
}

map_target_features <- function(features) {
  clean <- sub("\\..*$", "", features)
  direct_idx <- which(features %in% union_symbols)
  direct <- data.table(RowIndex = direct_idx, Feature = features[direct_idx], GeneSymbol = features[direct_idx], Mapping = "HGNC_symbol")

  ens_lookup <- unique(target_id_map[, .(EnsemblGene, GeneSymbol)])
  ens_lookup <- ens_lookup[!is.na(EnsemblGene) & nzchar(EnsemblGene)]
  ens_dt <- data.table(RowIndex = seq_along(features), Feature = features, EnsemblGene = clean)
  ens_dt <- merge(ens_dt, ens_lookup, by = "EnsemblGene", allow.cartesian = TRUE)
  if (nrow(ens_dt)) ens_dt[, Mapping := "Ensembl_to_symbol"]
  if (nrow(ens_dt)) ens_dt[, EnsemblGene := NULL]

  entrez_lookup <- unique(target_id_map[, .(EntrezGene, GeneSymbol)])
  entrez_lookup <- entrez_lookup[!is.na(EntrezGene) & nzchar(EntrezGene)]
  entrez_dt <- data.table(RowIndex = seq_along(features), Feature = features, EntrezGene = as.character(features))
  entrez_dt <- merge(entrez_dt, entrez_lookup, by = "EntrezGene", allow.cartesian = TRUE)
  if (nrow(entrez_dt)) entrez_dt[, Mapping := "Entrez_to_symbol"]
  if (nrow(entrez_dt)) entrez_dt[, EntrezGene := NULL]

  out <- rbindlist(list(direct, ens_dt, entrez_dt), use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(out)
  # If a direct symbol mapping exists for a gene, use it in preference to a
  # second identifier representation. Otherwise use Ensembl mapping.
  out[, MapPriority := fcase(
    Mapping == "HGNC_symbol", 1L,
    Mapping == "Ensembl_to_symbol", 2L,
    default = 3L
  )]
  setorder(out, GeneSymbol, MapPriority, RowIndex)
  out[, MapPriority := NULL]
  unique(out, by = c("RowIndex", "GeneSymbol"))
}

diagnostics <- list()
feature_examples <- list()

for (i in seq_along(tumor_samples)) {
  smp <- tumor_samples[i]
  msg(sprintf("Feature diagnostic %02d/%02d: %s", i, length(tumor_samples), smp))
  ly <- find_counts_layer(smp)
  if (is.na(ly)) stop_clean("No unique RNA counts layer found for sample: ", smp)
  mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  features <- rownames(mat)
  m <- map_target_features(features)
  available <- sort(unique(m$GeneSymbol))
  direct_n <- uniqueN(m[Mapping == "HGNC_symbol", GeneSymbol])
  ens_n <- uniqueN(m[Mapping == "Ensembl_to_symbol", GeneSymbol])
  entrez_n <- uniqueN(m[Mapping == "Entrez_to_symbol", GeneSymbol])
  clean_features <- sub("\\..*$", "", features)
  prdx4_mode <- if ("PRDX4" %in% features) {
    "HGNC_symbol"
  } else if (any(clean_features %in% prdx4_ens)) {
    "Ensembl_to_symbol"
  } else if (any(as.character(features) %in% prdx4_entrez)) {
    "Entrez_to_symbol"
  } else {
    "NOT_FOUND"
  }
  feature_type <- if (mean(grepl("^ENSG[0-9]+", clean_features)) > 0.5) {
    "predominantly_Ensembl"
  } else if (mean(grepl("^[0-9]+$", features)) > 0.5) {
    "predominantly_Entrez"
  } else {
    "predominantly_symbol_or_other"
  }

  diagnostics[[smp]] <- data.table(
    Sample_sc = smp,
    Layer = ly,
    FeatureN = length(features),
    LayerCellN = ncol(mat),
    FeatureType = feature_type,
    LactateUnionN = union_n,
    AvailableLactateGeneN = length(available),
    DirectSymbolTargetN = direct_n,
    EnsemblMappedTargetN = ens_n,
    EntrezMappedTargetN = entrez_n,
    CoverageFraction = length(available) / union_n,
    PRDX4Mapping = prdx4_mode
  )
  feature_examples[[smp]] <- data.table(
    Sample_sc = smp,
    FeatureExample = head(features, 20L)
  )
  rm(mat)
  gc(verbose = FALSE)
}

diag_dt <- rbindlist(diagnostics)
examples_dt <- rbindlist(feature_examples)
fwrite(diag_dt, file.path(out_dir, "Tables", "06_RNA_layer_geneID_diagnostic.csv"))
fwrite(examples_dt, file.path(out_dir, "Tables", "07_RNA_layer_feature_examples.csv"))

if (any(diag_dt$AvailableLactateGeneN == 0L)) {
  stop_clean("At least one T-sample layer still has zero usable lactate genes after symbol/Ensembl mapping. See 06_RNA_layer_geneID_diagnostic.csv.")
}
if (any(diag_dt$PRDX4Mapping == "NOT_FOUND")) {
  stop_clean("PRDX4 is not recoverable in at least one T-sample layer. See 06_RNA_layer_geneID_diagnostic.csv.")
}

# Use a common target-gene universe across all 16 tumor samples so every cell is
# scored from the same genes. This is essential for cross-sample comparability.
available_by_sample <- vector("list", length(tumor_samples))
names(available_by_sample) <- tumor_samples
for (smp in tumor_samples) {
  ly <- find_counts_layer(smp)
  mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  available_by_sample[[smp]] <- sort(unique(map_target_features(rownames(mat))$GeneSymbol))
  rm(mat)
}
common_lactate_genes <- Reduce(intersect, available_by_sample)
common_n <- length(common_lactate_genes)
fwrite(data.table(GeneSymbol = common_lactate_genes), file.path(out_dir, "Tables", "08_common_lactate_genes_all_16_tumors.csv"))
msg("Common lactate genes across all 16 T samples: ", common_n, " / ", union_n)
if (common_n < 20L) stop_clean("Too few common lactate genes across the 16 tumor samples: ", common_n)

# -----------------------------------------------------------------------------
# 4. Recalculate LactateScore from raw RNA counts.
#    Primary score = mean of gene-wise z-scored log-normalized expression.
#    Each gene therefore contributes on the same standardized scale.
# -----------------------------------------------------------------------------

extract_collapsed_counts <- function(mat, genes, mapping_dt) {
  mp <- mapping_dt[GeneSymbol %in% genes]
  if (!nrow(mp)) stop_clean("No requested genes could be mapped in a layer.")
  # Keep one mapping row for each feature/gene pair.
  mp <- unique(mp[, .(RowIndex, GeneSymbol)])
  gene_i <- match(mp$GeneSymbol, genes)
  feature_j <- mp$RowIndex
  collapse_matrix <- Matrix::sparseMatrix(
    i = gene_i,
    j = feature_j,
    x = 1,
    dims = c(length(genes), nrow(mat)),
    dimnames = list(genes, rownames(mat))
  )
  collapse_matrix %*% mat
}

score_one_sample <- function(smp) {
  ly <- find_counts_layer(smp)
  mat_all <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  needed_cells <- core[Sample_recalc == smp, CellID_recalc]
  idx <- match(needed_cells, colnames(mat_all))
  if (anyNA(idx)) {
    stop_clean("Cell-ID mismatch inside layer for ", smp, ": ", sum(is.na(idx)), " core cells missing.")
  }
  mat <- mat_all[, idx, drop = FALSE]
  libsize <- Matrix::colSums(mat)
  if (any(libsize <= 0)) stop_clean("Zero-library cells detected in ", smp)

  target_map <- map_target_features(rownames(mat))
  target_counts <- extract_collapsed_counts(mat, common_lactate_genes, target_map)
  lognorm <- target_counts %*% Matrix::Diagonal(x = 10000 / libsize)
  lognorm@x <- log1p(lognorm@x)
  lognorm_dense <- as.matrix(lognorm)
  gene_mean <- rowMeans(lognorm_dense)
  gene_sd <- apply(lognorm_dense, 1, stats::sd)
  good <- is.finite(gene_sd) & gene_sd > 0
  if (sum(good) < 20L) stop_clean("Fewer than 20 variable lactate genes in ", smp)
  z <- sweep(lognorm_dense[good, , drop = FALSE], 1, gene_mean[good], "-")
  z <- sweep(z, 1, gene_sd[good], "/")
  score <- colMeans(z)
  score_z <- as.numeric(scale(score))

  # PRDX4 extraction from either symbol or Ensembl row names.
  clean_features <- sub("\\..*$", "", rownames(mat))
  prdx4_rows <- which(
    rownames(mat) == "PRDX4" |
      clean_features %in% prdx4_ens |
      as.character(rownames(mat)) %in% prdx4_entrez
  )
  if (!length(prdx4_rows)) stop_clean("PRDX4 could not be extracted in ", smp)
  prdx4_counts <- Matrix::colSums(mat[prdx4_rows, , drop = FALSE])
  prdx4_lognorm <- log1p(prdx4_counts / libsize * 10000)

  # Deterministic, exact, non-overlapping top/bottom 30% within each sample.
  n <- length(score)
  n_each <- floor(0.30 * n)
  ord <- order(score, needed_cells, na.last = NA)
  if (length(ord) != n) stop_clean("Non-finite LactateScore detected in ", smp)
  state <- rep("LactateMid", n)
  state[ord[seq_len(n_each)]] <- "LactateLow"
  state[ord[(n - n_each + 1L):n]] <- "LactateHigh"

  dt <- data.table(
    CellID = needed_cells,
    Sample_sc = smp,
    LactateScore = score,
    LactateScore_sampleZ = score_z,
    LactateState = state,
    PRDX4_logNorm = as.numeric(prdx4_lognorm),
    LibrarySize = as.numeric(libsize),
    LactateGenesUsedN = sum(good)
  )

  qc <- data.table(
    Sample_sc = smp,
    CoreCellN = n,
    LactateGenesCommonN = length(common_lactate_genes),
    LactateGenesVariableN = sum(good),
    ScoreSD = stats::sd(score),
    UniqueScoreN = uniqueN(signif(score, 12)),
    LowN = sum(state == "LactateLow"),
    MidN = sum(state == "LactateMid"),
    HighN = sum(state == "LactateHigh"),
    HighLowOverlapN = 0L,
    PRDX4DetectedFraction = mean(prdx4_counts > 0)
  )

  rm(mat_all, mat, target_counts, lognorm, lognorm_dense, z)
  gc(verbose = FALSE)
  list(scores = dt, qc = qc)
}

scores_list <- vector("list", length(tumor_samples))
qc_list <- vector("list", length(tumor_samples))
for (i in seq_along(tumor_samples)) {
  smp <- tumor_samples[i]
  msg(sprintf("Scoring tumor sample %02d/%02d: %s", i, length(tumor_samples), smp))
  ans <- score_one_sample(smp)
  scores_list[[i]] <- ans$scores
  qc_list[[i]] <- ans$qc
}

scores <- rbindlist(scores_list)
score_qc <- rbindlist(qc_list)
setorder(scores, Sample_sc, CellID)
setorder(score_qc, Sample_sc)

fwrite(score_qc, file.path(out_dir, "Tables", "09_LactateScore_and_exact_state_QC_by_sample.csv"))
fwrite(scores, file.path(out_dir, "Tables", "10_core_tumorlike_16tumor_LactateScore_PRDX4.csv.gz"), compress = "gzip")
saveRDS(scores, file.path(out_dir, "Rdata", "10_core_tumorlike_16tumor_LactateScore_PRDX4.rds"), compress = TRUE)

overall <- scores[, .(
  CellN = .N,
  MeanLactateScore = mean(LactateScore),
  MeanLactateScoreSampleZ = mean(LactateScore_sampleZ),
  MeanPRDX4logNorm = mean(PRDX4_logNorm),
  PRDX4DetectedFraction = mean(PRDX4_logNorm > 0)
), by = LactateState]
fwrite(overall, file.path(out_dir, "Tables", "11_LactateState_overall_summary.csv"))

# PRDX4 association is descriptive at this stage; sample-aware inference is
# intentionally deferred to the proliferation-controlled paired analysis.
rho <- suppressWarnings(stats::cor(scores$PRDX4_logNorm, scores$LactateScore, method = "spearman"))
prdx4_by_sample <- scores[, .(
  PRDX4_Mean_High = mean(PRDX4_logNorm[LactateState == "LactateHigh"]),
  PRDX4_Mean_Low = mean(PRDX4_logNorm[LactateState == "LactateLow"]),
  Delta_HighMinusLow = mean(PRDX4_logNorm[LactateState == "LactateHigh"]) - mean(PRDX4_logNorm[LactateState == "LactateLow"])
), by = Sample_sc]
fwrite(prdx4_by_sample, file.path(out_dir, "Tables", "12_PRDX4_HighLow_descriptive_by_sample.csv"))

summary_dt <- data.table(
  Item = c(
    "MSigDB_version",
    "Five_sets_found",
    "MSigDB_union_unique_genes",
    "Common_lactate_genes_across_16_tumors",
    "Primary_tumor_samples",
    "Core_tumorlike_cells",
    "LactateHigh_cells",
    "LactateLow_cells",
    "LactateMid_cells",
    "Any_sample_score_SD_zero",
    "Any_sample_only_one_unique_score",
    "PRDX4_vs_LactateScore_Spearman_rho",
    "Samples_with_PRDX4_HighMinusLow_positive"
  ),
  Value = c(
    "2026.1.Hs",
    "5",
    as.character(union_n),
    as.character(common_n),
    as.character(length(tumor_samples)),
    as.character(nrow(scores)),
    as.character(sum(scores$LactateState == "LactateHigh")),
    as.character(sum(scores$LactateState == "LactateLow")),
    as.character(sum(scores$LactateState == "LactateMid")),
    as.character(any(score_qc$ScoreSD == 0 | !is.finite(score_qc$ScoreSD))),
    as.character(any(score_qc$UniqueScoreN <= 1)),
    format(rho, digits = 8),
    paste0(sum(prdx4_by_sample$Delta_HighMinusLow > 0), "/", nrow(prdx4_by_sample))
  )
)
fwrite(summary_dt, file.path(out_dir, "Tables", "13_RECALC_01_final_summary.csv"))

writeLines(
  c(
    "SC28_RECALC_01 COMPLETED SUCCESSFULLY",
    paste0("Finished: ", Sys.time()),
    paste0("Output: ", out_dir),
    paste0("MSigDB 2026.1 five-set union: ", union_n),
    paste0("Common score genes across 16 tumors: ", common_n),
    paste0("Core Tumor-like neural cells: ", nrow(scores)),
    paste0("High/Low/Mid: ", sum(scores$LactateState == "LactateHigh"), "/", sum(scores$LactateState == "LactateLow"), "/", sum(scores$LactateState == "LactateMid")),
    paste0("Any zero-SD sample: ", any(score_qc$ScoreSD == 0 | !is.finite(score_qc$ScoreSD))),
    paste0("PRDX4 High-Low positive samples: ", sum(prdx4_by_sample$Delta_HighMinusLow > 0), "/", nrow(prdx4_by_sample))
  ),
  file.path(out_dir, "RUN_STATUS.txt")
)

msg("============================================================")
msg("SC28_RECALC_01 COMPLETED SUCCESSFULLY")
msg("MSigDB five-set union: ", union_n)
msg("Common lactate genes used in all 16 tumors: ", common_n)
msg("Core tumor cells scored: ", nrow(scores))
msg("High / Low / Mid = ", sum(scores$LactateState == "LactateHigh"), " / ", sum(scores$LactateState == "LactateLow"), " / ", sum(scores$LactateState == "LactateMid"))
msg("Any zero-SD sample: ", any(score_qc$ScoreSD == 0 | !is.finite(score_qc$ScoreSD)))
msg("PRDX4 High-Low positive samples: ", sum(prdx4_by_sample$Delta_HighMinusLow > 0), "/", nrow(prdx4_by_sample))
msg("Please zip and return this output folder for Phase 3 proliferation-controlled analysis.")
msg("============================================================")
