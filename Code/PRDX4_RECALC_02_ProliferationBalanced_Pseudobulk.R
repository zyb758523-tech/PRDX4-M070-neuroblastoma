#!/usr/bin/env Rscript

# PRDX4 / neuroblastoma single-cell repair analysis
# Phase 3: proliferation-balanced re-screening and paired pseudobulk DEG.
#
# Input is the completed SC28_RECALC_01B result. This script does NOT use the
# historical SC28A/FIX2 labels or the historical SC28B pseudobulk table.

options(stringsAsFactors = FALSE, timeout = 3600)

project_root <- "C:/Users/Administrator/Desktop/NB_Lactate_Rebuild_2026"
seurat_path <- "C:/Users/Administrator/Desktop/NB_Lactate_Bulk_ServerInput_20260522_112238/GSE137804_seurat_final.rds"
phase1_dir <- file.path(project_root, "SC28_RECALC_01B_MSigDB_LactateState_20260806")
score_rds <- file.path(phase1_dir, "Rdata", "10_core_tumorlike_16tumor_LactateScore_PRDX4.rds")

out_dir <- file.path(project_root, "SC28_RECALC_02_ProliferationBalanced_Pseudobulk_20260806")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Rdata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "Logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "Logs", "SC28_RECALC_02_log.txt")
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

for (pkg in c("data.table", "Matrix", "ggplot2")) install_cran_if_missing(pkg)

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  stop_clean("Package 'SeuratObject' is required in the original single-cell R environment.")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in c("AnnotationDbi", "org.Hs.eg.db", "edgeR")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg("Installing Bioconductor package: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(SeuratObject)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(edgeR)
})

msg("============================================================")
msg("SC28_RECALC_02: proliferation-balanced paired pseudobulk")
msg("R version: ", R.version.string)
msg("edgeR version: ", as.character(utils::packageVersion("edgeR")))
msg("Output: ", out_dir)
msg("============================================================")

if (!file.exists(seurat_path)) stop_clean("Missing Seurat object: ", seurat_path)
if (!file.exists(score_rds)) stop_clean("Missing completed Phase-1 score RDS: ", score_rds)

# -----------------------------------------------------------------------------
# 1. Audited proliferation/cell-cycle definition.
# -----------------------------------------------------------------------------
# The explicit list and regex below are inherited from the historical SC28B
# code so that the biological definition is traceable. Unlike historical SC28B,
# the present analysis also balances proliferation BEFORE pseudobulk testing.

proliferation_genes <- unique(c(
  "MKI67", "TOP2A", "BIRC5", "CDK1", "UBE2C", "PCNA",
  "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
  "HMGB2", "TYMS", "CENPF", "AURKA", "AURKB",
  "CCNA2", "CCNB1", "CCNB2", "CCNE1", "CCNE2",
  "CDC20", "CDC25A", "CDC25B", "CDC25C",
  "CDCA2", "CDCA3", "CDCA5", "CDCA7", "CDCA8",
  "NUSAP1", "TPX2", "ASPM",
  "KIF2C", "KIF4A", "KIF11", "KIF14", "KIF15", "KIF20A", "KIF23",
  "RRM1", "RRM2", "MAD2L1", "BUB1", "BUB1B", "PLK1",
  "UBE2T", "UBE2S", "MELK", "DLGAP5", "MND1", "NCAPG", "NCAPH",
  "TUBA1B", "TUBB", "TUBB4B", "STMN1"
))

prolif_regex <- paste(c(
  "^MKI67$", "^TOP2A$", "^BIRC5$", "^CDK1$", "^UBE2C$", "^PCNA$",
  "^MCM[0-9]+", "^CENP", "^AURK", "^KIF", "^CCN[ABDE]", "^CDC",
  "^CDCA", "^HIST", "^HMGB2$", "^TYMS$", "^RRM", "^BUB", "^PLK1$",
  "^NUSAP1$", "^TPX2$", "^ASPM$", "^DLGAP5$", "^NCAP", "^TUBA",
  "^TUBB", "^STMN1$"
), collapse = "|")

is_proliferation_gene <- function(gene) {
  gene <- as.character(gene)
  gene %in% proliferation_genes | grepl(prolif_regex, gene, ignore.case = FALSE)
}

# Fixed analysis parameters.
n_prolif_strata <- 20L
fdr_cutoff <- 0.05
logfc_cutoff <- 0.10
min_matched_cells_per_state <- 20L

# -----------------------------------------------------------------------------
# 2. Load and validate the completed 01B LactateScore result.
# -----------------------------------------------------------------------------

scores <- as.data.table(readRDS(score_rds))
required_score_cols <- c(
  "CellID", "Sample_sc", "LactateScore", "LactateScore_sampleZ",
  "LactateState", "PRDX4_logNorm"
)
missing_score_cols <- setdiff(required_score_cols, names(scores))
if (length(missing_score_cols)) {
  stop_clean("Phase-1 score table is missing columns: ", paste(missing_score_cols, collapse = ", "))
}

scores[, CellID := as.character(CellID)]
scores[, Sample_sc := as.character(Sample_sc)]
scores[, LactateState := as.character(LactateState)]

if (anyDuplicated(scores$CellID)) stop_clean("Duplicate CellID values found in Phase-1 scores.")
if (nrow(scores) != 111132L) stop_clean("Expected 111132 core tumor cells from completed 01B; found ", nrow(scores))

tumor_samples <- sort(unique(scores$Sample_sc))
if (length(tumor_samples) != 16L) stop_clean("Expected 16 tumor samples; found ", length(tumor_samples))
if (any(!grepl("_T[0-9]+$", tumor_samples))) stop_clean("A non-T sample is present in Phase-1 scores.")

state_counts <- scores[, .N, by = LactateState]
if (state_counts[LactateState == "LactateHigh", N] != 33333L ||
    state_counts[LactateState == "LactateLow", N] != 33333L) {
  stop_clean("Phase-1 exact High/Low counts do not match the audited 33333/33333 result.")
}

msg("Validated Phase-1 input: 111132 cells; 16 T samples; High/Low = 33333/33333.")

# -----------------------------------------------------------------------------
# 3. Load Seurat object and build a consistent feature -> HGNC-symbol map.
# -----------------------------------------------------------------------------

msg("Loading original GSE137804 Seurat object (~4.5 GB)...")
obj <- readRDS(seurat_path)
if (!"RNA" %in% names(obj@assays)) stop_clean("RNA assay is missing from the Seurat object.")
layer_names <- SeuratObject::Layers(obj[["RNA"]])

find_counts_layer <- function(sample_name) {
  exact <- paste0("counts.", sample_name)
  if (exact %in% layer_names) return(exact)
  hits <- layer_names[grepl("^counts", layer_names) & endsWith(layer_names, sample_name)]
  if (length(hits) == 1L) return(hits)
  NA_character_
}

official_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
official_symbols_set <- unique(as.character(official_symbols))

safe_map_ids <- function(keys, keytype) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) return(character(0))
  out <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = keys,
    column = "SYMBOL",
    keytype = keytype,
    multiVals = "first"
  )
  out_chr <- as.character(out)
  names(out_chr) <- names(out)
  out_chr
}

map_all_features <- function(features) {
  features <- as.character(features)
  clean <- sub("\\..*$", "", features)
  symbol <- rep(NA_character_, length(features))

  direct <- features %in% official_symbols_set
  symbol[direct] <- features[direct]

  ens_idx <- which(is.na(symbol) & grepl("^ENSG[0-9]+$", clean))
  if (length(ens_idx)) {
    ens_keys <- clean[ens_idx]
    ens_map <- safe_map_ids(ens_keys, "ENSEMBL")
    symbol[ens_idx] <- unname(ens_map[ens_keys])
  }

  entrez_idx <- which(is.na(symbol) & grepl("^[0-9]+$", features))
  if (length(entrez_idx)) {
    entrez_keys <- features[entrez_idx]
    entrez_map <- safe_map_ids(entrez_keys, "ENTREZID")
    symbol[entrez_idx] <- unname(entrez_map[entrez_keys])
  }
  symbol
}

feature_maps <- vector("list", length(tumor_samples))
feature_ids <- vector("list", length(tumor_samples))
names(feature_maps) <- names(feature_ids) <- tumor_samples
map_qc <- vector("list", length(tumor_samples))

for (i in seq_along(tumor_samples)) {
  smp <- tumor_samples[i]
  msg(sprintf("Feature mapping %02d/%02d: %s", i, length(tumor_samples), smp))
  ly <- find_counts_layer(smp)
  if (is.na(ly)) stop_clean("No unique counts layer for ", smp)
  mat <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  feats <- rownames(mat)
  fmap <- map_all_features(feats)
  feature_ids[[smp]] <- feats
  feature_maps[[smp]] <- fmap
  map_qc[[smp]] <- data.table(
    Sample_sc = smp,
    Layer = ly,
    FeatureN = length(feats),
    MappedSymbolN = sum(!is.na(fmap)),
    UniqueMappedSymbolN = uniqueN(fmap[!is.na(fmap)]),
    ProliferationGeneN = uniqueN(intersect(fmap[!is.na(fmap)], proliferation_genes)),
    PRDX4Mapped = "PRDX4" %in% fmap
  )
  rm(mat)
  gc(verbose = FALSE)
}

map_qc_dt <- rbindlist(map_qc)
fwrite(map_qc_dt, file.path(out_dir, "Tables", "01_feature_mapping_QC.csv"))
if (any(!map_qc_dt$PRDX4Mapped)) stop_clean("PRDX4 mapping failed in at least one tumor sample.")

available_symbol_lists <- lapply(feature_maps, function(x) unique(x[!is.na(x)]))
common_expression_genes <- sort(Reduce(intersect, available_symbol_lists))
common_proliferation_genes <- sort(intersect(proliferation_genes, common_expression_genes))

if (length(common_expression_genes) < 5000L) {
  stop_clean("Too few common mapped genes across all 16 tumors: ", length(common_expression_genes))
}
if (length(common_proliferation_genes) < 20L) {
  stop_clean("Too few common proliferation genes across all 16 tumors: ", length(common_proliferation_genes))
}

fwrite(data.table(GeneSymbol = common_expression_genes), file.path(out_dir, "Tables", "02_common_expression_genes_all_16_tumors.csv"))
fwrite(data.table(GeneSymbol = common_proliferation_genes), file.path(out_dir, "Tables", "03_common_proliferation_genes_all_16_tumors.csv"))
msg("Common mapped genes: ", length(common_expression_genes))
msg("Common explicit proliferation genes: ", length(common_proliferation_genes), " / ", length(proliferation_genes))

# -----------------------------------------------------------------------------
# 4. Calculate a sample-wise proliferation program score from raw RNA counts.
# -----------------------------------------------------------------------------

collapse_selected_counts <- function(mat, feature_symbols, target_genes) {
  hit <- which(!is.na(feature_symbols) & feature_symbols %in% target_genes)
  if (!length(hit)) stop_clean("No requested target genes were mapped in a layer.")
  mp <- data.table(RowIndex = hit, GeneSymbol = feature_symbols[hit])
  mp <- unique(mp, by = c("RowIndex", "GeneSymbol"))
  gene_i <- match(mp$GeneSymbol, target_genes)
  collapse_matrix <- Matrix::sparseMatrix(
    i = gene_i,
    j = mp$RowIndex,
    x = 1,
    dims = c(length(target_genes), nrow(mat)),
    dimnames = list(target_genes, rownames(mat))
  )
  collapse_matrix %*% mat
}

scores[, ProliferationScore := NA_real_]
scores[, ProliferationScore_sampleZ := NA_real_]
prolif_score_qc <- vector("list", length(tumor_samples))

for (i in seq_along(tumor_samples)) {
  smp <- tumor_samples[i]
  msg(sprintf("Proliferation scoring %02d/%02d: %s", i, length(tumor_samples), smp))
  ly <- find_counts_layer(smp)
  mat_all <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  if (!identical(rownames(mat_all), feature_ids[[smp]])) stop_clean("Feature order changed unexpectedly for ", smp)

  cells <- scores[Sample_sc == smp, CellID]
  idx <- match(cells, colnames(mat_all))
  if (anyNA(idx)) stop_clean("Missing core cells in counts layer for ", smp, ": ", sum(is.na(idx)))
  mat <- mat_all[, idx, drop = FALSE]
  libsize <- Matrix::colSums(mat)
  if (any(libsize <= 0)) stop_clean("Zero-library cells in ", smp)

  pcounts <- collapse_selected_counts(mat, feature_maps[[smp]], common_proliferation_genes)
  pnorm <- pcounts %*% Matrix::Diagonal(x = 10000 / libsize)
  pnorm@x <- log1p(pnorm@x)
  pdense <- as.matrix(pnorm)
  pmean <- rowMeans(pdense)
  psd <- apply(pdense, 1, stats::sd)
  good <- is.finite(psd) & psd > 0
  if (sum(good) < 20L) stop_clean("Too few variable proliferation genes in ", smp)
  pz <- sweep(pdense[good, , drop = FALSE], 1, pmean[good], "-")
  pz <- sweep(pz, 1, psd[good], "/")
  pscore <- colMeans(pz)
  pscore_z <- as.numeric(scale(pscore))

  row_idx <- which(scores$Sample_sc == smp)
  ord <- match(scores$CellID[row_idx], cells)
  scores$ProliferationScore[row_idx] <- pscore[ord]
  scores$ProliferationScore_sampleZ[row_idx] <- pscore_z[ord]

  prolif_score_qc[[smp]] <- data.table(
    Sample_sc = smp,
    CellN = length(cells),
    ProliferationGenesUsedN = sum(good),
    ScoreSD = stats::sd(pscore),
    UniqueScoreN = uniqueN(signif(pscore, 12))
  )
  rm(mat_all, mat, pcounts, pnorm, pdense, pz)
  gc(verbose = FALSE)
}

prolif_score_qc_dt <- rbindlist(prolif_score_qc)
fwrite(prolif_score_qc_dt, file.path(out_dir, "Tables", "04_proliferation_score_QC.csv"))
if (any(!is.finite(scores$ProliferationScore))) stop_clean("Non-finite proliferation scores were generated.")

saveRDS(scores, file.path(out_dir, "Rdata", "05_phase1_scores_with_proliferation_score.rds"), compress = TRUE)

# -----------------------------------------------------------------------------
# 5. Within-sample proliferation balancing.
# -----------------------------------------------------------------------------
# Coarsened exact matching is performed within each tumor. Cells are divided
# into 20 quantile strata of the proliferation score; an equal number of
# LactateHigh and LactateLow cells is retained from every stratum. Selection
# within an over-represented state is deterministic and evenly spaced across the
# sorted proliferation-score distribution, so the result is reproducible.

smd <- function(x_high, x_low) {
  x_high <- as.numeric(x_high)
  x_low <- as.numeric(x_low)
  pooled <- sqrt((stats::var(x_high) + stats::var(x_low)) / 2)
  if (!is.finite(pooled) || pooled == 0) return(0)
  (mean(x_high) - mean(x_low)) / pooled
}

even_pick <- function(dt, n_pick) {
  if (nrow(dt) <= n_pick) return(dt)
  setorder(dt, ProliferationScore, CellID)
  pos <- unique(round(seq(1, nrow(dt), length.out = n_pick)))
  if (length(pos) < n_pick) {
    remaining <- setdiff(seq_len(nrow(dt)), pos)
    pos <- c(pos, head(remaining, n_pick - length(pos)))
  }
  dt[sort(pos[seq_len(n_pick)])]
}

match_one_sample <- function(smp) {
  d <- copy(scores[Sample_sc == smp & LactateState %in% c("LactateLow", "LactateHigh")])
  if (!nrow(d)) stop_clean("No High/Low cells for ", smp)
  if (d[LactateState == "LactateHigh", .N] == 0L) stop_clean("No High cells for ", smp)
  if (d[LactateState == "LactateLow", .N] == 0L) stop_clean("No Low cells for ", smp)

  q <- unique(as.numeric(stats::quantile(
    d$ProliferationScore,
    probs = seq(0, 1, length.out = n_prolif_strata + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )))
  if (length(q) >= 3L) {
    d[, MatchStratum := cut(ProliferationScore, breaks = q, include.lowest = TRUE, labels = FALSE)]
  } else {
    # Very unlikely fallback for extensive ties.
    d[, MatchStratum := pmin(
      n_prolif_strata,
      floor((frank(ProliferationScore, ties.method = "average") - 1) / .N * n_prolif_strata) + 1L
    )]
  }

  matched_parts <- list()
  counter <- 0L
  for (st in sort(unique(d$MatchStratum))) {
    ds <- d[MatchStratum == st]
    h <- copy(ds[LactateState == "LactateHigh"])
    l <- copy(ds[LactateState == "LactateLow"])
    n_take <- min(nrow(h), nrow(l))
    if (n_take <= 0L) next
    h <- even_pick(h, n_take)
    l <- even_pick(l, n_take)
    counter <- counter + 1L
    matched_parts[[counter]] <- rbindlist(list(h, l), use.names = TRUE)
  }
  if (!length(matched_parts)) stop_clean("No matched cells retained for ", smp)
  m <- rbindlist(matched_parts, use.names = TRUE)
  m[, MatchedForProliferation := TRUE]
  m
}

matched_list <- lapply(tumor_samples, match_one_sample)
matched <- rbindlist(matched_list, use.names = TRUE, fill = TRUE)

matching_qc <- rbindlist(lapply(tumor_samples, function(smp) {
  b <- scores[Sample_sc == smp & LactateState %in% c("LactateLow", "LactateHigh")]
  a <- matched[Sample_sc == smp]
  bh <- b[LactateState == "LactateHigh", ProliferationScore]
  bl <- b[LactateState == "LactateLow", ProliferationScore]
  ah <- a[LactateState == "LactateHigh", ProliferationScore]
  al <- a[LactateState == "LactateLow", ProliferationScore]
  data.table(
    Sample_sc = smp,
    HighN_Before = length(bh),
    LowN_Before = length(bl),
    HighN_After = length(ah),
    LowN_After = length(al),
    RetentionFraction = length(ah) / length(bh),
    ProlifMeanHigh_Before = mean(bh),
    ProlifMeanLow_Before = mean(bl),
    SMD_Before = smd(bh, bl),
    ProlifMeanHigh_After = mean(ah),
    ProlifMeanLow_After = mean(al),
    SMD_After = smd(ah, al)
  )
}))

if (any(matching_qc$HighN_After < min_matched_cells_per_state |
        matching_qc$LowN_After < min_matched_cells_per_state)) {
  stop_clean("At least one tumor has fewer than ", min_matched_cells_per_state, " matched cells per state.")
}

fwrite(matching_qc, file.path(out_dir, "Tables", "06_proliferation_matching_QC_by_sample.csv"))
fwrite(
  matched[, .(CellID, Sample_sc, LactateState, LactateScore, LactateScore_sampleZ,
              ProliferationScore, ProliferationScore_sampleZ, PRDX4_logNorm, MatchStratum)],
  file.path(out_dir, "Tables", "07_matched_cell_table.csv.gz"),
  compress = "gzip"
)
saveRDS(matched, file.path(out_dir, "Rdata", "07_matched_cell_table.rds"), compress = TRUE)

msg("Matched cells per state total: High = ", matched[LactateState == "LactateHigh", .N],
    "; Low = ", matched[LactateState == "LactateLow", .N])
msg("Maximum absolute proliferation SMD: before = ", signif(max(abs(matching_qc$SMD_Before)), 4),
    "; after = ", signif(max(abs(matching_qc$SMD_After)), 4))

# -----------------------------------------------------------------------------
# 6. Aggregate proliferation-balanced cells into 16 paired pseudobulk profiles.
# -----------------------------------------------------------------------------

pb_columns <- list()
pb_meta <- list()
pb_counter <- 0L

for (i in seq_along(tumor_samples)) {
  smp <- tumor_samples[i]
  msg(sprintf("Pseudobulk aggregation %02d/%02d: %s", i, length(tumor_samples), smp))
  ly <- find_counts_layer(smp)
  mat_all <- SeuratObject::LayerData(obj, assay = "RNA", layer = ly)
  if (!identical(rownames(mat_all), feature_ids[[smp]])) stop_clean("Feature order changed unexpectedly for ", smp)
  fmap <- feature_maps[[smp]]

  for (state in c("LactateLow", "LactateHigh")) {
    cells <- matched[Sample_sc == smp & LactateState == state, CellID]
    idx <- match(cells, colnames(mat_all))
    if (anyNA(idx)) stop_clean("Matched cells missing from layer for ", smp, " / ", state)
    feature_sum <- Matrix::rowSums(mat_all[, idx, drop = FALSE])
    tmp <- data.table(GeneSymbol = fmap, Count = as.numeric(feature_sum))
    tmp <- tmp[!is.na(GeneSymbol) & GeneSymbol %in% common_expression_genes,
               .(Count = sum(Count)), by = GeneSymbol]
    v <- tmp$Count[match(common_expression_genes, tmp$GeneSymbol)]
    if (anyNA(v)) stop_clean("A common expression gene disappeared during pseudobulk aggregation for ", smp)

    pb_counter <- pb_counter + 1L
    pb_id <- paste(smp, state, sep = "__")
    pb_columns[[pb_counter]] <- v
    pb_meta[[pb_counter]] <- data.table(
      PseudobulkID = pb_id,
      Sample_sc = smp,
      LactateState = state,
      MatchedCellN = length(cells),
      MeanProliferationScore = mean(matched[Sample_sc == smp & LactateState == state, ProliferationScore])
    )
  }
  rm(mat_all)
  gc(verbose = FALSE)
}

pb_counts <- do.call(cbind, pb_columns)
rownames(pb_counts) <- common_expression_genes
pb_meta_dt <- rbindlist(pb_meta)
colnames(pb_counts) <- pb_meta_dt$PseudobulkID

fwrite(pb_meta_dt, file.path(out_dir, "Tables", "08_pseudobulk_metadata.csv"))
saveRDS(pb_counts, file.path(out_dir, "Rdata", "09_pseudobulk_raw_counts_matrix.rds"), compress = TRUE)

# -----------------------------------------------------------------------------
# 7. Paired edgeR quasi-likelihood differential-expression analysis.
# -----------------------------------------------------------------------------

pb_meta_dt[, SampleFactor := factor(Sample_sc)]
pb_meta_dt[, StateFactor := factor(LactateState, levels = c("LactateLow", "LactateHigh"))]
design <- stats::model.matrix(~ SampleFactor + StateFactor, data = pb_meta_dt)
rownames(design) <- pb_meta_dt$PseudobulkID

y <- edgeR::DGEList(counts = pb_counts)
keep <- edgeR::filterByExpr(y, design = design)
if (!"PRDX4" %in% rownames(y)) stop_clean("PRDX4 is absent from the common pseudobulk matrix.")
if (!keep["PRDX4"]) stop_clean("PRDX4 did not pass edgeR filterByExpr; stop rather than force the gene into testing.")

y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")
y <- edgeR::estimateDisp(y, design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design, robust = TRUE)
coef_name <- "StateFactorLactateHigh"
coef_idx <- match(coef_name, colnames(design))
if (is.na(coef_idx)) stop_clean("Could not identify the LactateHigh coefficient in the paired design.")
qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
tab <- as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table, keep.rownames = "Gene")
setnames(tab, "FDR", "FDR_edgeR")
setnames(tab, "PValue", "PValue_edgeR")
tab[, Is_CellCycle_Proliferation := is_proliferation_gene(Gene)]
tab[, Direction := fcase(
  FDR_edgeR < fdr_cutoff & logFC >= logfc_cutoff, "High_up",
  FDR_edgeR < fdr_cutoff & logFC <= -logfc_cutoff, "Low_up",
  default = "NS"
)]
tab[, AbsLogFC_for_order := abs(logFC)]
setorder(tab, PValue_edgeR, -AbsLogFC_for_order)
tab[, AbsLogFC_for_order := NULL]
fwrite(tab, file.path(out_dir, "Tables", "10_paired_edgeR_DEG_all.csv"))

candidate_all_high <- tab[!Is_CellCycle_Proliferation & logFC > 0]
setorder(candidate_all_high, PValue_edgeR, -logFC)
candidate_all_high[, Rank_NonProliferative_HighUp := seq_len(.N)]

candidate_sig_high <- candidate_all_high[FDR_edgeR < fdr_cutoff & logFC >= logfc_cutoff]
candidate_sig_high[, Rank_Significant_NonProliferative_HighUp := seq_len(.N)]

fwrite(candidate_all_high, file.path(out_dir, "Tables", "11_nonproliferative_HighUp_all_ranked.csv"))
fwrite(candidate_sig_high, file.path(out_dir, "Tables", "12_nonproliferative_HighUp_significant_ranked.csv"))

prdx4 <- tab[Gene == "PRDX4"]
if (nrow(prdx4) != 1L) stop_clean("Expected exactly one PRDX4 row after edgeR; found ", nrow(prdx4))
prdx4_rank_all <- candidate_all_high[Gene == "PRDX4", Rank_NonProliferative_HighUp]
if (!length(prdx4_rank_all)) prdx4_rank_all <- NA_integer_
prdx4_rank_sig <- candidate_sig_high[Gene == "PRDX4", Rank_Significant_NonProliferative_HighUp]
if (!length(prdx4_rank_sig)) prdx4_rank_sig <- NA_integer_

# Normalized pseudobulk expression for an interpretable paired PRDX4 plot/test.
logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
prdx4_pb <- data.table(
  PseudobulkID = colnames(logcpm),
  PRDX4_logCPM = as.numeric(logcpm["PRDX4", ])
)
prdx4_pb <- merge(prdx4_pb, pb_meta_dt[, .(PseudobulkID, Sample_sc, LactateState, MatchedCellN, MeanProliferationScore)], by = "PseudobulkID")
setorder(prdx4_pb, Sample_sc, LactateState)
fwrite(prdx4_pb, file.path(out_dir, "Tables", "13_PRDX4_matched_pseudobulk_expression.csv"))

prdx4_wide <- dcast(prdx4_pb, Sample_sc ~ LactateState, value.var = "PRDX4_logCPM")
paired_wilcox <- stats::wilcox.test(
  prdx4_wide$LactateHigh,
  prdx4_wide$LactateLow,
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE
)
prdx4_positive_samples <- sum(prdx4_wide$LactateHigh > prdx4_wide$LactateLow)

prdx4_summary <- data.table(
  Item = c(
    "PRDX4_logFC_High_vs_Low_edgeR",
    "PRDX4_P_edgeR",
    "PRDX4_FDR_edgeR",
    "PRDX4_is_proliferation_filtered",
    "PRDX4_rank_all_nonproliferative_HighUp",
    "PRDX4_in_significant_nonproliferative_HighUp",
    "PRDX4_rank_significant_nonproliferative_HighUp",
    "PRDX4_paired_pseudobulk_Wilcoxon_P",
    "PRDX4_High_gt_Low_samples_after_matching",
    "Nonproliferative_HighUp_geneN",
    "Significant_nonproliferative_HighUp_geneN"
  ),
  Value = c(
    as.character(prdx4$logFC),
    as.character(prdx4$PValue_edgeR),
    as.character(prdx4$FDR_edgeR),
    as.character(prdx4$Is_CellCycle_Proliferation),
    as.character(prdx4_rank_all),
    as.character("PRDX4" %in% candidate_sig_high$Gene),
    as.character(prdx4_rank_sig),
    as.character(paired_wilcox$p.value),
    paste0(prdx4_positive_samples, "/", nrow(prdx4_wide)),
    as.character(nrow(candidate_all_high)),
    as.character(nrow(candidate_sig_high))
  )
)
fwrite(prdx4_summary, file.path(out_dir, "Tables", "14_PRDX4_rank_and_inference_summary.csv"))

# -----------------------------------------------------------------------------
# 8. Audit figures with publication-facing titles (no internal SC labels).
# -----------------------------------------------------------------------------

balance_long <- rbindlist(list(
  matching_qc[, .(Sample_sc, Stage = "Before matching", SMD = SMD_Before)],
  matching_qc[, .(Sample_sc, Stage = "After matching", SMD = SMD_After)]
))
balance_long[, Stage := factor(Stage, levels = c("Before matching", "After matching"))]
p_balance <- ggplot(balance_long, aes(x = Stage, y = SMD, group = Sample_sc)) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_line(alpha = 0.45, color = "grey55") +
  geom_point(aes(color = Stage), size = 2) +
  scale_color_manual(values = c("Before matching" = "#D95F02", "After matching" = "#1B9E77")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(
    title = "Proliferation-score balance before and after within-sample matching",
    x = NULL, y = "Standardized mean difference"
  )
ggsave(file.path(out_dir, "Figures", "01_proliferation_balance.png"), p_balance, width = 6.8, height = 5.4, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Figures", "01_proliferation_balance.pdf"), p_balance, width = 6.8, height = 5.4, bg = "white")

p_prdx4 <- ggplot(prdx4_pb, aes(x = LactateState, y = PRDX4_logCPM, group = Sample_sc)) +
  geom_line(color = "grey65", linewidth = 0.45, alpha = 0.8) +
  geom_point(aes(color = LactateState), size = 2.1) +
  scale_color_manual(values = c("LactateLow" = "#2166AC", "LactateHigh" = "#B2182B")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(
    title = "PRDX4 expression after proliferation-balanced matching",
    subtitle = "Paired pseudobulk profiles from 16 primary neuroblastomas",
    x = NULL, y = "Normalized pseudobulk logCPM"
  )
ggsave(file.path(out_dir, "Figures", "02_PRDX4_matched_pseudobulk.png"), p_prdx4, width = 6.8, height = 5.4, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Figures", "02_PRDX4_matched_pseudobulk.pdf"), p_prdx4, width = 6.8, height = 5.4, bg = "white")

volcano <- copy(tab)
volcano[, neglog10FDR := -log10(pmax(FDR_edgeR, 1e-300))]
volcano[, PlotGroup := fcase(
  Gene == "PRDX4", "PRDX4",
  Is_CellCycle_Proliferation & Direction == "High_up", "Proliferation/cell-cycle High-up",
  Direction == "High_up", "High-up",
  Direction == "Low_up", "Low-up",
  default = "Not significant"
)]
p_volcano <- ggplot(volcano, aes(x = logFC, y = neglog10FDR, color = PlotGroup)) +
  geom_point(size = 0.9, alpha = 0.75) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", linewidth = 0.3) +
  geom_text(data = volcano[Gene == "PRDX4"], aes(label = Gene), color = "#6A00A8", vjust = -0.8, size = 4, fontface = "bold") +
  scale_color_manual(values = c(
    "Not significant" = "grey78", "Low-up" = "#2166AC", "High-up" = "#B2182B",
    "Proliferation/cell-cycle High-up" = "#F4A582", "PRDX4" = "#6A00A8"
  )) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.title = element_blank()) +
  labs(
    title = "Differential expression after proliferation-balanced matching",
    subtitle = "Paired edgeR quasi-likelihood analysis across 16 primary tumors",
    x = "log2 fold change: LactateHigh vs LactateLow", y = "-log10(FDR)"
  )
ggsave(file.path(out_dir, "Figures", "03_matched_pseudobulk_volcano.png"), p_volcano, width = 8.2, height = 6.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "Figures", "03_matched_pseudobulk_volcano.pdf"), p_volcano, width = 8.2, height = 6.2, bg = "white")

# -----------------------------------------------------------------------------
# 9. Final audit summary.
# -----------------------------------------------------------------------------

final_summary <- data.table(
  Item = c(
    "Input_core_tumor_cells",
    "Input_primary_tumors",
    "Input_LactateHigh_cells",
    "Input_LactateLow_cells",
    "Common_expression_genes_all_16",
    "Common_explicit_proliferation_genes",
    "Proliferation_matching_strata",
    "Matched_High_cells",
    "Matched_Low_cells",
    "Minimum_retention_fraction",
    "Maximum_abs_SMD_before",
    "Maximum_abs_SMD_after",
    "edgeR_tested_geneN",
    "Significant_HighUp_geneN",
    "Significant_LowUp_geneN",
    "PRDX4_logFC",
    "PRDX4_P_edgeR",
    "PRDX4_FDR_edgeR",
    "PRDX4_rank_all_nonproliferative_HighUp",
    "PRDX4_rank_significant_nonproliferative_HighUp",
    "PRDX4_High_gt_Low_samples_after_matching"
  ),
  Value = c(
    nrow(scores),
    length(tumor_samples),
    scores[LactateState == "LactateHigh", .N],
    scores[LactateState == "LactateLow", .N],
    length(common_expression_genes),
    length(common_proliferation_genes),
    n_prolif_strata,
    matched[LactateState == "LactateHigh", .N],
    matched[LactateState == "LactateLow", .N],
    min(matching_qc$RetentionFraction),
    max(abs(matching_qc$SMD_Before)),
    max(abs(matching_qc$SMD_After)),
    nrow(tab),
    tab[Direction == "High_up", .N],
    tab[Direction == "Low_up", .N],
    prdx4$logFC,
    prdx4$PValue_edgeR,
    prdx4$FDR_edgeR,
    prdx4_rank_all,
    prdx4_rank_sig,
    paste0(prdx4_positive_samples, "/", nrow(prdx4_wide))
  )
)
fwrite(final_summary, file.path(out_dir, "Tables", "15_RECALC_02_final_summary.csv"))

writeLines(c(
  "SC28_RECALC_02 COMPLETED SUCCESSFULLY",
  paste0("Finished: ", Sys.time()),
  paste0("Output: ", out_dir),
  paste0("Matched High/Low cells: ", matched[LactateState == "LactateHigh", .N], "/", matched[LactateState == "LactateLow", .N]),
  paste0("Maximum |SMD| after proliferation matching: ", signif(max(abs(matching_qc$SMD_After)), 6)),
  paste0("edgeR tested genes: ", nrow(tab)),
  paste0("PRDX4 logFC: ", signif(prdx4$logFC, 6)),
  paste0("PRDX4 P/FDR: ", signif(prdx4$PValue_edgeR, 6), " / ", signif(prdx4$FDR_edgeR, 6)),
  paste0("PRDX4 rank among all non-proliferative High-up genes: ", prdx4_rank_all),
  paste0("PRDX4 rank among significant non-proliferative High-up genes: ", prdx4_rank_sig),
  paste0("PRDX4 High > Low samples after matching: ", prdx4_positive_samples, "/", nrow(prdx4_wide))
), file.path(out_dir, "RUN_STATUS.txt"))

capture.output(utils::sessionInfo(), file = file.path(out_dir, "Logs", "sessionInfo.txt"))

msg("============================================================")
msg("SC28_RECALC_02 COMPLETED SUCCESSFULLY")
msg("Matched High / Low cells = ", matched[LactateState == "LactateHigh", .N], " / ", matched[LactateState == "LactateLow", .N])
msg("Maximum |SMD| after matching = ", signif(max(abs(matching_qc$SMD_After)), 6))
msg("PRDX4 logFC = ", signif(prdx4$logFC, 6), "; P = ", signif(prdx4$PValue_edgeR, 6), "; FDR = ", signif(prdx4$FDR_edgeR, 6))
msg("PRDX4 rank among all non-proliferative High-up genes = ", prdx4_rank_all)
msg("PRDX4 rank among significant non-proliferative High-up genes = ", prdx4_rank_sig)
msg("PRDX4 High > Low samples after matching = ", prdx4_positive_samples, "/", nrow(prdx4_wide))
msg("Please zip and return the complete output folder.")
msg("============================================================")
