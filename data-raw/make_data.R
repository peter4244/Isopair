#!/usr/bin/env Rscript
#
# Build all package data files for Isopair
#
# Run from the package root: Rscript data-raw/make_data.R
#
# Requires: rtracklayer, Isopair (devtools::load_all("."))
#
# Inputs:
#   - GENCODE v49 GTF (path below, change if needed)
#
# Outputs:
#   - inst/extdata/test_real_genes.gtf   (8 real genes)
#   - inst/extdata/example_small.gtf     (5-gene subset)
#   - inst/extdata/example_expression.csv (mock CPM matrix)
#   - inst/extdata/example_de.csv        (mock DE results)
#   - inst/extdata/example_du.csv        (mock DU results)
#   - data/example_structures.rda
#   - data/example_pairs.rda
#   - data/example_profiles.rda

library(rtracklayer)

# ============================================================================
# Configuration
# ============================================================================

gencode_gtf <- "~/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.gtf"

target_genes <- c("MTCL1", "SRSF3", "SRSF1", "FAM13A", "AGER", "BCL2L1",
                  "TP63", "TTN")

# TTN isoforms to keep (avoid 365-exon variants)
ttn_keep_transcripts <- c("ENST00000589042.6",   # TTN-201 (canonical short)
                          "ENST00000342992.12",   # TTN-202 (moderate)
                          "ENST00000460472.6")     # TTN-208 (short)

# Example subset: 5 genes with ~5 isoforms each
example_genes <- c("SRSF3", "SRSF1", "BCL2L1", "TP63", "FAM13A")

# Max isoforms per gene in example_small (canonical + 4 others)
max_isoforms_example <- 5

# ============================================================================
# Step 1: Extract real gene test data
# ============================================================================

cat("Loading GENCODE v49 GTF...\n")
gtf <- import(gencode_gtf)

cat("Filtering to target genes...\n")
gene_gtf <- gtf[gtf$gene_name %in% target_genes]

# For TTN: further filter to named transcripts
ttn_features <- gene_gtf[gene_gtf$gene_name == "TTN"]
non_ttn <- gene_gtf[gene_gtf$gene_name != "TTN"]

# Keep gene-level and selected transcript features for TTN
ttn_gene_level <- ttn_features[ttn_features$type == "gene"]
ttn_tx_features <- ttn_features[
  !is.na(ttn_features$transcript_id) &
  ttn_features$transcript_id %in% ttn_keep_transcripts]
ttn_filtered <- c(ttn_gene_level, ttn_tx_features)

gene_gtf_final <- c(non_ttn, ttn_filtered)
gene_gtf_final <- sort(gene_gtf_final)

cat(sprintf("  %d features for %d genes\n",
            length(gene_gtf_final),
            length(unique(gene_gtf_final$gene_name))))

outpath1 <- "inst/extdata/test_real_genes.gtf"
export(gene_gtf_final, outpath1)
cat(sprintf("  Wrote: %s\n\n", outpath1))

# ============================================================================
# Step 2: Create example_small.gtf
# ============================================================================

cat("Creating example_small.gtf...\n")

example_gtf <- gene_gtf_final[gene_gtf_final$gene_name %in% example_genes]

# For each gene, select canonical + diverse isoforms
example_exons <- as.data.frame(example_gtf[example_gtf$type == "exon"])
tx_counts <- table(example_exons$transcript_id)

# Per-gene transcript selection
selected_tx <- character(0)
for (gn in example_genes) {
  gene_exons <- example_exons[example_exons$gene_name == gn, ]
  gene_txs <- unique(gene_exons$transcript_id)

  if (length(gene_txs) <= max_isoforms_example) {
    selected_tx <- c(selected_tx, gene_txs)
    next
  }

  # Sort by exon count (diversity proxy) and pick spread
  exon_counts <- sort(table(gene_exons$transcript_id), decreasing = TRUE)
  tx_names <- names(exon_counts)

  # Strategy: pick 1st (most exons), last (fewest exons), and 3 evenly spaced
  n <- length(tx_names)
  indices <- unique(round(seq(1, n, length.out = max_isoforms_example)))
  selected_tx <- c(selected_tx, tx_names[indices])
}

# Keep gene-level features + selected transcripts
example_gene_level <- example_gtf[example_gtf$type == "gene"]
example_tx_features <- example_gtf[
  !is.na(example_gtf$transcript_id) &
  example_gtf$transcript_id %in% selected_tx]
example_small <- c(example_gene_level, example_tx_features)
example_small <- sort(example_small)

cat(sprintf("  %d features, %d transcripts across %d genes\n",
            length(example_small),
            length(unique(example_small$transcript_id[
              !is.na(example_small$transcript_id)])),
            length(unique(example_small$gene_name))))

outpath2 <- "inst/extdata/example_small.gtf"
export(example_small, outpath2)
cat(sprintf("  Wrote: %s\n\n", outpath2))

# ============================================================================
# Step 3: Create example tabular data
# ============================================================================

cat("Creating example tabular data...\n")

# Get isoform IDs from example_small
example_exon_df <- as.data.frame(example_small[example_small$type == "exon"])
isoform_ids <- unique(example_exon_df$transcript_id)
gene_ids <- example_exon_df$gene_id[
  match(isoform_ids, example_exon_df$transcript_id)]

# Gene map
gene_map <- data.frame(
  isoform_id = isoform_ids,
  gene_id = gene_ids,
  stringsAsFactors = FALSE
)

set.seed(42)
n_iso <- length(isoform_ids)

# -- Expression matrix (mock CPM) --
sample_names <- paste0("S", 1:6)
expr_mat <- matrix(0, nrow = n_iso, ncol = 6,
                   dimnames = list(isoform_ids, sample_names))

for (gn in unique(gene_ids)) {
  iso_idx <- which(gene_ids == gn)
  n_g <- length(iso_idx)

  # Give one isoform dominant expression (>50% of gene total)
  dominant_idx <- iso_idx[1]
  base_dominant <- runif(1, 100, 500)
  expr_mat[dominant_idx, ] <- pmax(0.5,
    rnorm(6, mean = base_dominant, sd = base_dominant * 0.15))

  # Others get lower expression
  if (n_g > 1) {
    for (j in iso_idx[-1]) {
      base <- runif(1, 1, 50)
      expr_mat[j, ] <- pmax(0.5,
        rnorm(6, mean = base, sd = base * 0.3))
    }
  }
}

expr_df <- data.frame(isoform_id = isoform_ids, expr_mat,
                      check.names = FALSE, stringsAsFactors = FALSE)
outpath3 <- "inst/extdata/example_expression.csv"
write.csv(expr_df, outpath3, row.names = FALSE)
cat(sprintf("  Wrote: %s (%d isoforms x %d samples)\n",
            outpath3, n_iso, 6))

# -- DE results --
logfc_vals <- rnorm(n_iso, mean = 0, sd = 1.5)
pvals <- runif(n_iso, 0, 1)
# Make ~8 significant (mix up/down)
sig_idx <- sample(n_iso, min(8, n_iso))
logfc_vals[sig_idx] <- sample(c(-1, 1), length(sig_idx), replace = TRUE) *
                       runif(length(sig_idx), 1.2, 4)
pvals[sig_idx] <- runif(length(sig_idx), 1e-8, 0.04)

de_df <- data.frame(
  isoform_id = isoform_ids,
  gene_id = gene_ids,
  logFC = round(logfc_vals, 4),
  adj_p_val = signif(pvals, 4),
  stringsAsFactors = FALSE
)
outpath4 <- "inst/extdata/example_de.csv"
write.csv(de_df, outpath4, row.names = FALSE)
cat(sprintf("  Wrote: %s (%d rows, %d significant)\n",
            outpath4, n_iso, sum(pvals < 0.05)))

# -- DU results --
dpsi_vals <- rnorm(n_iso, mean = 0, sd = 0.15)
du_pvals <- runif(n_iso, 0, 1)
sig_du_idx <- sample(n_iso, min(8, n_iso))
dpsi_vals[sig_du_idx] <- sample(c(-1, 1), length(sig_du_idx),
                                replace = TRUE) *
                         runif(length(sig_du_idx), 0.15, 0.5)
du_pvals[sig_du_idx] <- runif(length(sig_du_idx), 1e-6, 0.04)

du_df <- data.frame(
  isoform_id = isoform_ids,
  gene_id = gene_ids,
  dPSI = round(dpsi_vals, 4),
  adj_p_val = signif(du_pvals, 4),
  stringsAsFactors = FALSE
)
outpath5 <- "inst/extdata/example_du.csv"
write.csv(du_df, outpath5, row.names = FALSE)
cat(sprintf("  Wrote: %s (%d rows, %d significant)\n\n",
            outpath5, n_iso, sum(du_pvals < 0.05)))

# ============================================================================
# Step 7: Pre-computed data objects
# ============================================================================

cat("Creating pre-computed data objects...\n")

devtools::load_all(".")

# Parse structures from example_small
example_structures <- parseIsoformStructures(
  system.file("extdata", "example_small.gtf", package = "Isopair"),
  verbose = FALSE
)
cat(sprintf("  example_structures: %d isoforms\n", nrow(example_structures)))

# Build union exons
ue_result <- buildUnionExons(example_structures, verbose = FALSE)

# Generate pairs (top_two)
example_pairs <- generatePairsExpression(
  expr_mat, gene_map, sample_names, method = "top_two"
)
cat(sprintf("  example_pairs: %d pairs\n", nrow(example_pairs)))

# Build profiles (strand auto-lookup)
example_profiles <- buildProfiles(
  example_pairs, example_structures,
  ue_result$union_exons, ue_result$isoform_union_mapping,
  verbose = FALSE
)
cat(sprintf("  example_profiles: %d profiles\n", nrow(example_profiles)))

# Save
usethis::use_data(example_structures, overwrite = TRUE)
usethis::use_data(example_pairs, overwrite = TRUE)
usethis::use_data(example_profiles, overwrite = TRUE)

cat("\nDone!\n")
