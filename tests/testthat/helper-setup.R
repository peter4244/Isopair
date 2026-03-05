# Shared test fixtures for Isopair tests
#
# Loads GTF, builds structures, builds union exons once for all tests.
# Available to all test files automatically via testthat helper loading.

# Path to test data
test_gtf_path <- system.file("extdata", "test_cases.gtf", package = "Isopair")
test_pairs_path <- system.file("extdata", "test_pairs.tsv", package = "Isopair")

# Parse structures from test GTF
test_structures <- parseIsoformStructures(test_gtf_path, verbose = FALSE)

# Build union exons
test_ue_result <- buildUnionExons(test_structures, verbose = FALSE)
test_union_exons <- test_ue_result$union_exons
test_isoform_union_mapping <- test_ue_result$isoform_union_mapping

# Parse test pairs
test_pairs_raw <- utils::read.delim(test_pairs_path, header = TRUE,
                                    stringsAsFactors = FALSE)

# Expand structures to per-exon rows (for event detection)
.expandStructuresForTest <- function(structures) {
  rows <- list()
  for (i in seq_len(nrow(structures))) {
    s <- structures[i, ]
    starts <- s$exon_starts[[1]]
    ends <- s$exon_ends[[1]]
    rows[[i]] <- data.frame(
      isoform_id = s$isoform_id,
      gene_id = s$gene_id,
      chr = s$chr,
      strand = s$strand,
      exon_number = seq_along(starts),
      exon_start = as.integer(starts),
      exon_end = as.integer(ends),
      transcript_id = s$isoform_id,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

test_exons <- .expandStructuresForTest(test_structures)
