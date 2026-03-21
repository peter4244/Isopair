# Tests for parseIsoformStructures

test_that("parseIsoformStructures returns correct columns", {
  result <- parseIsoformStructures(test_gtf_path, verbose = FALSE)
  expected_cols <- c("isoform_id", "gene_id", "chr", "strand", "n_exons",
                     "exon_starts", "exon_ends", "tx_start", "tx_end",
                     "n_junctions")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("parseIsoformStructures filters by isoform_ids", {
  all_result <- parseIsoformStructures(test_gtf_path, verbose = FALSE)
  # Pick first 5 isoforms
  subset_ids <- head(all_result$isoform_id, 5)
  filtered <- parseIsoformStructures(test_gtf_path, isoform_ids = subset_ids,
                                     verbose = FALSE)
  expect_equal(nrow(filtered), length(subset_ids))
  expect_true(all(filtered$isoform_id %in% subset_ids))
})

test_that("n_junctions equals n_exons - 1", {
  result <- parseIsoformStructures(test_gtf_path, verbose = FALSE)
  expect_true(all(result$n_junctions == result$n_exons - 1))
})

test_that("exon_starts and exon_ends are sorted within each isoform", {
  result <- parseIsoformStructures(test_gtf_path, verbose = FALSE)
  for (i in seq_len(nrow(result))) {
    starts <- result$exon_starts[[i]]
    if (length(starts) > 1) {
      expect_true(all(diff(starts) > 0),
                  info = sprintf("Isoform %s exon_starts not sorted",
                                  result$isoform_id[i]))
    }
  }
})


# ==========================================================================
# deduplicateStructures()
# ==========================================================================

test_that("deduplicateStructures: no duplicates when disjoint sets", {
  s1 <- test_structures[1:3, ]
  s2 <- test_structures[4:6, ]
  # Ensure no isoform_id overlap
  if (any(s1$isoform_id %in% s2$isoform_id)) skip("Overlap in test data")
  result <- deduplicateStructures(s1, s2, verbose = FALSE)
  expect_equal(result$n_duplicates, 0L)
  expect_equal(result$n_total, nrow(s1) + nrow(s2))
  expect_true("source" %in% names(result$structures))
})

test_that("deduplicateStructures: exact duplicates detected (same GTF twice)", {
  s <- test_structures[1:5, ]
  result <- deduplicateStructures(s, s, verbose = FALSE)
  expect_equal(result$n_duplicates, 5L)
  expect_equal(result$n_total, 5L)
  # All kept rows should be from source "a" (default prefer)
  expect_true(all(result$structures$source == "a"))
})

test_that("deduplicateStructures: prefer='b' keeps source b", {
  s <- test_structures[1:3, ]
  result <- deduplicateStructures(s, s, prefer = "b", verbose = FALSE)
  expect_equal(result$n_duplicates, 3L)
  expect_true(all(result$structures$source == "b"))
})

test_that("deduplicateStructures: partial overlap", {
  # s1 has isoforms 1-5, s2 has isoforms 3-7
  n <- nrow(test_structures)
  if (n < 7) skip("Need at least 7 isoforms in test data")
  s1 <- test_structures[1:5, ]
  s2 <- test_structures[3:7, ]
  result <- deduplicateStructures(s1, s2, verbose = FALSE)
  # Isoforms 3-5 are duplicates (3 duplicates removed)
  expect_equal(result$n_duplicates, 3L)
  expect_equal(result$n_total, 7L)
})

test_that("deduplicateStructures: id_mapping has correct columns", {
  s <- test_structures[1:2, ]
  result <- deduplicateStructures(s, s, verbose = FALSE)
  expect_true(all(c("kept_id", "removed_id", "gene_id") %in%
                    names(result$id_mapping)))
  expect_equal(nrow(result$id_mapping), 2L)
})

test_that("deduplicateStructures: empty inputs handled", {
  empty <- test_structures[0, ]
  result <- deduplicateStructures(empty, empty, verbose = FALSE)
  expect_equal(result$n_total, 0L)
  expect_equal(result$n_duplicates, 0L)
})

test_that("deduplicateStructures: one empty, one non-empty", {
  s1 <- test_structures[1:3, ]
  empty <- test_structures[0, ]
  result <- deduplicateStructures(s1, empty, verbose = FALSE)
  expect_equal(result$n_total, 3L)
  expect_equal(result$n_duplicates, 0L)
  expect_true(all(result$structures$source == "a"))
})

test_that("deduplicateStructures: strip_gene_version groups correctly", {
  # Create two structures with versioned vs unversioned gene_ids
  s1 <- test_structures[1:2, ]
  s2 <- s1
  s2$source <- NULL
  # Add version to gene_id in s2
  s2$gene_id <- paste0(s2$gene_id, ".1")
  # Without strip_gene_version: different genes → no duplicates
  result_no_strip <- deduplicateStructures(s1, s2, strip_gene_version = FALSE,
                                            verbose = FALSE)
  expect_equal(result_no_strip$n_duplicates, 0L)
  # With strip_gene_version: same genes → duplicates
  result_strip <- deduplicateStructures(s1, s2, strip_gene_version = TRUE,
                                         verbose = FALSE)
  expect_equal(result_strip$n_duplicates, 2L)
})

test_that("deduplicateStructures: same isoform_id, different structure warns", {
  s1 <- test_structures[1, ]
  s2 <- s1
  s2$source <- NULL
  # Modify exon structure to differ
  s2$exon_starts[[1]] <- s2$exon_starts[[1]] + 100L
  expect_warning(
    deduplicateStructures(s1, s2, verbose = FALSE),
    "different exon structures"
  )
})

test_that("deduplicateStructures: same isoform_id, same structure deduplicates normally", {
  s1 <- test_structures[1, ]
  result <- deduplicateStructures(s1, s1, verbose = FALSE)
  expect_equal(result$n_duplicates, 1L)
  expect_equal(result$n_total, 1L)
})

test_that("deduplicateStructures: source column preserved in output", {
  s1 <- test_structures[1:2, ]
  s2 <- test_structures[3:4, ]
  result <- deduplicateStructures(s1, s2, verbose = FALSE)
  expect_true("source" %in% names(result$structures))
  expect_true(all(result$structures$source %in% c("a", "b")))
  expect_equal(sum(result$structures$source == "a"), 2L)
  expect_equal(sum(result$structures$source == "b"), 2L)
})

test_that("deduplicateStructures: output has same columns as input plus source", {
  s1 <- test_structures[1:3, ]
  s2 <- test_structures[4:6, ]
  result <- deduplicateStructures(s1, s2, verbose = FALSE)
  expected_cols <- c(names(test_structures), "source")
  expect_true(all(expected_cols %in% names(result$structures)))
})
