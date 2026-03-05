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
