# Tests for reconstruction functions

test_that("reconstructDominant handles empty events", {
  # When no events: comparator is identical to reference
  exons <- data.frame(
    chr = "chr1", exon_start = c(100L, 300L, 500L),
    exon_end = c(200L, 400L, 600L), strand = "+",
    gene_id = "G1", transcript_id = "T1",
    stringsAsFactors = FALSE
  )
  events <- tibble::tibble()
  result <- reconstructDominant(exons, events)
  expect_equal(nrow(result), 3)
  expect_equal(result$exon_start, c(100L, 300L, 500L))
})

test_that("verifyReconstruction passes for identical exons", {
  exons <- data.frame(
    exon_start = c(100L, 300L), exon_end = c(200L, 400L),
    stringsAsFactors = FALSE
  )
  result <- verifyReconstruction(exons, exons, "+")
  expect_equal(result$status, "PASS")
})

test_that("verifyReconstruction fails for mismatched exon counts", {
  orig <- data.frame(exon_start = c(100L, 300L), exon_end = c(200L, 400L),
                     stringsAsFactors = FALSE)
  recon <- data.frame(exon_start = 100L, exon_end = 200L,
                      stringsAsFactors = FALSE)
  result <- verifyReconstruction(orig, recon, "+")
  expect_equal(result$status, "FAIL")
  expect_true(grepl("mismatch", result$reason))
})

test_that("verifyReconstruction allows TSS tolerance", {
  orig <- data.frame(exon_start = c(100L, 300L), exon_end = c(200L, 400L),
                     stringsAsFactors = FALSE)
  recon <- data.frame(exon_start = c(110L, 300L), exon_end = c(200L, 400L),
                      stringsAsFactors = FALSE)
  # Plus strand: first exon start is TSS, tolerance=20 allows 10bp diff
  result <- verifyReconstruction(orig, recon, "+", tss_tolerance = 20)
  expect_equal(result$status, "PASS")
})
