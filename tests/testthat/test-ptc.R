# Tests for computePtcStatus(), testPtcAssociation(), analyzeFrameDisruption()

test_cds_gtf_path <- system.file("extdata", "test_cds_cases.gtf",
                                  package = "Isopair")

# Parse structures and CDS once for PTC tests
cds_structures <- parseIsoformStructures(test_cds_gtf_path, verbose = FALSE)
cds_metadata <- extractCdsAnnotations(test_cds_gtf_path, verbose = FALSE)

# ==========================================================================
# computePtcStatus()
# ==========================================================================

test_that("computePtcStatus: C1 has_ptc=TRUE, distance=100", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  c1 <- ptc[ptc$isoform_id == "TX_C1_1", ]
  expect_equal(nrow(c1), 1L)
  expect_true(c1$has_ptc)
  expect_equal(c1$ptc_distance, 100)
})

test_that("computePtcStatus: C2 has_ptc=FALSE, distance=30", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  c2 <- ptc[ptc$isoform_id == "TX_C2_1", ]
  expect_false(c2$has_ptc)
  expect_equal(c2$ptc_distance, 30)
})

test_that("computePtcStatus: C3 stop_in_last_exon=TRUE, has_ptc=FALSE", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  c3 <- ptc[ptc$isoform_id == "TX_C3_1", ]
  expect_true(c3$stop_in_last_exon)
  expect_false(c3$has_ptc)
  # ptc_distance should be negative (stop is downstream of last EJC)
  expect_true(c3$ptc_distance < 0)
})

test_that("computePtcStatus: C4 single-exon has_ptc=FALSE", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  c4 <- ptc[ptc$isoform_id == "TX_C4_1", ]
  expect_false(c4$has_ptc)
  expect_true(is.na(c4$ptc_distance))
  expect_true(c4$stop_in_last_exon)
  expect_equal(c4$n_downstream_ejcs, 0L)
})

test_that("computePtcStatus: minus-strand gene F2 has_ptc=TRUE, distance=100", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  f2 <- ptc[ptc$isoform_id == "TX_F2_1", ]
  expect_true(f2$has_ptc)
  expect_equal(f2$ptc_distance, 100)
})

test_that("computePtcStatus: non-coding excluded from output", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  expect_false("TX_G1_1" %in% ptc$isoform_id)
})

test_that("computePtcStatus: empty input → empty output", {
  empty_struct <- cds_structures[0L, ]
  empty_cds <- cds_metadata[0L, ]
  ptc <- computePtcStatus(empty_struct, empty_cds)
  expect_equal(nrow(ptc), 0L)
  expect_true(all(c("isoform_id", "has_ptc", "ptc_distance") %in% names(ptc)))
})

test_that("computePtcStatus: output has correct columns", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  expected_cols <- c("isoform_id", "gene_id", "n_exons", "orf_length",
                     "ptc_distance", "has_ptc", "n_downstream_ejcs",
                     "stop_in_last_exon")
  expect_true(all(expected_cols %in% names(ptc)))
})

test_that("computePtcStatus: n_downstream_ejcs correct for C1", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  c1 <- ptc[ptc$isoform_id == "TX_C1_1", ]
  # Stop is in exon 2 (of 3), so 1 downstream EJC
  expect_equal(c1$n_downstream_ejcs, 1L)
})

# ==========================================================================
# testPtcAssociation()
# ==========================================================================

test_that("testPtcAssociation: basic McNemar test", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)

  # Create pairs where C1 (PTC+) is comparator, A1 (PTC-) is reference
  pairs <- data.frame(
    gene_id = "test_gene",
    reference_isoform_id = "TX_A1_1",
    comparator_isoform_id = "TX_C1_1",
    stringsAsFactors = FALSE
  )
  result <- testPtcAssociation(ptc, pairs)

  expect_true(is.list(result))
  expect_true("pair_summary" %in% names(result))
  expect_true("comp_ptc_rate" %in% names(result))
  expect_true("ref_ptc_rate" %in% names(result))
  expect_equal(result$comp_ptc_rate, 1.0)
})

test_that("testPtcAssociation: all-concordant pairs", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)

  # Both reference and comparator have PTC
  pairs <- data.frame(
    gene_id = "test_gene",
    reference_isoform_id = "TX_C1_1",
    comparator_isoform_id = "TX_C1_1",
    stringsAsFactors = FALSE
  )
  result <- testPtcAssociation(ptc, pairs)
  # McNemar should be NULL (no discordant pairs)
  expect_null(result$mcnemar_test)
})

test_that("testPtcAssociation: with group_col stratification", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)

  pairs <- data.frame(
    gene_id = c("g1", "g2"),
    reference_isoform_id = c("TX_A1_1", "TX_A1_1"),
    comparator_isoform_id = c("TX_C1_1", "TX_C2_1"),
    group = c("grp_a", "grp_b"),
    stringsAsFactors = FALSE
  )
  result <- testPtcAssociation(ptc, pairs, group_col = "group")
  expect_true("group_tests" %in% names(result))
})

# ==========================================================================
# analyzeFrameDisruption()
# ==========================================================================

# Helper to build minimal profiles with detailed_events
.make_fd_profiles <- function(structures, ref_id, comp_id, gene_id, strand) {
  ref_struct <- structures[structures$isoform_id == ref_id, ]
  comp_struct <- structures[structures$isoform_id == comp_id, ]

  ref_exons <- data.frame(
    exon_start = ref_struct$exon_starts[[1]],
    exon_end = ref_struct$exon_ends[[1]],
    chr = ref_struct$chr, strand = ref_struct$strand,
    gene_id = ref_struct$gene_id,
    transcript_id = ref_id,
    exon_number = seq_along(ref_struct$exon_starts[[1]])
  )
  comp_exons <- data.frame(
    exon_start = comp_struct$exon_starts[[1]],
    exon_end = comp_struct$exon_ends[[1]],
    chr = comp_struct$chr, strand = comp_struct$strand,
    gene_id = comp_struct$gene_id,
    transcript_id = comp_id,
    exon_number = seq_along(comp_struct$exon_starts[[1]])
  )

  events <- detectEvents(ref_exons, comp_exons,
                          gene_id, ref_id, comp_id, strand)

  tibble::tibble(
    gene_id = gene_id,
    reference_isoform_id = ref_id,
    comparator_isoform_id = comp_id,
    detailed_events = list(events)
  )
}

test_that("analyzeFrameDisruption: D1 not frame-disrupting (90bp, 3n)", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D1_ref", "TX_D1_comp",
                                 "GENE_FD_D1", "+")
  fd <- analyzeFrameDisruption(profiles, cds_metadata)

  expect_true(is.list(fd))
  expect_true("events" %in% names(fd))
  expect_true("profile_summary" %in% names(fd))

  summary <- fd$profile_summary
  expect_equal(summary$has_frame_disruption[1], FALSE)
})

test_that("analyzeFrameDisruption: D2 frame-disrupting (91bp, 3n+1)", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D2_ref", "TX_D2_comp",
                                 "GENE_FD_D2", "+")
  fd <- analyzeFrameDisruption(profiles, cds_metadata)

  summary <- fd$profile_summary
  expect_equal(summary$has_frame_disruption[1], TRUE)
  expect_true(summary$n_frame_disrupting[1] > 0L)
})

test_that("analyzeFrameDisruption: D3 A5SS event frame disruption", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D3_ref", "TX_D3_comp",
                                 "GENE_FD_D3", "+")
  fd <- analyzeFrameDisruption(profiles, cds_metadata)

  # Check that events include frame disruption info
  events <- fd$events
  expect_true("frame_disrupting" %in% names(events))
  expect_true("cds_overlap_bp" %in% names(events))
  expect_true("affects_cds" %in% names(events))

  summary <- fd$profile_summary
  expect_true(summary$has_frame_disruption[1])
})

test_that("analyzeFrameDisruption: non-CDS events → affects_cds=FALSE", {
  # Use a pair where ref is non-coding (no CDS metadata)
  profiles <- tibble::tibble(
    gene_id = "GENE_NC_G1",
    reference_isoform_id = "TX_G1_1",
    comparator_isoform_id = "TX_G1_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 42300L, three_prime = 42600L,
      bp_diff = 300L
    ))
  )
  fd <- analyzeFrameDisruption(profiles, cds_metadata)
  expect_false(fd$events$affects_cds[1])
  expect_false(fd$events$frame_disrupting[1])
})

test_that("analyzeFrameDisruption: empty events handled", {
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A1",
    reference_isoform_id = "TX_A1_1",
    comparator_isoform_id = "TX_A1_1",
    detailed_events = list(tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ))
  )
  fd <- analyzeFrameDisruption(profiles, cds_metadata)
  expect_equal(fd$profile_summary$n_cds_events[1], 0L)
  expect_false(fd$profile_summary$has_frame_disruption[1])
})

test_that("analyzeFrameDisruption: minus-strand F3 frame-disrupting", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_F3_ref", "TX_F3_comp",
                                 "GENE_FD_F3", "-")
  fd <- analyzeFrameDisruption(profiles, cds_metadata)
  summary <- fd$profile_summary
  expect_true(summary$has_frame_disruption[1])
})

test_that("analyzeFrameDisruption: IR events excluded from frameshift", {
  # Synthetic profile with an IR event that overlaps CDS
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A1",
    reference_isoform_id = "TX_A1_1",
    comparator_isoform_id = "TX_A1_1",
    detailed_events = list(tibble::tibble(
      event_type = "IR", direction = "GAIN",
      five_prime = 1200L, three_prime = 1500L,
      bp_diff = NA_integer_
    ))
  )
  fd <- analyzeFrameDisruption(profiles, cds_metadata)
  # IR is not in frameshift_types, so n_frame_disrupting = 0
  expect_equal(fd$profile_summary$n_frame_disrupting[1], 0L)
  # But affects_cds can still be TRUE
  expect_true(fd$events$affects_cds[1])
  expect_false(fd$events$frame_disrupting[1])
})

# ==========================================================================
# Integration: end-to-end with real gene GTF (SRSF3 if available)
# ==========================================================================

test_that("Integration: CDS extraction from test_real_genes.gtf works", {
  real_gtf <- system.file("extdata", "test_real_genes.gtf", package = "Isopair")
  if (nchar(real_gtf) > 0L && file.exists(real_gtf)) {
    cds <- extractCdsAnnotations(real_gtf, verbose = FALSE)
    # Should find at least some coding isoforms (TTN, SRSF3, etc.)
    if (nrow(cds) > 0L) {
      coding <- cds[cds$coding_status == "coding", ]
      expect_true(nrow(coding) > 0L)
      expect_true(all(coding$cds_start < coding$cds_stop))
      expect_true(all(coding$orf_length > 0L))
    }
  } else {
    skip("test_real_genes.gtf not available")
  }
})
