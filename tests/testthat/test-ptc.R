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
# analyzeFrameWalk()
# ==========================================================================

test_that("analyzeFrameWalk: H1 in-frame multi-exon skip", {
  # H1: LOSS 4bp + LOSS 2bp in same comparator intron → single splice group
  # Net CDS change: -6, 6%%3 = 0 → in-frame, no frameshift
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_H1_ref", "TX_H1_comp",
                                 "GENE_FW_H1", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_true(fw$profile_summary$frame_resolved[1])
  expect_equal(fw$profile_summary$n_frameshifts[1], 0L)
  expect_equal(fw$profile_summary$n_compensatory[1], 0L)
  # Single splice group
  expect_equal(length(unique(fw$events$splice_group)), 1L)
})

test_that("analyzeFrameWalk: H2 unresolved frameshift", {
  # H2: single 4bp LOSS, offset 2, no compensation
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_H2_ref", "TX_H2_comp",
                                 "GENE_FW_H2", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_false(fw$profile_summary$frame_resolved[1])
  expect_equal(fw$profile_summary$n_frameshifts[1], 1L)
  expect_equal(fw$profile_summary$n_compensatory[1], 0L)
  # frameshifted_bp should be positive
  expect_true(fw$profile_summary$total_frameshifted_cds_bp[1] > 0L)
  expect_true(fw$profile_summary$pct_orf_frameshifted[1] > 0)
})

test_that("analyzeFrameWalk: H3 multi-exon skip within single intron", {
  # H3: three Missing_Internal events all within one comparator intron
  # (same comp_junctions) → single splice group. Net CDS change: -4-2-5 = -11,
  # 11%%3 = 2 → one frameshift, zero compensatory.
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_H3_ref", "TX_H3_comp",
                                 "GENE_FW_H3", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_equal(fw$profile_summary$n_frameshifts[1], 1L)
  expect_equal(fw$profile_summary$n_compensatory[1], 0L)
  expect_false(fw$profile_summary$frame_resolved[1])
  # All events should share one splice group
  expect_equal(length(unique(fw$events$splice_group)), 1L)
})

test_that("analyzeFrameWalk: in-frame SE (6bp) no frameshift", {
  # D1: 90bp SE, 90%%3=0 → no frameshift
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D1_ref", "TX_D1_comp",
                                 "GENE_FD_D1", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  if (fw$profile_summary$n_cds_events[1] > 0L) {
    expect_equal(fw$profile_summary$n_frameshifts[1], 0L)
    expect_true(fw$profile_summary$frame_resolved[1])
  }
})

test_that("analyzeFrameWalk: D2 single frameshift (91bp)", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D2_ref", "TX_D2_comp",
                                 "GENE_FD_D2", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_equal(fw$profile_summary$n_frameshifts[1], 1L)
  expect_false(fw$profile_summary$frame_resolved[1])
})

test_that("analyzeFrameWalk: signed arithmetic prevents false compensation", {
  # Construct profile manually: GAIN of 4bp then LOSS of 5bp
  # Cumulative: +4 then +4-5 = -1, offset = ((-1%%3)+3)%%3 = 2
  # Unsigned would give 4+5=9, 9%%3=0 (incorrectly compensated)
  # These are independent splicing events (different junctions)
  profiles <- tibble::tibble(
    gene_id = "GENE_FD_D2",
    reference_isoform_id = "TX_D2_ref",
    comparator_isoform_id = "TX_D2_comp",
    detailed_events = list(tibble::tibble(
      event_type = c("SE", "A5SS"),
      direction = c("GAIN", "LOSS"),
      five_prime = c(25400L, 25495L),
      three_prime = c(25403L, 25490L),
      bp_diff = c(4L, 5L),
      ref_junctions = c("25300:25400", "25490:25600"),
      comp_junctions = c("25300:25403,25403:25500", "25495:25600")
    ))
  )
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  # Should NOT be compensated — final offset should be non-zero
  expect_false(fw$profile_summary$frame_resolved[1])
  # Two independent splice groups
  expect_equal(length(unique(fw$events$splice_group)), 2L)
})

test_that("analyzeFrameWalk: minus-strand ordering (H4)", {
  # H4: minus strand, events at 60600 (more 5') and 60400 (more 3')
  # 5'→3' on minus = descending genomic → 60600 first, then 60400
  # Both events share same comp intron → single splice group
  # Net: 2bp + 4bp = 6bp, 6%%3 = 0 → in-frame, no frameshift
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_H4_ref", "TX_H4_comp",
                                 "GENE_FW_H4", "-")
  fw <- analyzeFrameWalk(profiles, cds_metadata)

  if (nrow(fw$events) >= 2L) {
    # First event (5' on minus strand) should be the one at higher coords
    expect_true(fw$events$genomic_start[1] >= fw$events$genomic_start[2])
    # Same splice group
    expect_equal(length(unique(fw$events$splice_group)), 1L)
  }
  expect_true(fw$profile_summary$frame_resolved[1])
  expect_equal(fw$profile_summary$n_frameshifts[1], 0L)
})

test_that("analyzeFrameWalk: UTR-only events excluded", {
  # Events entirely in UTR → n_cds_events = 0
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 5010L, three_prime = 5050L,
      bp_diff = 40L
    ))
  )
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_equal(fw$profile_summary$n_cds_events[1], 0L)
})

test_that("analyzeFrameWalk: PTC cross-reference works", {
  ptc <- computePtcStatus(cds_structures, cds_metadata)
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D2_ref", "TX_D2_comp",
                                 "GENE_FD_D2", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata, ptc_table = ptc)
  # comp is TX_D2_comp; should have a non-NA value if it was in PTC table
  expect_false(is.na(fw$profile_summary$comparator_has_ptc[1]))
})

test_that("analyzeFrameWalk: PTC NA when no ptc_table", {
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_D2_ref", "TX_D2_comp",
                                 "GENE_FD_D2", "+")
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_true(is.na(fw$profile_summary$comparator_has_ptc[1]))
})

test_that("analyzeFrameWalk: non-coding reference → skip", {
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
  fw <- analyzeFrameWalk(profiles, cds_metadata)
  expect_equal(fw$profile_summary$n_cds_events[1], 0L)
  expect_true(is.na(fw$profile_summary$frame_resolved[1]))
})

test_that("analyzeFrameWalk: works without coding_status column", {
  # Strip coding_status to simulate manually constructed CDS metadata
  cds_no_status <- cds_metadata[cds_metadata$coding_status == "coding",
                                 c("isoform_id", "cds_start", "cds_stop",
                                   "orf_length", "strand")]
  profiles <- .make_fd_profiles(cds_structures,
                                 "TX_H2_ref", "TX_H2_comp",
                                 "GENE_FW_H2", "+")
  fw_with <- analyzeFrameWalk(profiles, cds_metadata)
  fw_without <- analyzeFrameWalk(profiles, cds_no_status)
  expect_equal(fw_with$profile_summary$n_frameshifts,
               fw_without$profile_summary$n_frameshifts)
  expect_equal(fw_with$profile_summary$frame_resolved,
               fw_without$profile_summary$frame_resolved)
  expect_equal(fw_with$events$cumulative_frame_offset,
               fw_without$events$cumulative_frame_offset)
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


# =============================================================================
# extractCdsExons tests
# =============================================================================

test_that("extractCdsExons returns per-exon CDS intervals", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)

  expect_true(nrow(cds_exons) > 0L)
  expect_true(all(c("isoform_id", "cds_exon_start", "cds_exon_end",
                     "strand", "cds_exon_rank") %in% names(cds_exons)))

  # TX_A1_1 has 3 CDS exons
  a1 <- cds_exons[cds_exons$isoform_id == "TX_A1_1", ]
  expect_equal(nrow(a1), 3L)
  expect_equal(a1$cds_exon_rank, c(1L, 2L, 3L))
})

test_that("extractCdsExons handles isoform_ids filter", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, isoform_ids = "TX_A1_1", verbose = FALSE)
  expect_equal(length(unique(cds_exons$isoform_id)), 1L)
  expect_equal(unique(cds_exons$isoform_id), "TX_A1_1")
})

test_that("extractCdsExons: minus strand rank is 5'→3'", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)

  # TX_J7_ref is minus strand, 3 CDS exons
  j7 <- cds_exons[cds_exons$isoform_id == "TX_J7_ref", ]
  expect_equal(nrow(j7), 3L)
  # Rank 1 should be the highest-coordinate exon (5' end on - strand)
  expect_equal(j7$cds_exon_rank[j7$cds_exon_start == 76700], 1L)
  expect_equal(j7$cds_exon_rank[j7$cds_exon_start == 76500], 2L)
  expect_equal(j7$cds_exon_rank[j7$cds_exon_start == 76000], 3L)
})


# =============================================================================
# compareIsoformFrames tests
# =============================================================================

test_that("compareIsoformFrames: same start, in-frame (J1)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J1",
    reference_isoform_id = "TX_J1_ref",
    comparator_isoform_id = "TX_J1_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(nrow(result$pair_summary), 1L)
  expect_equal(result$pair_summary$frame_category, "same_start_in_frame")
  expect_true(result$pair_summary$same_start_codon)
  expect_equal(result$pair_summary$pct_shared_cds_in_frame, 100)
  expect_equal(result$pair_summary$pct_shared_cds_frameshifted, 0)
  expect_true(all(result$region_detail$phase_match))
})

test_that("compareIsoformFrames: same start, frameshift (J2)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J2",
    reference_isoform_id = "TX_J2_ref",
    comparator_isoform_id = "TX_J2_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "same_start_frameshift")
  expect_true(result$pair_summary$same_start_codon)
  # First shared region should match (upstream of skip), second should not
  expect_equal(nrow(result$region_detail), 2L)
  expect_true(result$region_detail$phase_match[1])   # [71000-71100] both phase 0
  expect_false(result$region_detail$phase_match[2])   # [71500-71700] phase differs
})

test_that("compareIsoformFrames: diff start, frame preserved (J3)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J3",
    reference_isoform_id = "TX_J3_ref",
    comparator_isoform_id = "TX_J3_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "diff_start_same_frame")
  expect_false(result$pair_summary$same_start_codon)
  expect_equal(result$pair_summary$pct_shared_cds_in_frame, 100)
})

test_that("compareIsoformFrames: diff start, frame divergent (J4)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J4",
    reference_isoform_id = "TX_J4_ref",
    comparator_isoform_id = "TX_J4_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "diff_start_diff_frame")
  expect_false(result$pair_summary$same_start_codon)
  expect_true(result$pair_summary$pct_shared_cds_frameshifted > 0)
})

test_that("compareIsoformFrames: non-coding pair (J5)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf,
    isoform_ids = c("TX_J5_ref", "TX_J5_comp"), verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J5",
    reference_isoform_id = "TX_J5_ref",
    comparator_isoform_id = "TX_J5_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "non_coding")
})

test_that("compareIsoformFrames: no shared CDS (J6)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J6",
    reference_isoform_id = "TX_J6_ref",
    comparator_isoform_id = "TX_J6_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "no_shared_cds")
  expect_equal(result$pair_summary$n_shared_cds_regions, 0L)
})

test_that("compareIsoformFrames: minus strand frameshift (J7)", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J7",
    reference_isoform_id = "TX_J7_ref",
    comparator_isoform_id = "TX_J7_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(result$pair_summary$frame_category, "same_start_frameshift")
  expect_true(result$pair_summary$same_start_codon)
  expect_equal(nrow(result$region_detail), 2L)
  # First region (5' on - strand = higher coordinates) should match
  # Second region should not match
  high_region <- result$region_detail[result$region_detail$region_start == 76700, ]
  low_region <- result$region_detail[result$region_detail$region_start == 76000, ]
  expect_true(high_region$phase_match)
  expect_false(low_region$phase_match)
})

test_that("compareIsoformFrames: empty input", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = character(0),
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_equal(nrow(result$pair_summary), 0L)
  expect_equal(nrow(result$region_detail), 0L)
})

test_that("compareIsoformFrames: region_detail has correct columns and phase values", {
  gtf <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
  cds_exons <- extractCdsExons(gtf, verbose = FALSE)
  cds <- extractCdsAnnotations(gtf, verbose = FALSE)

  pairs <- data.frame(
    gene_id = "GENE_FRAME_J2",
    reference_isoform_id = "TX_J2_ref",
    comparator_isoform_id = "TX_J2_comp",
    stringsAsFactors = FALSE
  )
  result <- compareIsoformFrames(pairs, cds_exons, cds)

  expect_true(all(c("region_start", "region_end", "region_bp",
                     "ref_phase_at_start", "comp_phase_at_start",
                     "phase_match") %in% names(result$region_detail)))
  expect_true(all(result$region_detail$ref_phase_at_start %in% 0:2))
  expect_true(all(result$region_detail$comp_phase_at_start %in% 0:2))
})
