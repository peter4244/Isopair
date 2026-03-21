# Tests for testPositionalBias(), classifyTopology(), testProximity()

# Helper: build profiles with detailed events for spatial tests
.make_spatial_profiles <- function() {
  # Gene spans 1000-5000, + strand
  # Create profiles with various event types and positions
  profiles <- list()

  # Profile 1: A5SS at position ~1500 + A3SS at ~2500 (non-adjacent → Distal)
  profiles[[1]] <- tibble::tibble(
    gene_id = "GENE_S1",
    reference_isoform_id = "REF_S1",
    comparator_isoform_id = "COMP_S1",
    n_a5ss = 1L, n_a3ss = 1L, n_se = 0L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 0L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "A3SS"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1480L, 2480L),
      three_prime = c(1520L, 2520L),
      bp_diff = c(40L, 40L)
    ))
  )

  # Profile 2: A5SS at exon 2 end (~2000) + A3SS at exon 3 start (~2200) → F2F
  # (consecutive exons on + strand)
  profiles[[2]] <- tibble::tibble(
    gene_id = "GENE_S2",
    reference_isoform_id = "REF_S2",
    comparator_isoform_id = "COMP_S2",
    n_a5ss = 1L, n_a3ss = 1L, n_se = 0L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 0L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "A3SS"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1995L, 2205L),
      three_prime = c(2005L, 2195L),
      bp_diff = c(10L, 10L)
    ))
  )

  # Profile 3: A5SS and A3SS on same exon (exon 1: 1000-1200) → B2B
  # A5SS midpoint=1200 → nearest exon_end=1200 (exon 1, idx 1)
  # A3SS midpoint=1000 → nearest exon_start=1000 (exon 1, idx 1)
  profiles[[3]] <- tibble::tibble(
    gene_id = "GENE_S3",
    reference_isoform_id = "REF_S3",
    comparator_isoform_id = "COMP_S3",
    n_a5ss = 1L, n_a3ss = 1L, n_se = 0L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 0L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "A3SS"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1195L, 1005L),
      three_prime = c(1205L, 995L),
      bp_diff = c(10L, 10L)
    ))
  )

  # Profile 4: Only SE events (no A5SS/A3SS) — for proximity
  profiles[[4]] <- tibble::tibble(
    gene_id = "GENE_S4",
    reference_isoform_id = "REF_S4",
    comparator_isoform_id = "COMP_S4",
    n_a5ss = 0L, n_a3ss = 0L, n_se = 2L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 0L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("SE", "SE"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1500L, 1600L),
      three_prime = c(1550L, 1650L),
      bp_diff = c(50L, 50L)
    ))
  )

  dplyr::bind_rows(profiles)
}

# Structures matching the spatial profiles
.make_spatial_structures <- function() {
  # Each gene has 4 exons spanning 1000-5000
  genes <- c("GENE_S1", "GENE_S2", "GENE_S3", "GENE_S4")
  ref_ids <- c("REF_S1", "REF_S2", "REF_S3", "REF_S4")

  tibble::tibble(
    isoform_id = ref_ids,
    gene_id = genes,
    chr = "chrT",
    strand = "+",
    n_exons = 4L,
    exon_starts = list(
      c(1000L, 1500L, 2200L, 3000L),
      c(1000L, 1500L, 2200L, 3000L),
      c(1000L, 1500L, 2200L, 3000L),
      c(1000L, 1500L, 2200L, 3000L)
    ),
    exon_ends = list(
      c(1200L, 2000L, 2800L, 4000L),
      c(1200L, 2000L, 2800L, 4000L),
      c(1200L, 2000L, 2800L, 4000L),
      c(1200L, 2000L, 2800L, 4000L)
    ),
    tx_start = 1000L,
    tx_end = 4000L,
    n_junctions = 3L
  )
}

# ==========================================================================
# testPositionalBias()
# ==========================================================================

test_that("testPositionalBias returns correct columns", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- testPositionalBias(profiles, structures, min_events = 1L)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("event_type", "n_events", "median_position",
                     "mean_position", "ks_statistic", "p_value",
                     "adj_p_value", "significant",
                     "bias_direction") %in% names(result)))
})

test_that("testPositionalBias: relative positions in [0,1]", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  events <- Isopair:::.extractEventsWithPositions(profiles, structures)
  expect_true(all(events$relative_pos >= 0))
  expect_true(all(events$relative_pos <= 1))
})

test_that("testPositionalBias: min_events filtering works", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  # With min_events = 100, should skip all types
  result <- testPositionalBias(profiles, structures, min_events = 100L)
  expect_equal(nrow(result), 0L)
})

test_that("testPositionalBias: FDR correction applied", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- testPositionalBias(profiles, structures, min_events = 1L)
  if (nrow(result) > 0L) {
    expect_true(all(result$adj_p_value >= result$p_value, na.rm = TRUE))
  }
})

# ==========================================================================
# classifyTopology()
# ==========================================================================

test_that("classifyTopology returns correct structure", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_true(is.list(result))
  expect_true("classifications" %in% names(result))
  expect_true("summary" %in% names(result))
})

test_that("classifyTopology: F2F classification correct", {
  # Profile 2: A5SS near exon 2 end, A3SS near exon 3 start → F2F
  profiles <- .make_spatial_profiles()[2, ]
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(nrow(result$classifications), 1L)
  expect_equal(result$classifications$topology[1], "F2F")
})

test_that("classifyTopology: B2B classification correct", {
  # Profile 3: A5SS and A3SS both near exon 1 → B2B
  profiles <- .make_spatial_profiles()[3, ]
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(nrow(result$classifications), 1L)
  expect_equal(result$classifications$topology[1], "B2B")
})

test_that("classifyTopology: Distal classification correct", {
  # Profile 1: A5SS on exon 2, A3SS on exon 3 → but non-adjacent = Distal
  # A5SS at 1500 (near exon 2 end=2000), A3SS at 2500 (near exon 3 start=2200)
  # For + strand: A5SS nearest exon_end = 2000 (exon 2, idx 2),
  #               A3SS nearest exon_start = 2200 (exon 3, idx 3)
  # a5ss_exon_idx=2, a3ss_exon_idx=3, 2 == 3-1 → F2F actually!
  # Let me check again... Profile 1: five_prime=1480, three_prime=1520 → mid=1500
  # nearest exon_end to 1500: |1200-1500|=300, |2000-1500|=500, |2800-1500|=1300
  # → exon 1 (idx 1)
  # A3SS: five_prime=2480, three_prime=2520 → mid=2500
  # nearest exon_start to 2500: |1000-2500|=1500, |1500-2500|=1000, |2200-2500|=300, |3000-2500|=500
  # → exon 3 (idx 3)
  # a5ss_idx=1, a3ss_idx=3, not same, not consecutive → Distal ✓
  profiles <- .make_spatial_profiles()[1, ]
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(result$classifications$topology[1], "Distal")
})

test_that("classifyTopology: skips profiles without both A5SS and A3SS", {
  # Profile 4 has only SE events
  profiles <- .make_spatial_profiles()[4, ]
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(nrow(result$classifications), 0L)
})

test_that("classifyTopology: summary proportions sum to 1", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  if (nrow(result$summary) > 0L) {
    expect_equal(sum(result$summary$proportion), 1.0, tolerance = 1e-10)
  }
})

test_that("classifyTopology: permutation test runs", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = TRUE, n_perm = 50L)
  expect_true("permutation_p_value" %in% names(result$summary))
  expect_true(all(result$summary$permutation_p_value >= 0))
  expect_true(all(result$summary$permutation_p_value <= 1))
})

test_that("classifyTopology: empty profiles → empty result", {
  profiles <- tibble::tibble(
    gene_id = character(0), reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    detailed_events = list()
  )
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(nrow(result$classifications), 0L)
})

test_that("classifyTopology: Partial_IR_5 + Partial_IR_3 classified with expanded types", {
  # Reuse GENE_S2 structure (4 exons, + strand, exon2 end=2000, exon3 start=2200)
  # Partial_IR_5 near exon 2 end → F2F with Partial_IR_3 near exon 3 start
  profiles <- tibble::tibble(
    gene_id = "GENE_S2",
    reference_isoform_id = "REF_S2",
    comparator_isoform_id = "COMP_S2_PIR",
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 2L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("Partial_IR_5", "Partial_IR_3"),
      direction = c("GAIN", "GAIN"),
      five_prime = c(1995L, 2205L),
      three_prime = c(2005L, 2195L),
      bp_diff = c(10L, 10L)
    ))
  )
  structures <- .make_spatial_structures()

  # Default types (A5SS/A3SS only) should skip this profile
  result_default <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(nrow(result_default$classifications), 0L)

  # Expanded types should classify it as F2F
  result_expanded <- classifyTopology(profiles, structures,
    five_prime_types = c("A5SS", "Partial_IR_5"),
    three_prime_types = c("A3SS", "Partial_IR_3"),
    test = FALSE)
  expect_equal(nrow(result_expanded$classifications), 1L)
  expect_equal(result_expanded$classifications$topology[1], "F2F")
  expect_equal(result_expanded$classifications$event_a_type[1], "Partial_IR_5")
  expect_equal(result_expanded$classifications$event_b_type[1], "Partial_IR_3")
})

test_that("classifyTopology: mixed A5SS + Partial_IR_3 classified correctly", {
  # A5SS near exon 2 end + Partial_IR_3 near exon 3 start → F2F
  profiles <- tibble::tibble(
    gene_id = "GENE_S2",
    reference_isoform_id = "REF_S2",
    comparator_isoform_id = "COMP_S2_MIX",
    n_a5ss = 1L, n_a3ss = 0L, n_se = 0L, n_ir = 0L, n_ir_diff = 0L,
    n_partial_ir = 1L, n_missing_internal = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L, n_events = 2L,
    tss_changed = FALSE, tes_changed = FALSE,
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "Partial_IR_3"),
      direction = c("LOSS", "GAIN"),
      five_prime = c(1995L, 2205L),
      three_prime = c(2005L, 2195L),
      bp_diff = c(10L, 10L)
    ))
  )
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures,
    five_prime_types = c("A5SS", "Partial_IR_5"),
    three_prime_types = c("A3SS", "Partial_IR_3"),
    test = FALSE)
  expect_equal(nrow(result$classifications), 1L)
  expect_equal(result$classifications$topology[1], "F2F")
  expect_equal(result$classifications$event_a_type[1], "A5SS")
  expect_equal(result$classifications$event_b_type[1], "Partial_IR_3")
})

test_that("classifyTopology: actual event types stored in output", {
  profiles <- .make_spatial_profiles()[2, ]  # A5SS + A3SS → F2F
  structures <- .make_spatial_structures()
  result <- classifyTopology(profiles, structures, test = FALSE)
  expect_equal(result$classifications$event_a_type[1], "A5SS")
  expect_equal(result$classifications$event_b_type[1], "A3SS")
})

# ==========================================================================
# testProximity()
# ==========================================================================

test_that("testProximity returns correct structure", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- testProximity(profiles, structures, n_perm = 50L)
  expect_true(is.list(result))
  expect_true(all(c("z_score", "p_value", "observed_mean_distance",
                     "expected_mean_distance", "n_profiles_used",
                     "n_profiles_excluded") %in% names(result)))
})

test_that("testProximity: uses only internal events", {
  # Profile with only terminal events should be excluded
  profiles <- tibble::tibble(
    gene_id = "GENE_T",
    reference_isoform_id = "REF_T",
    comparator_isoform_id = "COMP_T",
    detailed_events = list(tibble::tibble(
      event_type = c("Alt_TSS", "Alt_TES"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1000L, 4000L),
      three_prime = c(1050L, 3950L),
      bp_diff = c(50L, 50L)
    ))
  )
  structures <- tibble::tibble(
    isoform_id = "REF_T", gene_id = "GENE_T", chr = "chrT",
    strand = "+", n_exons = 2L,
    exon_starts = list(c(1000L, 3000L)),
    exon_ends = list(c(2000L, 4000L)),
    tx_start = 1000L, tx_end = 4000L, n_junctions = 1L
  )
  result <- testProximity(profiles, structures, n_perm = 50L)
  expect_equal(result$n_profiles_used, 0L)
})

test_that("testProximity: clustered events produce negative z-score", {
  # Two SE events very close together in many profiles
  n <- 20
  profiles <- tibble::tibble(
    gene_id = paste0("G", seq_len(n)),
    reference_isoform_id = paste0("R", seq_len(n)),
    comparator_isoform_id = paste0("C", seq_len(n)),
    detailed_events = replicate(n, tibble::tibble(
      event_type = c("SE", "SE"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(2000L, 2010L),
      three_prime = c(2050L, 2060L),
      bp_diff = c(50L, 50L)
    ), simplify = FALSE)
  )
  structures <- tibble::tibble(
    isoform_id = paste0("R", seq_len(n)),
    gene_id = paste0("G", seq_len(n)),
    chr = "chrT", strand = "+", n_exons = 4L,
    exon_starts = replicate(n, c(1000L, 1500L, 2200L, 3000L), simplify = FALSE),
    exon_ends = replicate(n, c(1200L, 2000L, 2800L, 4000L), simplify = FALSE),
    tx_start = 1000L, tx_end = 4000L, n_junctions = 3L
  )
  set.seed(123)
  result <- testProximity(profiles, structures, n_perm = 200L)
  expect_true(result$z_score < 0)  # Clustered = closer than expected
  expect_true(result$n_profiles_used > 0L)
})

test_that("testProximity: empty profiles handled", {
  profiles <- tibble::tibble(
    gene_id = character(0), reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    detailed_events = list()
  )
  structures <- .make_spatial_structures()
  result <- testProximity(profiles, structures, n_perm = 10L)
  expect_true(is.na(result$z_score))
  expect_equal(result$n_profiles_used, 0L)
})

# ==========================================================================
# .extractEventsWithPositions()
# ==========================================================================

test_that(".extractEventsWithPositions: correct gene_start and gene_end", {
  profiles <- .make_spatial_profiles()[1, ]
  structures <- .make_spatial_structures()
  events <- Isopair:::.extractEventsWithPositions(profiles, structures)
  expect_true(all(events$gene_start == 1000L))
  expect_true(all(events$gene_end == 4000L))
})

# ==========================================================================
# summarizeBoundaryLengths()
# ==========================================================================

test_that("summarizeBoundaryLengths: only boundary types extracted", {
  # Profile with mix of boundary and non-boundary events
  profiles <- tibble::tibble(
    gene_id = "G1",
    reference_isoform_id = "R1",
    comparator_isoform_id = "C1",
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "A3SS", "SE", "Partial_IR_5", "Partial_IR_3", "IR"),
      direction = rep("LOSS", 6),
      five_prime = c(100L, 200L, 300L, 400L, 500L, 600L),
      three_prime = c(110L, 230L, 380L, 550L, 700L, 900L),
      bp_diff = c(10L, 30L, 80L, 150L, 200L, NA_integer_)
    ))
  )
  result <- summarizeBoundaryLengths(profiles)
  expect_equal(nrow(result$events), 4L)
  expect_true(all(result$events$event_type %in%
    c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")))
})

test_that("summarizeBoundaryLengths: summary_stats correct on known data", {
  profiles <- tibble::tibble(
    gene_id = c("G1", "G2"),
    reference_isoform_id = c("R1", "R2"),
    comparator_isoform_id = c("C1", "C2"),
    detailed_events = list(
      tibble::tibble(
        event_type = c("A5SS", "A5SS"),
        direction = c("LOSS", "LOSS"),
        five_prime = c(100L, 200L),
        three_prime = c(110L, 230L),
        bp_diff = c(10L, 30L)
      ),
      tibble::tibble(
        event_type = "A5SS",
        direction = "LOSS",
        five_prime = 300L,
        three_prime = 320L,
        bp_diff = 20L
      )
    )
  )
  result <- summarizeBoundaryLengths(profiles)
  a5ss <- result$summary_stats[result$summary_stats$event_type == "A5SS", ]
  expect_equal(a5ss$n, 3L)
  expect_equal(a5ss$median, 20)
  expect_equal(a5ss$mean, 20)
  expect_equal(a5ss$min, 10)
  expect_equal(a5ss$max, 30)
})

test_that("summarizeBoundaryLengths: combined_stats merges correctly", {
  profiles <- tibble::tibble(
    gene_id = "G1",
    reference_isoform_id = "R1",
    comparator_isoform_id = "C1",
    detailed_events = list(tibble::tibble(
      event_type = c("A5SS", "Partial_IR_5", "A3SS", "Partial_IR_3"),
      direction = rep("LOSS", 4),
      five_prime = c(100L, 200L, 300L, 400L),
      three_prime = c(110L, 350L, 330L, 650L),
      bp_diff = c(10L, 150L, 30L, 250L)
    ))
  )
  result <- summarizeBoundaryLengths(profiles)
  expect_true("5_prime" %in% result$combined_stats$event_type)
  expect_true("3_prime" %in% result$combined_stats$event_type)
  five_row <- result$combined_stats[result$combined_stats$event_type == "5_prime", ]
  expect_equal(five_row$n, 2L)  # A5SS + Partial_IR_5
})

# ==========================================================================
# classifyTopologyExpanded()
# ==========================================================================

test_that("classifyTopologyExpanded: returns list with expected names", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- classifyTopologyExpanded(profiles, structures,
                                      test = FALSE)
  expect_true(is.list(result))
  expect_true(all(c("per_type_pair", "collapsed", "comparison") %in%
    names(result)))
  expect_s3_class(result$comparison, "tbl_df")
})

test_that("classifyTopologyExpanded: per_type_pair sub-results correct", {
  # Profile with A5SS + A3SS → F2F (Profile 2)
  profiles <- .make_spatial_profiles()[2, ]
  structures <- .make_spatial_structures()
  result <- classifyTopologyExpanded(profiles, structures, test = FALSE)
  # Should have A5SS_x_A3SS in per_type_pair
  expect_true("A5SS_x_A3SS" %in% names(result$per_type_pair))
  a5a3 <- result$per_type_pair$A5SS_x_A3SS
  expect_equal(a5a3$classifications$topology[1], "F2F")
})

test_that("classifyTopologyExpanded: comparison F2F proportions match", {
  profiles <- .make_spatial_profiles()
  structures <- .make_spatial_structures()
  result <- classifyTopologyExpanded(profiles, structures, test = FALSE)
  # collapsed row should exist
  collapsed_row <- result$comparison[result$comparison$level == "collapsed", ]
  expect_equal(nrow(collapsed_row), 1L)
  # n_classified should be consistent
  if (!is.null(result$collapsed) && nrow(result$collapsed$summary) > 0L) {
    expect_equal(collapsed_row$n_classified,
                 sum(result$collapsed$summary$count))
  }
})
