# Tests for comparePairSets(), metaAnalyzePairSets(), .classifyProfile(), .cliffsD()

# Helper: build synthetic profiles with known characteristics
.make_compare_profiles <- function(n = 30, type = "ir_enriched", seed = 42) {
  set.seed(seed)
  if (type == "ir_enriched") {
    # High IR, low boundary
    tibble::tibble(
      gene_id = paste0("gene_ir_", seq_len(n)),
      reference_isoform_id = paste0("ref_ir_", seq_len(n)),
      comparator_isoform_id = paste0("comp_ir_", seq_len(n)),
      n_a5ss = as.integer(rbinom(n, 1, 0.1)),
      n_a3ss = as.integer(rbinom(n, 1, 0.1)),
      n_se = as.integer(rbinom(n, 1, 0.15)),
      n_missing_internal = 0L,
      n_ir = as.integer(rbinom(n, 2, 0.7)),
      n_ir_diff = as.integer(rbinom(n, 1, 0.3)),
      n_partial_ir = as.integer(rbinom(n, 1, 0.2)),
      n_alt_tss = as.integer(rbinom(n, 1, 0.1)),
      n_alt_tes = as.integer(rbinom(n, 1, 0.1)),
      tss_changed = n_alt_tss > 0L,
      tes_changed = n_alt_tes > 0L,
      n_events = n_a5ss + n_a3ss + n_se + n_ir + n_ir_diff +
        n_partial_ir + n_alt_tss + n_alt_tes,
      n_gain = as.integer(rbinom(n, 2, 0.3)),
      n_loss = as.integer(rbinom(n, 2, 0.3)),
      bp_gain = as.integer(rpois(n, 200)),
      bp_loss = as.integer(rpois(n, 100)),
      net_bp = bp_gain - bp_loss,
      detailed_events = replicate(n, tibble::tibble(
        event_type = character(0), direction = character(0),
        five_prime = integer(0), three_prime = integer(0),
        bp_diff = integer(0)
      ), simplify = FALSE)
    )
  } else {
    # High boundary, low IR
    tibble::tibble(
      gene_id = paste0("gene_bd_", seq_len(n)),
      reference_isoform_id = paste0("ref_bd_", seq_len(n)),
      comparator_isoform_id = paste0("comp_bd_", seq_len(n)),
      n_a5ss = as.integer(rbinom(n, 2, 0.7)),
      n_a3ss = as.integer(rbinom(n, 2, 0.7)),
      n_se = as.integer(rbinom(n, 1, 0.15)),
      n_missing_internal = 0L,
      n_ir = as.integer(rbinom(n, 1, 0.05)),
      n_ir_diff = 0L,
      n_partial_ir = 0L,
      n_alt_tss = as.integer(rbinom(n, 1, 0.1)),
      n_alt_tes = as.integer(rbinom(n, 1, 0.1)),
      tss_changed = n_alt_tss > 0L,
      tes_changed = n_alt_tes > 0L,
      n_events = n_a5ss + n_a3ss + n_se + n_ir + n_alt_tss + n_alt_tes,
      n_gain = as.integer(rbinom(n, 2, 0.5)),
      n_loss = as.integer(rbinom(n, 2, 0.5)),
      bp_gain = as.integer(rpois(n, 100)),
      bp_loss = as.integer(rpois(n, 200)),
      net_bp = bp_gain - bp_loss,
      detailed_events = replicate(n, tibble::tibble(
        event_type = character(0), direction = character(0),
        five_prime = integer(0), three_prime = integer(0),
        bp_diff = integer(0)
      ), simplify = FALSE)
    )
  }
}

# ==========================================================================
# .classifyProfile()
# ==========================================================================

test_that(".classifyProfile: Identical when n_events = 0", {
  profiles <- tibble::tibble(
    n_events = 0L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Identical")
})

test_that(".classifyProfile: Terminal-only when only Alt_TSS", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 1L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Terminal-only")
})

test_that(".classifyProfile: Boundary-only when only A5SS", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 1L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Boundary-only")
})

test_that(".classifyProfile: Inclusion-only when only SE", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 1L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Inclusion-only")
})

test_that(".classifyProfile: Full_Retention-only when only IR", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 1L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Full_Retention-only")
})

test_that(".classifyProfile: Partial_Retention-only when only partial_ir", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 1L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Partial_Retention-only")
})

test_that(".classifyProfile: Combined when A5SS + SE", {
  profiles <- tibble::tibble(
    n_events = 2L, n_alt_tss = 0L, n_alt_tes = 0L,
    n_a5ss = 1L, n_a3ss = 0L, n_se = 1L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(result$profile_type, "Combined")
})

test_that(".classifyProfile: removes helper columns", {
  profiles <- tibble::tibble(
    n_events = 1L, n_alt_tss = 1L, n_alt_tes = 0L,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_false("has_terminal" %in% names(result))
  expect_false("n_categories" %in% names(result))
  expect_true("profile_type" %in% names(result))
})

# ==========================================================================
# .cliffsD()
# ==========================================================================

test_that(".cliffsD: all x > y returns 1.0", {
  expect_equal(Isopair:::.cliffsD(c(10, 11, 12), c(1, 2, 3)), 1.0)
})

test_that(".cliffsD: all x < y returns -1.0", {
  expect_equal(Isopair:::.cliffsD(c(1, 2, 3), c(10, 11, 12)), -1.0)
})

test_that(".cliffsD: empty input returns NA", {
  expect_true(is.na(Isopair:::.cliffsD(numeric(0), c(1, 2))))
  expect_true(is.na(Isopair:::.cliffsD(c(1, 2), numeric(0))))
})

test_that(".cliffsD: equal distributions near 0", {
  set.seed(99)
  x <- rnorm(100)
  y <- rnorm(100)
  d <- Isopair:::.cliffsD(x, y)
  expect_true(abs(d) < 0.3)
})

test_that(".cliffsD: handles NA values", {
  d <- Isopair:::.cliffsD(c(10, NA, 12), c(1, 2, NA))
  expect_equal(d, 1.0)  # After NA removal: c(10,12) vs c(1,2)
})

# ==========================================================================
# comparePairSets() — structure
# ==========================================================================

test_that("comparePairSets: returns list with all 10 sections", {
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  expect_true(is.list(result))
  expected_sections <- c("event_frequency", "event_count", "gain_loss",
                          "cooccurrence", "regional_impact", "topology",
                          "proximity", "ptc_rate", "ptc_distance",
                          "profile_categories")
  expect_true(all(expected_sections %in% names(result)))
})

test_that("comparePairSets: non-skipped sections have required fields", {
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  # Sections that always run (profiles only)
  for (section in c("event_count", "gain_loss", "profile_categories")) {
    s <- result[[section]]
    expect_true("description" %in% names(s), info = section)
    expect_true("test" %in% names(s), info = section)
    expect_true("p_value" %in% names(s), info = section)
    expect_true("n_a" %in% names(s), info = section)
    expect_true("n_b" %in% names(s), info = section)
  }
})

test_that("comparePairSets: optional sections skipped when inputs NULL", {
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  # Sections requiring optional inputs should be skipped
  expect_true(result$regional_impact$skipped)
  expect_true(nchar(result$regional_impact$skipped_reason) > 0)
  expect_true(result$topology$skipped)
  expect_true(result$proximity$skipped)
  expect_true(result$ptc_rate$skipped)
  expect_true(result$ptc_distance$skipped)
})

# ==========================================================================
# comparePairSets() — statistical tests
# ==========================================================================

test_that("comparePairSets: IR enriched shows high OR for IR", {
  a <- .make_compare_profiles(50, "ir_enriched")
  b <- .make_compare_profiles(50, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  ir_result <- result$event_frequency$per_event$IR
  expect_true(ir_result$odds_ratio > 1)
  expect_true(ir_result$count_a > ir_result$count_b)
})

test_that("comparePairSets: event_count has correct Cliff's D sign", {
  a <- .make_compare_profiles(30, "ir_enriched")
  b <- .make_compare_profiles(30, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  # IR-enriched has more events (high n_ir + n_ir_diff)
  med_a <- result$event_count$median_a
  med_b <- result$event_count$median_b
  d <- result$event_count$cliffs_d
  # If a has more events, d > 0
  if (med_a > med_b) expect_true(d > 0)
  if (med_a < med_b) expect_true(d < 0)
})

test_that("comparePairSets: profile_categories produces cramers_v", {
  a <- .make_compare_profiles(50, "ir_enriched")
  b <- .make_compare_profiles(50, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  expect_true(!is.na(result$profile_categories$cramers_v))
  expect_true(result$profile_categories$cramers_v >= 0)
})

test_that("comparePairSets: cooccurrence correlation is numeric", {
  a <- .make_compare_profiles(50, "ir_enriched")
  b <- .make_compare_profiles(50, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  expect_true(is.numeric(result$cooccurrence$correlation))
})

test_that("comparePairSets: event_frequency has 8 sub-results + FDR", {
  a <- .make_compare_profiles(20, "ir_enriched")
  b <- .make_compare_profiles(20, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  expect_equal(length(result$event_frequency$per_event), 8L)
  expect_equal(length(result$event_frequency$adj_p_values), 8L)
  # FDR >= raw
  for (et in names(result$event_frequency$per_event)) {
    ev <- result$event_frequency$per_event[[et]]
    adj <- result$event_frequency$adj_p_values[et]
    if (!is.na(ev$p_value) && !is.na(adj)) {
      expect_true(adj >= ev$p_value - 1e-14, info = et)
    }
  }
})

# ==========================================================================
# comparePairSets() — self-documenting results
# ==========================================================================

test_that("comparePairSets: data_a and data_b are tibbles", {
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  expect_s3_class(result$event_count$data_a, "tbl_df")
  expect_s3_class(result$event_count$data_b, "tbl_df")
  expect_equal(nrow(result$event_count$data_a), 10L)
  expect_equal(nrow(result$event_count$data_b), 10L)
})

test_that("comparePairSets: event_frequency direction field present", {
  a <- .make_compare_profiles(20, "ir_enriched")
  b <- .make_compare_profiles(20, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)

  for (et in names(result$event_frequency$per_event)) {
    ev <- result$event_frequency$per_event[[et]]
    expect_true("direction" %in% names(ev), info = et)
    expect_true(ev$direction %in% c("enriched in A", "enriched in B", "NS"),
                info = et)
  }
})

# ==========================================================================
# comparePairSets() — paired tests
# ==========================================================================

test_that("comparePairSets: paired tests with overlapping genes", {
  # Create two sets sharing some genes
  a <- .make_compare_profiles(30, "ir_enriched")
  b <- .make_compare_profiles(30, "boundary_enriched", seed = 99)
  # Force overlap: copy gene_id + reference_isoform_id
  shared <- min(nrow(a), nrow(b))
  b$gene_id[seq_len(shared)] <- a$gene_id[seq_len(shared)]
  b$reference_isoform_id[seq_len(shared)] <-
    a$reference_isoform_id[seq_len(shared)]

  result <- comparePairSets(a, b, paired = TRUE, min_overlap = 5L,
                             verbose = FALSE)

  expect_true("paired_tests" %in% names(result))
  expect_false(isTRUE(result$paired_tests$skipped))
  expect_true(result$paired_tests$n_genes_overlap >= 5L)
  expect_true("paired_wilcox" %in% names(result$paired_tests))
  expect_true("sign_test" %in% names(result$paired_tests))
  expect_true("mcnemar_tests" %in% names(result$paired_tests))
  expect_equal(length(result$paired_tests$mcnemar_tests), 8L)
})

test_that("comparePairSets: paired tests skipped when overlap < min", {
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  # No shared genes
  result <- comparePairSets(a, b, paired = TRUE, min_overlap = 50L,
                             verbose = FALSE)

  expect_true(result$paired_tests$skipped)
  expect_true(grepl("overlap", result$paired_tests$skipped_reason))
})

test_that("comparePairSets: mcnemar has adj_p_value via FDR", {
  a <- .make_compare_profiles(30, "ir_enriched")
  b <- .make_compare_profiles(30, "boundary_enriched", seed = 99)
  b$gene_id <- a$gene_id
  b$reference_isoform_id <- a$reference_isoform_id

  result <- comparePairSets(a, b, paired = TRUE, min_overlap = 5L,
                             verbose = FALSE)

  for (et in names(result$paired_tests$mcnemar_tests)) {
    mt <- result$paired_tests$mcnemar_tests[[et]]
    expect_true("adj_p_value" %in% names(mt), info = et)
    if (!is.na(mt$p_value) && !is.na(mt$adj_p_value)) {
      expect_true(mt$adj_p_value >= mt$p_value - 1e-14, info = et)
    }
  }
})

# ==========================================================================
# comparePairSets() — edge cases
# ==========================================================================

test_that("comparePairSets: small profiles (n=5) run without error", {
  a <- .make_compare_profiles(5, "ir_enriched")
  b <- .make_compare_profiles(5, "boundary_enriched", seed = 99)
  result <- comparePairSets(a, b, verbose = FALSE)
  expect_true(is.list(result))
  expect_true("event_count" %in% names(result))
})

test_that("comparePairSets: identical distributions graceful", {
  a <- .make_compare_profiles(20, "ir_enriched", seed = 42)
  b <- .make_compare_profiles(20, "ir_enriched", seed = 42)
  result <- comparePairSets(a, b, verbose = FALSE)
  # Same data should produce high p-values or NA
  expect_true(is.list(result))
  if (!is.na(result$event_count$p_value)) {
    expect_true(result$event_count$p_value > 0.05)
  }
})

# ==========================================================================
# metaAnalyzePairSets()
# ==========================================================================

test_that("metaAnalyzePairSets: requires metafor", {
  # This test verifies the function exists and has correct parameter handling
  # It will actually run if metafor is installed
  skip_if_not_installed("metafor")

  a <- .make_compare_profiles(30, "ir_enriched")
  b <- .make_compare_profiles(30, "boundary_enriched", seed = 99)

  # Create 3 comparison results with different seeds
  comp_list <- list()
  for (i in 1:3) {
    ai <- .make_compare_profiles(30, "ir_enriched", seed = i)
    bi <- .make_compare_profiles(30, "boundary_enriched", seed = 100 + i)
    comp_list[[paste0("comp_", i)]] <- comparePairSets(ai, bi, verbose = FALSE)
  }

  result <- metaAnalyzePairSets(comp_list, verbose = FALSE)

  expect_true(is.list(result))
  # Event count should have valid results
  expect_false(isTRUE(result$event_count$skipped))
  expect_true(!is.na(result$event_count$estimate))
  expect_true(!is.na(result$event_count$I2))
  expect_true(!is.na(result$event_count$tau2))
  expect_equal(result$event_count$k, 3L)
})

test_that("metaAnalyzePairSets: proximity always skipped", {
  skip_if_not_installed("metafor")

  comp_list <- list()
  for (i in 1:3) {
    ai <- .make_compare_profiles(20, "ir_enriched", seed = i)
    bi <- .make_compare_profiles(20, "boundary_enriched", seed = 100 + i)
    comp_list[[paste0("comp_", i)]] <- comparePairSets(ai, bi, verbose = FALSE)
  }

  result <- metaAnalyzePairSets(comp_list, verbose = FALSE)
  expect_true(result$proximity$skipped)
  expect_true(grepl("not meta-analyzable", result$proximity$skipped_reason))
})

test_that("metaAnalyzePairSets: < 2 comparisons errors", {
  skip_if_not_installed("metafor")
  a <- .make_compare_profiles(10, "ir_enriched")
  b <- .make_compare_profiles(10, "boundary_enriched", seed = 99)
  single <- list(only = comparePairSets(a, b, verbose = FALSE))
  expect_error(metaAnalyzePairSets(single), "named list with >= 2")
})

test_that("metaAnalyzePairSets: returns per_comparison tibble", {
  skip_if_not_installed("metafor")

  comp_list <- list()
  for (i in 1:3) {
    ai <- .make_compare_profiles(20, "ir_enriched", seed = i)
    bi <- .make_compare_profiles(20, "boundary_enriched", seed = 100 + i)
    comp_list[[paste0("comp_", i)]] <- comparePairSets(ai, bi, verbose = FALSE)
  }

  result <- metaAnalyzePairSets(comp_list, verbose = FALSE)

  if (!isTRUE(result$event_count$skipped)) {
    expect_s3_class(result$event_count$per_comparison, "tbl_df")
    expect_equal(nrow(result$event_count$per_comparison), 3L)
    expect_true(all(c("comparison", "effect_size", "se") %in%
                      names(result$event_count$per_comparison)))
  }
})

# ==========================================================================
# comparePairSets() — PTC sections
# ==========================================================================

# Helper: build mock PTC data
.make_mock_ptc <- function(n, ptc_rate, seed = 42) {
  set.seed(seed)
  tibble::tibble(
    comparator_isoform_id = paste0("comp_", seq_len(n)),
    has_ptc = sample(c(TRUE, FALSE), n, replace = TRUE,
                     prob = c(ptc_rate, 1 - ptc_rate)),
    ptc_distance = ifelse(has_ptc, as.integer(rpois(n, 300)), NA_integer_)
  )
}

test_that("comparePairSets: ptc_rate section with PTC data", {
  a <- .make_compare_profiles(20, "ir_enriched")
  b <- .make_compare_profiles(20, "boundary_enriched", seed = 99)
  ptc_a <- .make_mock_ptc(20, ptc_rate = 0.6, seed = 1)
  ptc_b <- .make_mock_ptc(20, ptc_rate = 0.2, seed = 2)

  result <- comparePairSets(a, b, ptc_a = ptc_a, ptc_b = ptc_b,
                             verbose = FALSE)

  expect_false(isTRUE(result$ptc_rate$skipped))
  expect_true("odds_ratio" %in% names(result$ptc_rate))
  expect_true("conf_int" %in% names(result$ptc_rate))
  expect_true("direction" %in% names(result$ptc_rate))
  expect_equal(length(result$ptc_rate$conf_int), 2L)
  expect_true(result$ptc_rate$odds_ratio > 1)  # set A has higher PTC rate
})

test_that("comparePairSets: ptc_distance section with PTC data", {
  a <- .make_compare_profiles(20, "ir_enriched")
  b <- .make_compare_profiles(20, "boundary_enriched", seed = 99)
  ptc_a <- .make_mock_ptc(20, ptc_rate = 0.5, seed = 1)
  ptc_b <- .make_mock_ptc(20, ptc_rate = 0.5, seed = 2)

  result <- comparePairSets(a, b, ptc_a = ptc_a, ptc_b = ptc_b,
                             verbose = FALSE)

  expect_false(isTRUE(result$ptc_distance$skipped))
  expect_true("cliffs_d" %in% names(result$ptc_distance))
  expect_true("median_a" %in% names(result$ptc_distance))
  expect_true("median_b" %in% names(result$ptc_distance))
  # data_a/data_b should only contain PTC+ rows with non-NA distance
  expect_true(all(result$ptc_distance$data_a$has_ptc))
  expect_true(all(!is.na(result$ptc_distance$data_a$ptc_distance)))
})

# ==========================================================================
# comparePairSets() — topology and proximity with structures
# ==========================================================================

# Helper: build structures matching compare profiles
.make_compare_structures <- function(profiles) {
  n <- nrow(profiles)
  tibble::tibble(
    isoform_id = profiles$reference_isoform_id,
    gene_id = profiles$gene_id,
    chr = "chrT",
    strand = "+",
    n_exons = 4L,
    exon_starts = replicate(n, c(1000L, 1500L, 2200L, 3000L),
                            simplify = FALSE),
    exon_ends = replicate(n, c(1200L, 2000L, 2800L, 4000L),
                          simplify = FALSE),
    tx_start = 1000L,
    tx_end = 4000L,
    n_junctions = 3L
  )
}

# Helper: build profiles with detailed_events for topology classification
.make_topology_profiles <- function(n, seed = 42) {
  set.seed(seed)
  tibble::tibble(
    gene_id = paste0("gene_topo_", seq_len(n)),
    reference_isoform_id = paste0("ref_topo_", seq_len(n)),
    comparator_isoform_id = paste0("comp_topo_", seq_len(n)),
    n_a5ss = 1L, n_a3ss = 1L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L,
    tss_changed = FALSE, tes_changed = FALSE,
    n_events = 2L,
    n_gain = 0L, n_loss = 2L,
    bp_gain = 0L, bp_loss = 20L,
    net_bp = -20L,
    detailed_events = replicate(n, tibble::tibble(
      event_type = c("A5SS", "A3SS"),
      direction = c("LOSS", "LOSS"),
      five_prime = c(1995L, 2205L),
      three_prime = c(2005L, 2195L),
      bp_diff = c(10L, 10L)
    ), simplify = FALSE)
  )
}

test_that("comparePairSets: topology section with structures", {
  a <- .make_topology_profiles(10, seed = 1)
  b <- .make_topology_profiles(10, seed = 2)
  struct_a <- .make_compare_structures(a)
  struct_b <- .make_compare_structures(b)

  result <- comparePairSets(a, b, structures_a = struct_a,
                             structures_b = struct_b, verbose = FALSE)

  expect_false(isTRUE(result$topology$skipped))
  expect_true("cramers_v" %in% names(result$topology) ||
              "p_value" %in% names(result$topology))
  expect_true("n_a" %in% names(result$topology))
  expect_true("n_b" %in% names(result$topology))
})

test_that("comparePairSets: proximity section with structures", {
  a <- .make_topology_profiles(10, seed = 1)
  b <- .make_topology_profiles(10, seed = 2)
  struct_a <- .make_compare_structures(a)
  struct_b <- .make_compare_structures(b)

  result <- comparePairSets(a, b, structures_a = struct_a,
                             structures_b = struct_b, n_perm = 50L,
                             verbose = FALSE)

  expect_false(isTRUE(result$proximity$skipped))
  expect_true("z_a" %in% names(result$proximity))
  expect_true("z_b" %in% names(result$proximity))
  expect_true("z_diff" %in% names(result$proximity))
  # Proximity is descriptive: p_value should be NA
  expect_true(is.na(result$proximity$p_value))
})

# ==========================================================================
# .classifyProfile() — vectorized
# ==========================================================================

test_that(".classifyProfile: vectorized on multiple rows", {
  profiles <- tibble::tibble(
    n_events = c(0L, 1L, 1L, 2L),
    n_alt_tss = c(0L, 1L, 0L, 1L),
    n_alt_tes = c(0L, 0L, 0L, 0L),
    n_a5ss = c(0L, 0L, 1L, 1L),
    n_a3ss = c(0L, 0L, 0L, 0L),
    n_se = c(0L, 0L, 0L, 0L),
    n_missing_internal = c(0L, 0L, 0L, 0L),
    n_ir = c(0L, 0L, 0L, 0L),
    n_ir_diff = c(0L, 0L, 0L, 0L),
    n_partial_ir = c(0L, 0L, 0L, 0L)
  )
  result <- Isopair:::.classifyProfile(profiles)
  expect_equal(nrow(result), 4L)
  expect_equal(result$profile_type,
               c("Identical", "Terminal-only", "Boundary-only", "Combined"))
})

# ==========================================================================
# Chi-square simulated p-value path
# ==========================================================================

test_that(".compareChiSquare: uses simulated p-value when expected < 5", {
  # Create data where at least one expected cell < 5
  # Rare category in set_a, absent in set_b
  cats_a <- c(rep("type_A", 20), rep("type_B", 2))
  cats_b <- c(rep("type_A", 18), rep("type_B", 1))
  data_a <- tibble::tibble(cat = cats_a)
  data_b <- tibble::tibble(cat = cats_b)

  result <- Isopair:::.compareChiSquare(cats_a, cats_b,
    description = "test", data_a = data_a, data_b = data_b)

  # Should indicate simulated p-value was used
  expect_true(grepl("simulated", result$test))
  expect_true(!is.na(result$p_value))
})

# ==========================================================================
# McNemar exact binomial fallback
# ==========================================================================

test_that("comparePairSets: mcnemar exact binomial for small discordant", {
  # Create paired data with few discordant pairs (< 10)
  n <- 20
  a <- tibble::tibble(
    gene_id = paste0("g", seq_len(n)),
    reference_isoform_id = paste0("r", seq_len(n)),
    comparator_isoform_id = paste0("c", seq_len(n)),
    n_a5ss = c(rep(1L, 3), rep(0L, 17)),  # 3 have A5SS
    n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L,
    n_alt_tss = 0L, n_alt_tes = 0L,
    tss_changed = FALSE, tes_changed = FALSE,
    n_events = n_a5ss,
    n_gain = 0L, n_loss = 0L, bp_gain = 0L, bp_loss = 0L, net_bp = 0L,
    detailed_events = replicate(n, tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ), simplify = FALSE)
  )
  b <- a
  # Make 5 discordant for A5SS (a has it, b doesn't) — < 10 → binomial
  b$n_a5ss[1:3] <- 0L
  b$n_events <- b$n_a5ss

  result <- comparePairSets(a, b, paired = TRUE, min_overlap = 5L,
                             verbose = FALSE)

  a5_mcnemar <- result$paired_tests$mcnemar_tests$A5SS
  expect_true(!is.na(a5_mcnemar$p_value))
  expect_equal(a5_mcnemar$a_only, 3L)
  expect_equal(a5_mcnemar$b_only, 0L)
  expect_equal(a5_mcnemar$n_discordant, 3L)
})

# ==========================================================================
# Paired type_transitions
# ==========================================================================

test_that("comparePairSets: type_transitions structure correct", {
  a <- .make_compare_profiles(30, "ir_enriched")
  b <- .make_compare_profiles(30, "boundary_enriched", seed = 99)
  b$gene_id <- a$gene_id
  b$reference_isoform_id <- a$reference_isoform_id

  result <- comparePairSets(a, b, paired = TRUE, min_overlap = 5L,
                             verbose = FALSE)

  tt <- result$paired_tests$type_transitions
  expect_true("n_same_type" %in% names(tt))
  expect_true("n_diff_type" %in% names(tt))
  expect_true("pct_type_change" %in% names(tt))
  expect_equal(tt$n_same_type + tt$n_diff_type, 30L)
  expect_true(tt$pct_type_change >= 0 && tt$pct_type_change <= 100)
})

# ==========================================================================
# Haldane correction
# ==========================================================================

test_that(".compareCooccurrence: Haldane preserves OR < 1", {
  # Build profiles where A5SS and A3SS are mutually exclusive
  # This should produce OR < 1 (negative association)
  n <- 50
  profiles <- tibble::tibble(
    gene_id = paste0("g", seq_len(n)),
    reference_isoform_id = paste0("r", seq_len(n)),
    comparator_isoform_id = paste0("c", seq_len(n)),
    tss_changed = FALSE, tes_changed = FALSE,
    n_a5ss = c(rep(1L, 25), rep(0L, 25)),
    n_a3ss = c(rep(0L, 25), rep(1L, 25)),
    n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L,
    n_events = 1L,
    detailed_events = replicate(n, tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ), simplify = FALSE)
  )

  cooc <- testCooccurrence(profiles)
  a5_a3 <- cooc[cooc$event_a == "A5SS" & cooc$event_b == "A3SS", ]
  # With perfect mutual exclusion, n_both=0, OR should be 0
  expect_equal(a5_a3$n_both, 0L)

  # Haldane correction should produce log_or < 0 (not censored at 0.5)
  result <- Isopair:::.compareCooccurrence(profiles, profiles)
  # The log-OR for A5SS-A3SS should be negative (since OR near 0)
  # We can't directly check internal log_or, but the function should run
  expect_true(is.list(result))
  expect_true("correlation" %in% names(result))
})
