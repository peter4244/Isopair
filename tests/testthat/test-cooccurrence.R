# Tests for testCooccurrence()

# Helper: build minimal profiles with known event counts
.make_cooccurrence_profiles <- function(n = 100) {
  set.seed(42)
  # Create profiles with known co-occurrence patterns
  tibble::tibble(
    gene_id = paste0("gene_", seq_len(n)),
    reference_isoform_id = paste0("ref_", seq_len(n)),
    comparator_isoform_id = paste0("comp_", seq_len(n)),
    tss_changed = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.3, 0.7)),
    tes_changed = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.2, 0.8)),
    n_a5ss = as.integer(rbinom(n, 3, 0.3)),
    n_a3ss = as.integer(rbinom(n, 3, 0.3)),
    n_se = as.integer(rbinom(n, 2, 0.2)),
    n_missing_internal = as.integer(rbinom(n, 1, 0.1)),
    n_ir = as.integer(rbinom(n, 2, 0.15)),
    n_ir_diff = as.integer(rbinom(n, 1, 0.1)),
    n_partial_ir = as.integer(rbinom(n, 1, 0.1)),
    n_events = 1L,
    detailed_events = replicate(n, tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ), simplify = FALSE)
  )
}

test_that("testCooccurrence returns correct structure", {
  profiles <- .make_cooccurrence_profiles(50)
  result <- testCooccurrence(profiles)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("event_a", "event_b", "odds_ratio", "p_value",
                     "adj_p_value", "conf_low", "conf_high",
                     "n_both", "n_a_only", "n_b_only", "n_neither",
                     "significant") %in% names(result)))
})

test_that("testCooccurrence produces 28 pairwise tests", {
  profiles <- .make_cooccurrence_profiles(50)
  result <- testCooccurrence(profiles)
  expect_equal(nrow(result), 28L)  # C(8,2)
})

test_that("testCooccurrence: contingency table sums to n", {
  profiles <- .make_cooccurrence_profiles(100)
  result <- testCooccurrence(profiles)
  totals <- result$n_both + result$n_a_only + result$n_b_only + result$n_neither
  expect_true(all(totals == 100L))
})

test_that("testCooccurrence: odds_ratio >= 0 when computable", {
  profiles <- .make_cooccurrence_profiles(200)
  result <- testCooccurrence(profiles)
  non_na <- result$odds_ratio[!is.na(result$odds_ratio)]
  expect_true(all(non_na >= 0))
})

test_that("testCooccurrence: perfect co-occurrence detected", {
  # A5SS and A3SS always co-occur
  profiles <- tibble::tibble(
    gene_id = paste0("g", 1:50),
    reference_isoform_id = paste0("r", 1:50),
    comparator_isoform_id = paste0("c", 1:50),
    tss_changed = FALSE, tes_changed = FALSE,
    n_a5ss = rep(1L, 50), n_a3ss = rep(1L, 50),
    n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L,
    n_events = 1L,
    detailed_events = replicate(50, tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ), simplify = FALSE)
  )
  result <- testCooccurrence(profiles)
  a5_a3 <- result[result$event_a == "A5SS" & result$event_b == "A3SS", ]
  expect_equal(a5_a3$n_both, 50L)
  expect_equal(a5_a3$n_a_only, 0L)
  expect_equal(a5_a3$n_b_only, 0L)
})

test_that("testCooccurrence: mutual exclusion detected", {
  # A5SS and SE never co-occur
  n <- 100
  profiles <- tibble::tibble(
    gene_id = paste0("g", 1:n),
    reference_isoform_id = paste0("r", 1:n),
    comparator_isoform_id = paste0("c", 1:n),
    tss_changed = FALSE, tes_changed = FALSE,
    n_a5ss = c(rep(1L, 50), rep(0L, 50)),
    n_a3ss = 0L,
    n_se = c(rep(0L, 50), rep(1L, 50)),
    n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 0L, n_partial_ir = 0L,
    n_events = 1L,
    detailed_events = replicate(n, tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    ), simplify = FALSE)
  )
  result <- testCooccurrence(profiles)
  a5_se <- result[result$event_a == "A5SS" & result$event_b == "SE", ]
  expect_equal(a5_se$n_both, 0L)
})

test_that("testCooccurrence: IR and IR_diff merged", {
  profiles <- .make_cooccurrence_profiles(50)
  result <- testCooccurrence(profiles)
  expect_true("IR" %in% c(result$event_a, result$event_b))
  expect_false("IR_diff" %in% c(result$event_a, result$event_b))
})

test_that("testCooccurrence: FDR correction applied", {
  profiles <- .make_cooccurrence_profiles(200)
  result <- testCooccurrence(profiles)
  finite <- is.finite(result$adj_p_value) & is.finite(result$p_value)
  # Allow machine epsilon tolerance for floating point
  expect_true(all(result$adj_p_value[finite] >= result$p_value[finite] - 1e-14))
})

test_that("testCooccurrence: stratified analysis works", {
  profiles <- .make_cooccurrence_profiles(100)
  profiles$group <- rep(c("A", "B"), each = 50)
  result <- testCooccurrence(profiles, stratify_by = "group")
  expect_true("stratum" %in% names(result))
  expect_equal(length(unique(result$stratum)), 2L)
  expect_equal(nrow(result), 56L)  # 28 tests * 2 strata
})

test_that("testCooccurrence: error on invalid stratify_by", {
  profiles <- .make_cooccurrence_profiles(10)
  expect_error(testCooccurrence(profiles, stratify_by = "nonexistent"),
               "not found")
})

test_that(".addEventBinary: IR merges ir and ir_diff", {
  profiles <- tibble::tibble(
    tss_changed = FALSE, tes_changed = FALSE,
    n_a5ss = 0L, n_a3ss = 0L, n_se = 0L, n_missing_internal = 0L,
    n_ir = 0L, n_ir_diff = 1L, n_partial_ir = 0L
  )
  result <- Isopair:::.addEventBinary(profiles)
  expect_equal(result$has_ir, 1L)
})
