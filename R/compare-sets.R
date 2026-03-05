#' Compare Two Sets of Isoform Pair Profiles
#'
#' Runs a battery of 10 statistical tests comparing two sets of splicing
#' choice profiles. Each test section is self-documenting: the result
#' includes the test description, method name, p-value, effect size, and
#' the exact data that entered the test. Optional inputs (structures,
#' event_regions, ptc) enable additional analyses; sections requiring
#' missing inputs are gracefully skipped.
#'
#' @param profiles_a,profiles_b Tibbles from \code{\link{buildProfiles}()}.
#' @param structures_a,structures_b Optional tibbles from
#'   \code{\link{parseIsoformStructures}()}. Required for topology and
#'   proximity comparisons.
#' @param event_regions_a,event_regions_b Optional tibbles from
#'   \code{\link{mapEventsToRegions}()}. Required for regional impact.
#' @param ptc_a,ptc_b Optional tibbles from
#'   \code{\link{computePtcStatus}()}. Required for PTC comparisons.
#' @param paired Logical; if TRUE, run gene-matched paired tests on genes
#'   present in both sets (matched on gene_id + reference_isoform_id).
#' @param n_perm Integer; number of permutations for proximity tests.
#'   Default 1000.
#' @param min_overlap Integer; minimum gene overlap required for paired
#'   tests. Default 50.
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A named list with 10 test sections (event_frequency,
#'   event_count, gain_loss, cooccurrence, regional_impact, topology,
#'   proximity, ptc_rate, ptc_distance, profile_categories) and optionally
#'   \code{paired_tests}. Each section is a self-documenting list with
#'   \code{description}, \code{test}, \code{p_value}, effect size,
#'   \code{n_a}, \code{n_b}, and attached data. Skipped sections have
#'   \code{skipped = TRUE} and \code{skipped_reason}.
#' @examples
#' data(example_profiles)
#' result <- comparePairSets(example_profiles, example_profiles,
#'   min_overlap = 2L, verbose = FALSE)
#' names(result)
#' @export
#' @importFrom stats wilcox.test chisq.test cor.test binom.test
#'   fisher.test mcnemar.test p.adjust
#' @importFrom dplyr bind_rows case_when inner_join mutate
#' @importFrom tibble tibble
comparePairSets <- function(profiles_a, profiles_b,
                             structures_a = NULL, structures_b = NULL,
                             event_regions_a = NULL, event_regions_b = NULL,
                             ptc_a = NULL, ptc_b = NULL,
                             paired = FALSE, n_perm = 1000L,
                             min_overlap = 50L, verbose = TRUE) {

  n_a <- nrow(profiles_a)
  n_b <- nrow(profiles_b)
  if (verbose) message(sprintf("Comparing %d (set A) vs %d (set B) profiles", n_a, n_b))

  results <- list()

  # --- Section 1: Event frequency (Fisher's exact x 8 + FDR) ---
  results$event_frequency <- .compareEventFrequency(profiles_a, profiles_b)

  # --- Section 2: Event count (Wilcoxon + Cliff's D) ---
  results$event_count <- .compareDistribution(
    profiles_a$n_events, profiles_b$n_events,
    description = "Event complexity: total events per pair",
    n_a = n_a, n_b = n_b,
    data_a = profiles_a, data_b = profiles_b
  )

  # --- Section 3: Gain/loss balance ---
  results$gain_loss <- .compareDistribution(
    profiles_a$net_bp, profiles_b$net_bp,
    description = "Net basepair change: set A vs set B",
    n_a = n_a, n_b = n_b,
    data_a = profiles_a, data_b = profiles_b
  )

  # --- Section 4: Co-occurrence correlation ---
  results$cooccurrence <- .compareCooccurrence(profiles_a, profiles_b)

  # --- Section 5: Regional impact ---
  if (!is.null(event_regions_a) && !is.null(event_regions_b)) {
    results$regional_impact <- .compareChiSquare(
      event_regions_a$region_type, event_regions_b$region_type,
      description = "Event distribution across genomic regions",
      data_a = event_regions_a, data_b = event_regions_b
    )
  } else {
    results$regional_impact <- list(
      skipped = TRUE,
      skipped_reason = "event_regions_a and/or event_regions_b not provided"
    )
  }

  # --- Section 6: Topology ---
  if (!is.null(structures_a) && !is.null(structures_b)) {
    results$topology <- .compareTopology(
      profiles_a, profiles_b, structures_a, structures_b
    )
  } else {
    results$topology <- list(
      skipped = TRUE,
      skipped_reason = "structures_a and/or structures_b not provided"
    )
  }

  # --- Section 7: Proximity (descriptive) ---
  if (!is.null(structures_a) && !is.null(structures_b)) {
    results$proximity <- .compareProximity(
      profiles_a, profiles_b, structures_a, structures_b, n_perm
    )
  } else {
    results$proximity <- list(
      skipped = TRUE,
      skipped_reason = "structures_a and/or structures_b not provided"
    )
  }

  # --- Section 8: PTC rate ---
  if (!is.null(ptc_a) && !is.null(ptc_b)) {
    results$ptc_rate <- .comparePtcRate(ptc_a, ptc_b)
  } else {
    results$ptc_rate <- list(
      skipped = TRUE,
      skipped_reason = "ptc_a and/or ptc_b not provided"
    )
  }

  # --- Section 9: PTC distance ---
  if (!is.null(ptc_a) && !is.null(ptc_b)) {
    results$ptc_distance <- .comparePtcDistance(ptc_a, ptc_b)
  } else {
    results$ptc_distance <- list(
      skipped = TRUE,
      skipped_reason = "ptc_a and/or ptc_b not provided"
    )
  }

  # --- Section 10: Profile categories ---
  results$profile_categories <- .compareProfileCategories(
    profiles_a, profiles_b
  )

  # --- Paired tests ---
  if (paired) {
    results$paired_tests <- .runPairedTests(
      profiles_a, profiles_b, min_overlap, verbose
    )
  }

  results
}


#' Meta-Analyze Multiple Pair-Set Comparisons
#'
#' Pools effect sizes from multiple \code{\link{comparePairSets}()} results
#' using random-effects meta-analysis via the \pkg{metafor} package.
#' Each test section is meta-analyzed independently. Proximity is skipped
#' (no reliable SE available).
#'
#' @param comparison_list A named list of \code{\link{comparePairSets}()}
#'   results. Names are comparison labels (e.g., cell type names).
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A named list parallel to the \code{comparePairSets()} output.
#'   Each element contains: \code{estimate}, \code{se}, \code{p_value},
#'   \code{ci_lower}, \code{ci_upper}, \code{Q}, \code{I2}, \code{tau2},
#'   \code{k}, \code{per_comparison}. Tests with < 2 valid comparisons
#'   are skipped with \code{skipped_reason}. Note: Cramér's V
#'   meta-analysis is interpretive (V is non-negative; pooled estimates
#'   indicate effect magnitude, not direction).
#' @examples
#' \donttest{
#' if (requireNamespace("metafor", quietly = TRUE)) {
#'   data(example_profiles)
#'   comp1 <- comparePairSets(example_profiles, example_profiles,
#'     verbose = FALSE)
#'   comp2 <- comparePairSets(example_profiles, example_profiles,
#'     verbose = FALSE)
#'   meta <- metaAnalyzePairSets(list(set1 = comp1, set2 = comp2),
#'     verbose = FALSE)
#'   names(meta)
#' }
#' }
#' @export
#' @importFrom tibble tibble
metaAnalyzePairSets <- function(comparison_list, verbose = TRUE) {

  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for meta-analysis. ",
         "Install with: install.packages('metafor')",
         call. = FALSE)
  }

  comp_names <- names(comparison_list)
  if (is.null(comp_names) || length(comp_names) < 2L) {
    stop("comparison_list must be a named list with >= 2 elements",
         call. = FALSE)
  }

  results <- list()

  # --- Event count ---
  results$event_count <- .metaCliffD(
    comparison_list, comp_names, "event_count", verbose
  )

  # --- Gain/loss ---
  results$gain_loss <- .metaCliffD(
    comparison_list, comp_names, "gain_loss", verbose
  )

  # --- Event frequency (per event type) ---
  results$event_frequency <- .metaEventFrequency(
    comparison_list, comp_names, verbose
  )

  # --- Co-occurrence ---
  results$cooccurrence <- .metaCorrelation(
    comparison_list, comp_names, verbose
  )

  # --- Regional impact ---
  results$regional_impact <- .metaCramersV(
    comparison_list, comp_names, "regional_impact", verbose
  )

  # --- Topology ---
  results$topology <- .metaCramersV(
    comparison_list, comp_names, "topology", verbose
  )

  # --- Proximity (not meta-analyzable) ---
  results$proximity <- list(
    skipped = TRUE,
    skipped_reason = "Proximity z-scores do not have reliable SE; not meta-analyzable"
  )

  # --- PTC rate ---
  results$ptc_rate <- .metaLogOR(
    comparison_list, comp_names, "ptc_rate", verbose
  )

  # --- PTC distance ---
  results$ptc_distance <- .metaCliffD(
    comparison_list, comp_names, "ptc_distance", verbose
  )

  # --- Profile categories ---
  results$profile_categories <- .metaCramersV(
    comparison_list, comp_names, "profile_categories", verbose
  )

  results
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Classify profiles into pattern categories
#'
#' Adds a \code{profile_type} column categorizing each profile into one
#' of 7 types based on which event categories are present.
#'
#' @param profiles A tibble from \code{buildProfiles()}.
#' @return The input tibble with added \code{profile_type} column.
#' @keywords internal
.classifyProfile <- function(profiles) {
  profiles$has_terminal <- (profiles$n_alt_tss > 0L) |
    (profiles$n_alt_tes > 0L)
  profiles$has_boundary <- (profiles$n_a5ss > 0L) |
    (profiles$n_a3ss > 0L)
  profiles$has_inclusion <- (profiles$n_se > 0L) |
    (profiles$n_missing_internal > 0L)
  profiles$has_partial_ret <- profiles$n_partial_ir > 0L
  profiles$has_full_ret <- (profiles$n_ir + profiles$n_ir_diff) > 0L

  profiles$n_categories <- profiles$has_terminal + profiles$has_boundary +
    profiles$has_inclusion + profiles$has_partial_ret + profiles$has_full_ret

  profiles$profile_type <- dplyr::case_when(
    profiles$n_events == 0L ~ "Identical",
    profiles$n_categories == 1L & profiles$has_terminal ~ "Terminal-only",
    profiles$n_categories == 1L & profiles$has_boundary ~ "Boundary-only",
    profiles$n_categories == 1L & profiles$has_inclusion ~ "Inclusion-only",
    profiles$n_categories == 1L & profiles$has_partial_ret ~
      "Partial_Retention-only",
    profiles$n_categories == 1L & profiles$has_full_ret ~
      "Full_Retention-only",
    profiles$n_categories >= 2L ~ "Combined",
    TRUE ~ "Unclassified"
  )

  # Remove helper columns
  profiles$has_terminal <- NULL
  profiles$has_boundary <- NULL
  profiles$has_inclusion <- NULL
  profiles$has_partial_ret <- NULL
  profiles$has_full_ret <- NULL
  profiles$n_categories <- NULL

  profiles
}


#' Compute Cliff's delta (nonparametric effect size)
#'
#' Uses a rank-based O((n+m) log(n+m)) algorithm instead of the naive
#' O(n*m) outer product. For each element in x, counts how many in y
#' are less vs greater using sorted ranks, avoiding the n*m matrix.
#'
#' @param x,y Numeric vectors.
#' @return Cliff's delta in the range -1 to 1, or NA_real_ if either is empty.
#' @keywords internal
.cliffsD <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx == 0L || ny == 0L) return(NA_real_)

  # Rank-based approach: combine, rank, sum ranks of x
  # Cliff's D = (2 * sum_rank_x - nx*(nx+1)) / (nx*ny) - 1
  # where sum_rank_x is sum of ranks of x in the combined vector,
  # using average ranks for ties
  combined <- c(x, y)
  ranks <- rank(combined, ties.method = "average")
  sum_rank_x <- sum(ranks[seq_len(nx)])

  # Mann-Whitney U = sum_rank_x - nx*(nx+1)/2
  # Cliff's D = (2*U)/(nx*ny) - 1
  u <- sum_rank_x - nx * (nx + 1) / 2
  (2 * u) / (nx * ny) - 1
}


# --- Section helpers ---

#' @keywords internal
.compareEventFrequency <- function(profiles_a, profiles_b) {
  a <- .addEventBinary(profiles_a)
  b <- .addEventBinary(profiles_b)

  event_flags <- c("has_alt_tss", "has_alt_tes", "has_a5ss", "has_a3ss",
                    "has_se", "has_missing_internal", "has_ir",
                    "has_partial_ir")
  event_names <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS",
                    "SE", "Missing_Internal", "IR", "Partial_IR")

  n_a <- nrow(a)
  n_b <- nrow(b)
  per_event <- list()
  raw_p <- numeric(length(event_flags))

  for (i in seq_along(event_flags)) {
    flag <- event_flags[i]
    count_a <- sum(a[[flag]])
    count_b <- sum(b[[flag]])

    tbl <- matrix(
      c(count_a, n_a - count_a, count_b, n_b - count_b),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("set_a", "set_b"), c("present", "absent"))
    )

    ft <- tryCatch(stats::fisher.test(tbl), error = function(e) NULL)

    if (!is.null(ft)) {
      or <- as.numeric(ft$estimate)
      direction <- if (ft$p.value >= 0.05) {
        "NS"
      } else if (or > 1) {
        "enriched in A"
      } else {
        "enriched in B"
      }
      per_event[[event_names[i]]] <- list(
        event_type = event_names[i],
        odds_ratio = or,
        conf_int = as.numeric(ft$conf.int),
        p_value = ft$p.value,
        direction = direction,
        count_a = count_a, count_b = count_b,
        n_a = n_a, n_b = n_b
      )
      raw_p[i] <- ft$p.value
    } else {
      per_event[[event_names[i]]] <- list(
        event_type = event_names[i],
        odds_ratio = NA_real_,
        conf_int = c(NA_real_, NA_real_),
        p_value = NA_real_,
        direction = "NS",
        count_a = count_a, count_b = count_b,
        n_a = n_a, n_b = n_b
      )
      raw_p[i] <- NA_real_
    }
  }

  adj_p <- stats::p.adjust(raw_p, method = "fdr")
  for (i in seq_along(event_flags)) {
    per_event[[event_names[i]]]$adj_p_value <- adj_p[i]
  }

  list(
    description = "Event type prevalence: Fisher's exact per type with FDR",
    test = "Fisher's exact test (8 tests, BH-adjusted)",
    per_event = per_event,
    adj_p_values = stats::setNames(adj_p, event_names),
    n_a = n_a, n_b = n_b,
    data_a = profiles_a, data_b = profiles_b
  )
}


#' @keywords internal
.compareDistribution <- function(x_a, x_b, description, n_a, n_b,
                                  data_a, data_b) {
  wt <- tryCatch(
    suppressWarnings(stats::wilcox.test(x_a, x_b)),
    error = function(e) NULL
  )
  cd <- .cliffsD(x_a, x_b)

  list(
    description = description,
    test = "Wilcoxon rank-sum test",
    p_value = if (!is.null(wt)) wt$p.value else NA_real_,
    cliffs_d = cd,
    median_a = stats::median(x_a, na.rm = TRUE),
    median_b = stats::median(x_b, na.rm = TRUE),
    n_a = n_a, n_b = n_b,
    data_a = data_a, data_b = data_b
  )
}


#' @keywords internal
.compareCooccurrence <- function(profiles_a, profiles_b) {
  cooc_a <- testCooccurrence(profiles_a)
  cooc_b <- testCooccurrence(profiles_b)

  # Haldane correction: add 0.5 to all cells when any cell is zero,

  # then recompute OR. This preserves strong negative associations (OR < 1)
  # that pmax(OR, 0.5) would censor.
  .haldane_log_or <- function(cooc) {
    a <- cooc$n_both
    b <- cooc$n_a_only
    c <- cooc$n_b_only
    d <- cooc$n_neither
    needs_correction <- (a == 0L | b == 0L | c == 0L | d == 0L)
    a[needs_correction] <- a[needs_correction] + 0.5
    b[needs_correction] <- b[needs_correction] + 0.5
    c[needs_correction] <- c[needs_correction] + 0.5
    d[needs_correction] <- d[needs_correction] + 0.5
    log((a * d) / (b * c))
  }
  cooc_a$log_or <- .haldane_log_or(cooc_a)
  cooc_b$log_or <- .haldane_log_or(cooc_b)

  # Match by event pair
  merged <- dplyr::inner_join(
    cooc_a[, c("event_a", "event_b", "log_or")],
    cooc_b[, c("event_a", "event_b", "log_or")],
    by = c("event_a", "event_b"),
    suffix = c("_a", "_b")
  )

  # Filter to finite values
  valid <- is.finite(merged$log_or_a) & is.finite(merged$log_or_b)
  merged <- merged[valid, ]
  n_pairs <- nrow(merged)

  if (n_pairs >= 3L) {
    ct <- tryCatch(
      suppressWarnings(
        stats::cor.test(merged$log_or_a, merged$log_or_b, method = "pearson")
      ),
      error = function(e) NULL
    )
    r <- if (!is.null(ct)) as.numeric(ct$estimate) else NA_real_
    p <- if (!is.null(ct)) ct$p.value else NA_real_
  } else {
    r <- NA_real_
    p <- NA_real_
  }

  list(
    description = "Co-occurrence pattern similarity (log-odds correlation)",
    test = "Pearson correlation of log-odds ratios",
    p_value = p,
    correlation = r,
    n_pairs_compared = n_pairs,
    n_a = nrow(profiles_a), n_b = nrow(profiles_b),
    data_a = cooc_a, data_b = cooc_b
  )
}


#' @keywords internal
.compareChiSquare <- function(cats_a, cats_b, description, data_a, data_b) {
  all_levels <- sort(unique(c(cats_a, cats_b)))
  tbl_a <- table(factor(cats_a, levels = all_levels))
  tbl_b <- table(factor(cats_b, levels = all_levels))
  ctable <- rbind(set_a = as.integer(tbl_a), set_b = as.integer(tbl_b))
  colnames(ctable) <- all_levels

  # Remove zero-sum columns
  keep <- colSums(ctable) > 0L
  ctable <- ctable[, keep, drop = FALSE]

  n_total <- sum(ctable)
  k <- min(nrow(ctable), ncol(ctable))

  if (ncol(ctable) < 2L || n_total == 0L) {
    return(list(
      description = description,
      test = "Chi-square test",
      p_value = NA_real_,
      cramers_v = NA_real_,
      n_a = sum(tbl_a), n_b = sum(tbl_b),
      table = ctable,
      data_a = data_a, data_b = data_b
    ))
  }

  # Use simulated p-value when expected counts < 5
  expected <- tryCatch(
    suppressWarnings(stats::chisq.test(ctable)$expected),
    error = function(e) NULL
  )
  use_sim <- !is.null(expected) && min(expected) < 5

  chi <- tryCatch(
    if (use_sim) {
      stats::chisq.test(ctable, simulate.p.value = TRUE, B = 10000L)
    } else {
      stats::chisq.test(ctable)
    },
    error = function(e) NULL
  )

  cramers_v <- if (!is.null(chi) && k > 1L) {
    as.numeric(sqrt(chi$statistic / (n_total * (k - 1L))))
  } else {
    NA_real_
  }

  list(
    description = description,
    test = if (use_sim) "Chi-square test (simulated p-value)" else
      "Chi-square test",
    p_value = if (!is.null(chi)) chi$p.value else NA_real_,
    cramers_v = cramers_v,
    n_a = as.integer(sum(tbl_a)),
    n_b = as.integer(sum(tbl_b)),
    table = ctable,
    data_a = data_a, data_b = data_b
  )
}


#' @keywords internal
.compareTopology <- function(profiles_a, profiles_b,
                              structures_a, structures_b) {
  topo_a <- classifyTopology(profiles_a, structures_a, test = FALSE)
  topo_b <- classifyTopology(profiles_b, structures_b, test = FALSE)

  cats_a <- topo_a$classifications$topology
  cats_b <- topo_b$classifications$topology

  if (length(cats_a) == 0L && length(cats_b) == 0L) {
    return(list(
      description = "Boundary event topology (F2F/B2B/Distal)",
      test = "Chi-square test",
      p_value = NA_real_,
      cramers_v = NA_real_,
      n_a = 0L, n_b = 0L,
      table = matrix(nrow = 0, ncol = 0),
      data_a = topo_a$classifications,
      data_b = topo_b$classifications
    ))
  }

  result <- .compareChiSquare(
    cats_a, cats_b,
    description = "Boundary event topology (F2F/B2B/Distal)",
    data_a = topo_a$classifications,
    data_b = topo_b$classifications
  )
  result
}


#' @keywords internal
.compareProximity <- function(profiles_a, profiles_b,
                               structures_a, structures_b, n_perm) {
  prox_a <- testProximity(profiles_a, structures_a, n_perm = n_perm)
  prox_b <- testProximity(profiles_b, structures_b, n_perm = n_perm)

  z_diff <- if (!is.na(prox_a$z_score) && !is.na(prox_b$z_score)) {
    prox_a$z_score - prox_b$z_score
  } else {
    NA_real_
  }

  list(
    description = "Spatial clustering of events (descriptive comparison)",
    test = "Permutation-based z-score comparison (no formal test for difference)",
    p_value = NA_real_,
    z_diff = z_diff,
    z_a = prox_a$z_score,
    z_b = prox_b$z_score,
    p_a = prox_a$p_value,
    p_b = prox_b$p_value,
    n_a = prox_a$n_profiles_used,
    n_b = prox_b$n_profiles_used,
    data_a = prox_a,
    data_b = prox_b
  )
}


#' @keywords internal
.comparePtcRate <- function(ptc_a, ptc_b) {
  # PTC status for comparator isoforms
  a_ptc <- sum(ptc_a$has_ptc, na.rm = TRUE)
  a_total <- sum(!is.na(ptc_a$has_ptc))
  b_ptc <- sum(ptc_b$has_ptc, na.rm = TRUE)
  b_total <- sum(!is.na(ptc_b$has_ptc))

  if (a_total == 0L || b_total == 0L) {
    return(list(
      description = "PTC prevalence in comparator isoforms",
      test = "Fisher's exact test",
      p_value = NA_real_,
      odds_ratio = NA_real_,
      conf_int = c(NA_real_, NA_real_),
      direction = "NS",
      n_a = a_total, n_b = b_total,
      data_a = ptc_a, data_b = ptc_b
    ))
  }

  tbl <- matrix(
    c(a_ptc, a_total - a_ptc, b_ptc, b_total - b_ptc),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("set_a", "set_b"), c("PTC", "no_PTC"))
  )

  ft <- tryCatch(stats::fisher.test(tbl), error = function(e) NULL)

  or <- if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_
  p <- if (!is.null(ft)) ft$p.value else NA_real_
  ci <- if (!is.null(ft)) as.numeric(ft$conf.int) else c(NA_real_, NA_real_)

  direction <- if (is.na(p) || p >= 0.05) {
    "NS"
  } else if (or > 1) {
    "enriched in A"
  } else {
    "enriched in B"
  }

  list(
    description = "PTC prevalence in comparator isoforms",
    test = "Fisher's exact test",
    p_value = p,
    odds_ratio = or,
    conf_int = ci,
    direction = direction,
    n_a = a_total, n_b = b_total,
    data_a = ptc_a, data_b = ptc_b
  )
}


#' @keywords internal
.comparePtcDistance <- function(ptc_a, ptc_b) {
  dist_a <- ptc_a$ptc_distance[ptc_a$has_ptc & !is.na(ptc_a$ptc_distance)]
  dist_b <- ptc_b$ptc_distance[ptc_b$has_ptc & !is.na(ptc_b$ptc_distance)]

  .compareDistribution(
    dist_a, dist_b,
    description = "PTC-to-last-EJC distance for PTC+ isoforms",
    n_a = length(dist_a), n_b = length(dist_b),
    data_a = ptc_a[ptc_a$has_ptc & !is.na(ptc_a$ptc_distance), ],
    data_b = ptc_b[ptc_b$has_ptc & !is.na(ptc_b$ptc_distance), ]
  )
}


#' @keywords internal
.compareProfileCategories <- function(profiles_a, profiles_b) {
  a <- .classifyProfile(profiles_a)
  b <- .classifyProfile(profiles_b)

  .compareChiSquare(
    a$profile_type, b$profile_type,
    description = "Splicing pattern category distribution",
    data_a = a, data_b = b
  )
}


#' @keywords internal
.runPairedTests <- function(profiles_a, profiles_b, min_overlap, verbose) {
  # Match on gene_id + reference_isoform_id
  overlap <- dplyr::inner_join(
    profiles_a[, c("gene_id", "reference_isoform_id")],
    profiles_b[, c("gene_id", "reference_isoform_id")],
    by = c("gene_id", "reference_isoform_id")
  )
  overlap <- dplyr::distinct(overlap)

  n_overlap <- nrow(overlap)

  if (n_overlap < min_overlap) {
    if (verbose) {
      message(sprintf("Paired tests skipped: %d gene overlap < %d minimum",
                       n_overlap, min_overlap))
    }
    return(list(
      skipped = TRUE,
      skipped_reason = sprintf(
        "Gene overlap (%d) below minimum (%d)", n_overlap, min_overlap
      ),
      n_genes_overlap = n_overlap
    ))
  }

  if (verbose) message(sprintf("Running paired analysis (%d genes)", n_overlap))

  # Get paired profiles
  a_paired <- dplyr::inner_join(
    profiles_a, overlap, by = c("gene_id", "reference_isoform_id")
  )
  b_paired <- dplyr::inner_join(
    profiles_b, overlap, by = c("gene_id", "reference_isoform_id")
  )

  # Ensure 1:1 by taking first match per gene
  a_paired <- a_paired[!duplicated(
    paste(a_paired$gene_id, a_paired$reference_isoform_id)
  ), ]
  b_paired <- b_paired[!duplicated(
    paste(b_paired$gene_id, b_paired$reference_isoform_id)
  ), ]

  # Re-align rows
  key_a <- paste(a_paired$gene_id, a_paired$reference_isoform_id)
  key_b <- paste(b_paired$gene_id, b_paired$reference_isoform_id)
  common_keys <- intersect(key_a, key_b)
  a_paired <- a_paired[match(common_keys, key_a), ]
  b_paired <- b_paired[match(common_keys, key_b), ]
  n_paired <- length(common_keys)

  # Paired Wilcoxon signed-rank on n_events
  diff_events <- a_paired$n_events - b_paired$n_events
  pwt <- tryCatch(
    suppressWarnings(
      stats::wilcox.test(a_paired$n_events, b_paired$n_events, paired = TRUE)
    ),
    error = function(e) NULL
  )

  # Sign test
  n_more <- sum(diff_events > 0L)
  n_fewer <- sum(diff_events < 0L)
  n_tied <- sum(diff_events == 0L)
  sign_p <- if (n_more + n_fewer > 0L) {
    stats::binom.test(n_more, n_more + n_fewer, p = 0.5)$p.value
  } else {
    NA_real_
  }

  # McNemar per event type
  a_bin <- .addEventBinary(a_paired)
  b_bin <- .addEventBinary(b_paired)

  event_flags <- c("has_alt_tss", "has_alt_tes", "has_a5ss", "has_a3ss",
                    "has_se", "has_missing_internal", "has_ir",
                    "has_partial_ir")
  event_names <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS",
                    "SE", "Missing_Internal", "IR", "Partial_IR")

  mcnemar_tests <- list()
  mcnemar_p <- numeric(length(event_flags))

  for (i in seq_along(event_flags)) {
    flag <- event_flags[i]
    a_yes_b_no <- sum(a_bin[[flag]] == 1L & b_bin[[flag]] == 0L)
    a_no_b_yes <- sum(a_bin[[flag]] == 0L & b_bin[[flag]] == 1L)
    discordant <- a_yes_b_no + a_no_b_yes

    if (discordant >= 10L) {
      m <- matrix(c(
        sum(a_bin[[flag]] == 1L & b_bin[[flag]] == 1L), a_yes_b_no,
        a_no_b_yes, sum(a_bin[[flag]] == 0L & b_bin[[flag]] == 0L)
      ), nrow = 2)
      mt <- tryCatch(stats::mcnemar.test(m), error = function(e) NULL)
      mp <- if (!is.null(mt)) mt$p.value else NA_real_
    } else if (discordant > 0L) {
      mp <- stats::binom.test(a_yes_b_no, discordant, p = 0.5)$p.value
    } else {
      mp <- NA_real_
    }

    mcnemar_p[i] <- mp
    mcnemar_tests[[event_names[i]]] <- list(
      event_type = event_names[i],
      p_value = mp,
      a_only = a_yes_b_no,
      b_only = a_no_b_yes,
      n_discordant = discordant
    )
  }

  # FDR correction for McNemar tests
  adj_mcnemar <- stats::p.adjust(mcnemar_p, method = "fdr")
  for (i in seq_along(event_flags)) {
    mcnemar_tests[[event_names[i]]]$adj_p_value <- adj_mcnemar[i]
  }

  # Profile type transitions
  a_typed <- .classifyProfile(a_paired)
  b_typed <- .classifyProfile(b_paired)
  n_same <- sum(a_typed$profile_type == b_typed$profile_type)
  n_diff <- sum(a_typed$profile_type != b_typed$profile_type)

  list(
    n_genes_overlap = n_paired,
    paired_wilcox = list(
      description = "Paired event count comparison",
      test = "Wilcoxon signed-rank test",
      p_value = if (!is.null(pwt)) pwt$p.value else NA_real_,
      median_diff = stats::median(diff_events)
    ),
    sign_test = list(
      description = "Direction of event count difference",
      test = "Binomial sign test",
      p_value = sign_p,
      n_a_more = n_more,
      n_b_more = n_fewer,
      n_tied = n_tied,
      pct_a_more = if (n_more + n_fewer > 0L) {
        100 * n_more / (n_more + n_fewer)
      } else {
        NA_real_
      }
    ),
    mcnemar_tests = mcnemar_tests,
    type_transitions = list(
      n_same_type = n_same,
      n_diff_type = n_diff,
      pct_type_change = 100 * n_diff / n_paired
    )
  )
}


# ============================================================================
# Meta-analysis helpers
# ============================================================================

#' @keywords internal
.runMetaAnalysis <- function(yi, sei, comp_names, section_name, verbose) {
  valid <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[valid]
  sei <- sei[valid]
  names_valid <- comp_names[valid]

  if (length(yi) < 2L) {
    return(list(
      skipped = TRUE,
      skipped_reason = sprintf(
        "%s: fewer than 2 valid comparisons (%d)", section_name, length(yi)
      )
    ))
  }

  re <- tryCatch(
    metafor::rma(yi = yi, sei = sei, method = "REML"),
    error = function(e) NULL
  )

  if (is.null(re)) {
    return(list(
      skipped = TRUE,
      skipped_reason = sprintf("%s: metafor::rma() failed", section_name)
    ))
  }

  if (verbose) {
    message(sprintf("  %s: pooled = %.3f (%.3f, %.3f), I2 = %.1f%%",
                     section_name, re$beta, re$ci.lb, re$ci.ub, re$I2))
  }

  list(
    estimate = as.numeric(re$beta),
    se = re$se,
    p_value = re$pval,
    ci_lower = re$ci.lb,
    ci_upper = re$ci.ub,
    Q = re$QE,
    I2 = re$I2,
    tau2 = re$tau2,
    k = re$k,
    per_comparison = tibble::tibble(
      comparison = names_valid,
      effect_size = yi,
      se = sei
    )
  )
}


#' @keywords internal
.metaCliffD <- function(comparison_list, comp_names, section, verbose) {
  yi <- numeric(length(comp_names))
  sei <- numeric(length(comp_names))

  for (i in seq_along(comp_names)) {
    res <- comparison_list[[comp_names[i]]][[section]]
    if (is.null(res) || isTRUE(res$skipped)) {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
      next
    }
    d <- res$cliffs_d
    n_a <- res$n_a
    n_b <- res$n_b
    n_eff <- 2 * n_a * n_b / (n_a + n_b)
    yi[i] <- if (!is.na(d)) d else NA_real_
    sei[i] <- if (!is.na(d)) {
      sqrt(max(1 - d^2, 0.01) / max(n_eff - 1, 1))
    } else {
      NA_real_
    }
  }

  .runMetaAnalysis(yi, sei, comp_names, section, verbose)
}


#' @keywords internal
.metaEventFrequency <- function(comparison_list, comp_names, verbose) {
  event_names <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS",
                    "SE", "Missing_Internal", "IR", "Partial_IR")

  results <- list()
  for (et in event_names) {
    yi <- numeric(length(comp_names))
    sei <- numeric(length(comp_names))

    for (i in seq_along(comp_names)) {
      res <- comparison_list[[comp_names[i]]]$event_frequency
      if (is.null(res) || isTRUE(res$skipped)) {
        yi[i] <- NA_real_
        sei[i] <- NA_real_
        next
      }
      ev <- res$per_event[[et]]
      if (is.null(ev) || is.na(ev$odds_ratio) || ev$odds_ratio <= 0) {
        yi[i] <- NA_real_
        sei[i] <- NA_real_
        next
      }
      yi[i] <- log(ev$odds_ratio)
      ci <- ev$conf_int
      if (all(is.finite(ci)) && ci[1] > 0 && ci[2] > 0) {
        sei[i] <- (log(ci[2]) - log(ci[1])) / (2 * 1.96)
      } else {
        sei[i] <- NA_real_
      }
    }

    results[[et]] <- .runMetaAnalysis(yi, sei, comp_names, et, verbose)

    # Back-transform pooled estimate to OR scale
    if (!isTRUE(results[[et]]$skipped)) {
      results[[et]]$pooled_or <- exp(results[[et]]$estimate)
      results[[et]]$pooled_or_ci <- c(
        exp(results[[et]]$ci_lower),
        exp(results[[et]]$ci_upper)
      )
    }
  }

  results
}


#' @keywords internal
.metaLogOR <- function(comparison_list, comp_names, section, verbose) {
  yi <- numeric(length(comp_names))
  sei <- numeric(length(comp_names))

  for (i in seq_along(comp_names)) {
    res <- comparison_list[[comp_names[i]]][[section]]
    if (is.null(res) || isTRUE(res$skipped)) {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
      next
    }
    or <- res$odds_ratio
    ci <- res$conf_int
    if (!is.na(or) && or > 0 && all(is.finite(ci)) && ci[1] > 0) {
      yi[i] <- log(or)
      sei[i] <- (log(ci[2]) - log(ci[1])) / (2 * 1.96)
    } else {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
    }
  }

  result <- .runMetaAnalysis(yi, sei, comp_names, section, verbose)
  if (!isTRUE(result$skipped)) {
    result$pooled_or <- exp(result$estimate)
    result$pooled_or_ci <- c(exp(result$ci_lower), exp(result$ci_upper))
  }
  result
}


#' @keywords internal
.metaCramersV <- function(comparison_list, comp_names, section, verbose) {
  yi <- numeric(length(comp_names))
  sei <- numeric(length(comp_names))

  for (i in seq_along(comp_names)) {
    res <- comparison_list[[comp_names[i]]][[section]]
    if (is.null(res) || isTRUE(res$skipped)) {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
      next
    }
    v <- res$cramers_v
    n_total <- res$n_a + res$n_b
    yi[i] <- if (!is.na(v)) v else NA_real_
    sei[i] <- if (n_total > 0L) max(sqrt(2 / n_total), 0.01) else NA_real_
  }

  .runMetaAnalysis(yi, sei, comp_names, section, verbose)
}


#' @keywords internal
.metaCorrelation <- function(comparison_list, comp_names, verbose) {
  yi <- numeric(length(comp_names))
  sei <- numeric(length(comp_names))

  for (i in seq_along(comp_names)) {
    res <- comparison_list[[comp_names[i]]]$cooccurrence
    if (is.null(res) || isTRUE(res$skipped)) {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
      next
    }
    r <- res$correlation
    n_pairs <- res$n_pairs_compared
    if (!is.na(r) && n_pairs >= 4L) {
      # Fisher's z-transform
      yi[i] <- 0.5 * log((1 + r) / (1 - r))
      sei[i] <- 1 / sqrt(n_pairs - 3)
    } else {
      yi[i] <- NA_real_
      sei[i] <- NA_real_
    }
  }

  .runMetaAnalysis(yi, sei, comp_names, "cooccurrence", verbose)
}
