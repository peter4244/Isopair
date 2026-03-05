#' Test Event Co-occurrence
#'
#' Performs Fisher's exact test for each pair of event types to assess
#' whether they co-occur more or less often than expected by chance.
#' Uses 8 binary event indicators (C(8,2) = 28 pairwise tests). IR and
#' IR_diff are merged into a single indicator because they represent the
#' same biological mechanism.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with event
#'   count columns (n_a5ss, n_a3ss, n_se, etc.) and/or tss_changed,
#'   tes_changed.
#' @param stratify_by Optional character; column name in \code{profiles} to
#'   stratify tests by. When provided, tests are run within each stratum.
#' @return A tibble with columns: \code{event_a}, \code{event_b},
#'   \code{odds_ratio}, \code{p_value}, \code{adj_p_value}, \code{conf_low},
#'   \code{conf_high}, \code{n_both}, \code{n_a_only}, \code{n_b_only},
#'   \code{n_neither}, \code{significant}. Under stratification, adds
#'   \code{stratum} column.
#' @examples
#' data(example_profiles)
#' cooc <- testCooccurrence(example_profiles)
#' cooc[, c("event_a", "event_b", "odds_ratio", "p_value")]
#' @export
#' @importFrom stats fisher.test p.adjust
#' @importFrom dplyr bind_rows mutate
#' @importFrom tibble tibble
testCooccurrence <- function(profiles, stratify_by = NULL) {

  # Add binary event indicators
  profiles <- .addEventBinary(profiles)

  event_flags <- c("has_alt_tss", "has_alt_tes", "has_a5ss", "has_a3ss",
                    "has_se", "has_missing_internal", "has_ir", "has_partial_ir")
  event_names <- c("Alt_TSS", "Alt_TES", "A5SS", "A3SS",
                    "SE", "Missing_Internal", "IR", "Partial_IR")
  names(event_names) <- event_flags

  # Generate all pairwise combinations
  pairs <- utils::combn(seq_along(event_flags), 2, simplify = FALSE)

  if (is.null(stratify_by)) {
    results <- .runCooccurrenceTests(profiles, pairs, event_flags, event_names)
    results$adj_p_value <- stats::p.adjust(results$p_value, method = "fdr")
    results$significant <- results$adj_p_value < 0.05
  } else {
    if (!stratify_by %in% names(profiles)) {
      stop(sprintf("Column '%s' not found in profiles", stratify_by))
    }
    strata <- unique(profiles[[stratify_by]])
    result_list <- lapply(strata, function(s) {
      subset_data <- profiles[profiles[[stratify_by]] == s, ]
      if (nrow(subset_data) < 2L) return(NULL)
      res <- .runCooccurrenceTests(subset_data, pairs, event_flags, event_names)
      res$stratum <- s
      res$adj_p_value <- stats::p.adjust(res$p_value, method = "fdr")
      res$significant <- res$adj_p_value < 0.05
      res
    })
    results <- dplyr::bind_rows(result_list)
  }

  results
}


#' Add Binary Event Presence Indicators to Profiles
#'
#' Converts per-type event counts to binary (0/1) indicators. IR and
#' IR_diff are merged into a single \code{has_ir} indicator.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()}.
#' @return The input tibble with added binary columns: \code{has_alt_tss},
#'   \code{has_alt_tes}, \code{has_a5ss}, \code{has_a3ss}, \code{has_se},
#'   \code{has_missing_internal}, \code{has_ir}, \code{has_partial_ir}.
#' @keywords internal
.addEventBinary <- function(profiles) {
  profiles$has_alt_tss <- as.integer(profiles$tss_changed)
  profiles$has_alt_tes <- as.integer(profiles$tes_changed)
  profiles$has_a5ss <- as.integer(profiles$n_a5ss > 0L)
  profiles$has_a3ss <- as.integer(profiles$n_a3ss > 0L)
  profiles$has_se <- as.integer(profiles$n_se > 0L)
  profiles$has_missing_internal <- as.integer(profiles$n_missing_internal > 0L)
  profiles$has_ir <- as.integer((profiles$n_ir + profiles$n_ir_diff) > 0L)
  profiles$has_partial_ir <- as.integer(profiles$n_partial_ir > 0L)
  profiles
}


#' Run pairwise co-occurrence Fisher's tests
#' @keywords internal
.runCooccurrenceTests <- function(data, pairs, event_flags, event_names) {
  results <- vector("list", length(pairs))
  for (k in seq_along(pairs)) {
    idx <- pairs[[k]]
    flag_a <- event_flags[idx[1]]
    flag_b <- event_flags[idx[2]]

    tbl <- table(
      factor(data[[flag_a]], levels = c(0L, 1L)),
      factor(data[[flag_b]], levels = c(0L, 1L))
    )

    ft <- tryCatch(
      stats::fisher.test(tbl),
      error = function(e) NULL
    )

    if (is.null(ft)) {
      results[[k]] <- tibble::tibble(
        event_a = event_names[flag_a],
        event_b = event_names[flag_b],
        odds_ratio = NA_real_, p_value = NA_real_,
        conf_low = NA_real_, conf_high = NA_real_,
        n_both = as.integer(tbl[2, 2]),
        n_a_only = as.integer(tbl[2, 1]),
        n_b_only = as.integer(tbl[1, 2]),
        n_neither = as.integer(tbl[1, 1])
      )
    } else {
      results[[k]] <- tibble::tibble(
        event_a = event_names[flag_a],
        event_b = event_names[flag_b],
        odds_ratio = as.numeric(ft$estimate),
        p_value = ft$p.value,
        conf_low = ft$conf.int[1],
        conf_high = ft$conf.int[2],
        n_both = as.integer(tbl[2, 2]),
        n_a_only = as.integer(tbl[2, 1]),
        n_b_only = as.integer(tbl[1, 2]),
        n_neither = as.integer(tbl[1, 1])
      )
    }
  }

  dplyr::bind_rows(results)
}
