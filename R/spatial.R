#' Test Positional Bias of Splicing Events
#'
#' For each event type, tests whether events are uniformly distributed along
#' the gene span using a Kolmogorov-Smirnov test against Uniform(0,1).
#' Event positions are computed as relative positions within the reference
#' isoform's genomic span.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with a
#'   \code{detailed_events} list column.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}
#'   with \code{isoform_id}, \code{tx_start}, \code{tx_end}.
#' @param min_events Integer; minimum number of events per type to run the
#'   test. Default 30.
#' @return A tibble with columns: \code{event_type}, \code{n_events},
#'   \code{median_position}, \code{mean_position}, \code{ks_statistic},
#'   \code{p_value}, \code{adj_p_value}, \code{significant},
#'   \code{bias_direction}.
#' @examples
#' data(example_profiles)
#' data(example_structures)
#' # Lower min_events for small demo data
#' bias <- testPositionalBias(example_profiles, example_structures,
#'   min_events = 1L)
#' bias
#' @export
#' @importFrom stats ks.test p.adjust
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
testPositionalBias <- function(profiles, structures, min_events = 30L) {

  # Extract events with genomic coordinates
  events <- .extractEventsWithPositions(profiles, structures)

  if (nrow(events) == 0L) {
    return(tibble::tibble(
      event_type = character(0), n_events = integer(0),
      median_position = numeric(0), mean_position = numeric(0),
      ks_statistic = numeric(0), p_value = numeric(0),
      adj_p_value = numeric(0), significant = logical(0),
      bias_direction = character(0)
    ))
  }

  event_types <- unique(events$event_type)
  results <- list()

  for (et in event_types) {
    positions <- events$relative_pos[events$event_type == et]
    n <- length(positions)

    if (n < min_events) next

    ks <- suppressWarnings(stats::ks.test(positions, "punif", 0, 1))

    med <- stats::median(positions)
    bias <- if (ks$p.value >= 0.05) {
      "None"
    } else if (med < 0.3) {
      "5' biased"
    } else if (med > 0.7) {
      "3' biased"
    } else {
      "Central/other"
    }

    results[[length(results) + 1L]] <- tibble::tibble(
      event_type = et,
      n_events = n,
      median_position = med,
      mean_position = mean(positions),
      ks_statistic = as.numeric(ks$statistic),
      p_value = ks$p.value,
      adj_p_value = NA_real_,
      significant = NA,
      bias_direction = bias
    )
  }

  if (length(results) == 0L) {
    return(tibble::tibble(
      event_type = character(0), n_events = integer(0),
      median_position = numeric(0), mean_position = numeric(0),
      ks_statistic = numeric(0), p_value = numeric(0),
      adj_p_value = numeric(0), significant = logical(0),
      bias_direction = character(0)
    ))
  }

  result <- dplyr::bind_rows(results)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = "fdr")
  result$significant <- result$adj_p_value < 0.05

  # Update bias_direction based on adjusted p-value
  result$bias_direction <- ifelse(
    !result$significant, "None", result$bias_direction
  )

  result
}


#' Classify Boundary Event Topology
#'
#' For profiles containing boundary events at both the 5' and 3' splice sites,
#' classifies each 5'-type x 3'-type pair into topology categories:
#' \itemize{
#'   \item \strong{F2F} (Face-to-Face): events flank the same intron
#'   \item \strong{B2B} (Back-to-Back): events are on the same exon
#'   \item \strong{Distal}: events are on non-adjacent structural elements
#' }
#'
#' By default, only A5SS and A3SS events are classified. Set
#' \code{five_prime_types} and \code{three_prime_types} to include additional
#' boundary event types such as Partial_IR_5 and Partial_IR_3.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param five_prime_types Character vector of event types affecting the 5'
#'   splice site (exon donor boundary). Default \code{"A5SS"}. Partial_IR_5
#'   also affects this boundary.
#' @param three_prime_types Character vector of event types affecting the 3'
#'   splice site (exon acceptor boundary). Default \code{"A3SS"}. Partial_IR_3
#'   also affects this boundary.
#' @param test Logical; if TRUE (default), run permutation test for F2F
#'   enrichment.
#' @param n_perm Integer; number of permutations for the topology test.
#'   Default 1000.
#' @return A named list with:
#'   \describe{
#'     \item{classifications}{Tibble with \code{gene_id},
#'       \code{reference_isoform_id}, \code{comparator_isoform_id},
#'       \code{event_a_type}, \code{event_b_type}, \code{topology}.}
#'     \item{summary}{Tibble with \code{topology}, \code{count},
#'       \code{proportion}, and optionally \code{permutation_p_value}.}
#'   }
#' @examples
#' data(example_profiles)
#' data(example_structures)
#' topo <- classifyTopology(example_profiles, example_structures,
#'   test = FALSE)
#' topo$summary
#'
#' # Include partial intron retention events
#' topo2 <- classifyTopology(example_profiles, example_structures,
#'   five_prime_types = c("A5SS", "Partial_IR_5"),
#'   three_prime_types = c("A3SS", "Partial_IR_3"),
#'   test = FALSE)
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
classifyTopology <- function(profiles, structures,
                              five_prime_types = "A5SS",
                              three_prime_types = "A3SS",
                              test = TRUE, n_perm = 1000L) {

  all_classifications <- list()
  perm_data <- list()

  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next

    a5ss <- events[events$event_type %in% five_prime_types, , drop = FALSE]
    a3ss <- events[events$event_type %in% three_prime_types, , drop = FALSE]

    if (nrow(a5ss) == 0L || nrow(a3ss) == 0L) next

    ref_id <- profiles$reference_isoform_id[i]
    gene_id <- profiles$gene_id[i]

    # Look up reference isoform structure
    ref_struct <- structures[structures$isoform_id == ref_id, ]
    if (nrow(ref_struct) == 0L) next

    exon_starts <- ref_struct$exon_starts[[1]]
    exon_ends <- ref_struct$exon_ends[[1]]
    gene_strand <- ref_struct$strand

    n_exons <- length(exon_starts)
    if (n_exons < 2L) next

    # Classify each A5SS x A3SS pair
    for (a in seq_len(nrow(a5ss))) {
      for (b in seq_len(nrow(a3ss))) {
        topo <- .classifyEventPairTopology(
          a5ss[a, ], a3ss[b, ], exon_starts, exon_ends, gene_strand
        )
        all_classifications[[length(all_classifications) + 1L]] <-
          tibble::tibble(
            gene_id = gene_id,
            reference_isoform_id = ref_id,
            comparator_isoform_id = profiles$comparator_isoform_id[i],
            event_a_type = a5ss$event_type[a],
            event_b_type = a3ss$event_type[b],
            topology = topo
          )
      }
    }

    # Store data for permutation test
    if (test) {
      perm_data[[length(perm_data) + 1L]] <- list(
        n_a5ss = nrow(a5ss),
        n_a3ss = nrow(a3ss),
        n_exons = n_exons,
        gene_strand = gene_strand
      )
    }
  }

  classifications <- dplyr::bind_rows(all_classifications)

  if (nrow(classifications) == 0L) {
    summary_tbl <- tibble::tibble(
      topology = character(0), count = integer(0), proportion = numeric(0)
    )
    return(list(classifications = classifications, summary = summary_tbl))
  }

  # Summary counts
  topo_counts <- table(classifications$topology)
  total <- sum(topo_counts)
  summary_tbl <- tibble::tibble(
    topology = names(topo_counts),
    count = as.integer(topo_counts),
    proportion = as.numeric(topo_counts) / total
  )

  # Permutation test for F2F enrichment
  if (test && length(perm_data) > 0L) {
    observed_f2f_prop <- if ("F2F" %in% names(topo_counts)) {
      topo_counts["F2F"] / total
    } else {
      0
    }

    perm_f2f_props <- numeric(n_perm)
    for (p in seq_len(n_perm)) {
      perm_topos <- character(0)
      for (pd in perm_data) {
        a5_idx <- sample.int(pd$n_exons, pd$n_a5ss, replace = TRUE)
        a3_idx <- sample.int(pd$n_exons, pd$n_a3ss, replace = TRUE)
        for (ai in a5_idx) {
          for (bi in a3_idx) {
            if (ai == bi) {
              perm_topos <- c(perm_topos, "B2B")
            } else if (pd$gene_strand == "+" && ai == bi - 1L) {
              perm_topos <- c(perm_topos, "F2F")
            } else if (pd$gene_strand == "-" && ai == bi + 1L) {
              perm_topos <- c(perm_topos, "F2F")
            } else {
              perm_topos <- c(perm_topos, "Distal")
            }
          }
        }
      }
      n_perm_total <- length(perm_topos)
      perm_f2f_props[p] <- if (n_perm_total > 0L) {
        sum(perm_topos == "F2F") / n_perm_total
      } else {
        0
      }
    }

    perm_p <- mean(perm_f2f_props >= observed_f2f_prop)
    summary_tbl$permutation_p_value <- perm_p
  }

  list(classifications = classifications, summary = summary_tbl)
}


#' Summarize Boundary Event Length Distributions
#'
#' Extracts boundary events (A5SS, A3SS, Partial_IR_5, Partial_IR_3) and
#' computes length distributions (\code{abs(bp_diff)}). Returns per-type and
#' collapsed 5'/3' summary statistics.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with a
#'   \code{detailed_events} list column.
#' @return A list with three elements:
#'   \describe{
#'     \item{events}{Tibble of all boundary events with \code{event_type} and
#'       \code{length_bp} (= \code{abs(bp_diff)}).}
#'     \item{summary_stats}{Tibble per event type: \code{n}, \code{median},
#'       \code{mean}, \code{q25}, \code{q75}, \code{min}, \code{max}.}
#'     \item{combined_stats}{Same but with A5SS+Partial_IR_5 merged as
#'       \code{"5_prime"} and A3SS+Partial_IR_3 merged as \code{"3_prime"}.}
#'   }
#' @examples
#' data(example_profiles)
#' bl <- summarizeBoundaryLengths(example_profiles)
#' bl$summary_stats
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
summarizeBoundaryLengths <- function(profiles) {

  boundary_types <- c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")

  all_events <- list()
  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next
    boundary <- events[events$event_type %in% boundary_types, , drop = FALSE]
    if (nrow(boundary) == 0L) next
    all_events[[length(all_events) + 1L]] <- tibble::tibble(
      event_type = boundary$event_type,
      length_bp = abs(boundary$bp_diff)
    )
  }

  if (length(all_events) == 0L) {
    empty_events <- tibble::tibble(event_type = character(0),
                                    length_bp = numeric(0))
    empty_stats <- tibble::tibble(
      event_type = character(0), n = integer(0), median = numeric(0),
      mean = numeric(0), q25 = numeric(0), q75 = numeric(0),
      min = numeric(0), max = numeric(0)
    )
    return(list(events = empty_events, summary_stats = empty_stats,
                combined_stats = empty_stats))
  }

  events_tbl <- dplyr::bind_rows(all_events)

  # Per-type summary
  .compute_stats <- function(df) {
    types <- unique(df$event_type)
    do.call(rbind, lapply(types, function(et) {
      vals <- df$length_bp[df$event_type == et]
      qs <- stats::quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
      tibble::tibble(
        event_type = et, n = length(vals),
        median = stats::median(vals, na.rm = TRUE),
        mean = mean(vals, na.rm = TRUE),
        q25 = as.numeric(qs[1L]), q75 = as.numeric(qs[2L]),
        min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE)
      )
    }))
  }

  summary_stats <- .compute_stats(events_tbl)

  # Combined: collapse 5' and 3' categories
  combined <- events_tbl
  combined$event_type <- ifelse(
    combined$event_type %in% c("A5SS", "Partial_IR_5"), "5_prime",
    ifelse(combined$event_type %in% c("A3SS", "Partial_IR_3"), "3_prime",
           combined$event_type)
  )
  combined_stats <- .compute_stats(combined)

  list(events = events_tbl, summary_stats = summary_stats,
       combined_stats = combined_stats)
}


#' Expanded Multi-Level Topology Classification
#'
#' Classifies boundary event topology for each 5'-type x 3'-type pair
#' combination (e.g., A5SS x A3SS, A5SS x Partial_IR_3, etc.) as well as
#' a collapsed all-5'-types x all-3'-types analysis. Builds on
#' \code{\link{classifyTopology}()}.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param test Logical; if TRUE (default), run permutation tests.
#' @param n_perm Integer; number of permutations. Default 1000.
#' @return A list with:
#'   \describe{
#'     \item{per_type_pair}{Named list of \code{\link{classifyTopology}()}
#'       results for each 5' x 3' type pair that has data.}
#'     \item{collapsed}{\code{\link{classifyTopology}()} result with all 5'
#'       and 3' types combined.}
#'     \item{comparison}{Tibble comparing F2F proportion and p-value across
#'       all levels.}
#'   }
#' @examples
#' data(example_profiles)
#' data(example_structures)
#' te <- classifyTopologyExpanded(example_profiles, example_structures,
#'   test = FALSE)
#' te$comparison
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
classifyTopologyExpanded <- function(profiles, structures,
                                      test = TRUE, n_perm = 1000L) {

  five_types <- c("A5SS", "Partial_IR_5")
  three_types <- c("A3SS", "Partial_IR_3")

  # Per-type pair results
  per_type_pair <- list()
  for (ft in five_types) {
    for (tt in three_types) {
      pair_name <- paste0(ft, "_x_", tt)
      res <- tryCatch(
        classifyTopology(profiles, structures,
                          five_prime_types = ft, three_prime_types = tt,
                          test = test, n_perm = n_perm),
        error = function(e) NULL
      )
      if (!is.null(res) && nrow(res$classifications) > 0L) {
        per_type_pair[[pair_name]] <- res
      }
    }
  }

  # Collapsed: all 5' types x all 3' types
  collapsed <- tryCatch(
    classifyTopology(profiles, structures,
                      five_prime_types = five_types,
                      three_prime_types = three_types,
                      test = test, n_perm = n_perm),
    error = function(e) NULL
  )

  # Build comparison tibble
  comp_rows <- list()

  .extract_f2f <- function(result, label) {
    if (is.null(result) || nrow(result$summary) == 0L) {
      return(tibble::tibble(
        level = label, n_classified = 0L, F2F_n = 0L, F2F_pct = NA_real_,
        perm_p = NA_real_
      ))
    }
    n_total <- sum(result$summary$count)
    f2f_row <- result$summary[result$summary$topology == "F2F", ]
    f2f_n <- if (nrow(f2f_row) > 0L) f2f_row$count[1L] else 0L
    f2f_pct <- if (n_total > 0L) 100 * f2f_n / n_total else NA_real_
    perm_p <- if ("permutation_p_value" %in% names(result$summary) &&
                  nrow(f2f_row) > 0L) {
      f2f_row$permutation_p_value[1L]
    } else {
      NA_real_
    }
    tibble::tibble(level = label, n_classified = as.integer(n_total),
                    F2F_n = as.integer(f2f_n), F2F_pct = f2f_pct,
                    perm_p = perm_p)
  }

  for (nm in names(per_type_pair)) {
    comp_rows[[length(comp_rows) + 1L]] <- .extract_f2f(per_type_pair[[nm]], nm)
  }
  comp_rows[[length(comp_rows) + 1L]] <- .extract_f2f(collapsed, "collapsed")

  comparison <- dplyr::bind_rows(comp_rows)

  list(per_type_pair = per_type_pair, collapsed = collapsed,
       comparison = comparison)
}


#' Test Event Proximity (Clustering)
#'
#' Tests whether internal splicing events cluster closer together than
#' expected by chance using a permutation-based approach. Terminal events
#' (Alt_TSS, Alt_TES) are excluded. Only profiles with >= 2 internal events
#' contribute.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}
#'   with \code{isoform_id}, \code{tx_start}, \code{tx_end}.
#' @param n_perm Integer; number of permutations. Default 1000.
#' @return A list with: \code{z_score}, \code{p_value},
#'   \code{observed_mean_distance}, \code{expected_mean_distance},
#'   \code{n_profiles_used}, \code{n_profiles_excluded}.
#' @examples
#' data(example_profiles)
#' data(example_structures)
#' prox <- testProximity(example_profiles, example_structures, n_perm = 50L)
#' prox$z_score
#' @export
#' @importFrom stats sd
testProximity <- function(profiles, structures, n_perm = 1000L) {

  # Extract internal events with genomic positions
  events <- .extractEventsWithPositions(profiles, structures)
  internal <- events[!events$event_type %in% c("Alt_TSS", "Alt_TES"), ]

  if (nrow(internal) == 0L) {
    return(list(z_score = NA_real_, p_value = NA_real_,
                observed_mean_distance = NA_real_,
                expected_mean_distance = NA_real_,
                n_profiles_used = 0L, n_profiles_excluded = nrow(profiles)))
  }

  # Group by profile, keep only those with >= 2 internal events
  profile_key <- paste(internal$reference_isoform_id,
                       internal$comparator_isoform_id, sep = "|")
  profile_counts <- table(profile_key)
  qualifying <- names(profile_counts[profile_counts >= 2L])

  n_used <- length(qualifying)
  n_excluded <- nrow(profiles) - n_used

  if (n_used == 0L) {
    return(list(z_score = NA_real_, p_value = NA_real_,
                observed_mean_distance = NA_real_,
                expected_mean_distance = NA_real_,
                n_profiles_used = 0L, n_profiles_excluded = nrow(profiles)))
  }

  internal$profile_key <- profile_key
  internal <- internal[internal$profile_key %in% qualifying, ]

  # Compute observed mean pairwise distance per profile
  profile_keys_unique <- unique(internal$profile_key)
  obs_dists <- numeric(length(profile_keys_unique))
  n_events_vec <- integer(length(profile_keys_unique))
  gene_starts <- numeric(length(profile_keys_unique))
  gene_ends <- numeric(length(profile_keys_unique))

  for (j in seq_along(profile_keys_unique)) {
    pk <- profile_keys_unique[j]
    midpoints <- internal$event_midpoint[internal$profile_key == pk]
    obs_dists[j] <- mean(stats::dist(midpoints))
    n_events_vec[j] <- length(midpoints)
    gene_starts[j] <- internal$gene_start[internal$profile_key == pk][1]
    gene_ends[j] <- internal$gene_end[internal$profile_key == pk][1]
  }

  observed_mean <- mean(obs_dists)

  # Permutation test: randomly place events within gene span
  perm_means <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    perm_dists <- numeric(length(profile_keys_unique))
    for (j in seq_along(profile_keys_unique)) {
      pts <- stats::runif(n_events_vec[j], min = gene_starts[j],
                           max = gene_ends[j])
      perm_dists[j] <- mean(stats::dist(pts))
    }
    perm_means[p] <- mean(perm_dists)
  }

  perm_sd <- stats::sd(perm_means)
  perm_mean <- mean(perm_means)
  z_score <- if (perm_sd > 0) (observed_mean - perm_mean) / perm_sd else 0
  # One-sided: events cluster CLOSER than expected (observed < permuted)
  p_value <- mean(perm_means <= observed_mean)

  list(
    z_score = z_score,
    p_value = p_value,
    observed_mean_distance = observed_mean,
    expected_mean_distance = perm_mean,
    n_profiles_used = n_used,
    n_profiles_excluded = n_excluded
  )
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Extract events with genomic positions from profiles
#'
#' Unnests detailed_events and computes genomic midpoints and relative
#' positions within the gene span.
#'
#' @param profiles A tibble from buildProfiles().
#' @param structures A tibble from parseIsoformStructures().
#' @return A tibble with event_type, genomic_start, genomic_end,
#'   event_midpoint, relative_pos, gene_start, gene_end,
#'   reference_isoform_id, comparator_isoform_id.
#' @keywords internal
.extractEventsWithPositions <- function(profiles, structures) {

  # Build gene span lookup from structures via reference isoforms
  span_lookup <- structures[, c("isoform_id", "tx_start", "tx_end")]

  all_events <- list()
  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next

    ref_id <- profiles$reference_isoform_id[i]
    ref_span <- span_lookup[span_lookup$isoform_id == ref_id, ]
    if (nrow(ref_span) == 0L) next

    gene_start <- ref_span$tx_start[1]
    gene_end <- ref_span$tx_end[1]
    gene_length <- gene_end - gene_start

    if (gene_length <= 0L) next

    genomic_start <- pmin(events$five_prime, events$three_prime)
    genomic_end <- pmax(events$five_prime, events$three_prime)
    midpoint <- (genomic_start + genomic_end) / 2
    relative_pos <- pmax(0, pmin(1, (midpoint - gene_start) / gene_length))

    all_events[[length(all_events) + 1L]] <- tibble::tibble(
      event_type = events$event_type,
      genomic_start = genomic_start,
      genomic_end = genomic_end,
      event_midpoint = midpoint,
      relative_pos = relative_pos,
      gene_start = gene_start,
      gene_end = gene_end,
      gene_id = profiles$gene_id[i],
      reference_isoform_id = ref_id,
      comparator_isoform_id = profiles$comparator_isoform_id[i]
    )
  }

  if (length(all_events) == 0L) {
    return(tibble::tibble(
      event_type = character(0), genomic_start = numeric(0),
      genomic_end = numeric(0), event_midpoint = numeric(0),
      relative_pos = numeric(0), gene_start = numeric(0),
      gene_end = numeric(0), gene_id = character(0),
      reference_isoform_id = character(0),
      comparator_isoform_id = character(0)
    ))
  }

  dplyr::bind_rows(all_events)
}


#' Classify topology of one 5'-boundary x 3'-boundary event pair
#'
#' @param a5ss_event One-row data frame for a 5' boundary event (A5SS or Partial_IR_5).
#' @param a3ss_event One-row data frame for a 3' boundary event (A3SS or Partial_IR_3).
#' @param exon_starts Integer vector of reference exon start coordinates.
#' @param exon_ends Integer vector of reference exon end coordinates.
#' @param gene_strand Character; "+" or "-".
#' @return Character: "F2F", "B2B", or "Distal".
#' @keywords internal
.classifyEventPairTopology <- function(a5ss_event, a3ss_event,
                                        exon_starts, exon_ends,
                                        gene_strand) {
  # Event midpoints in genomic space
  a5ss_mid <- (pmin(a5ss_event$five_prime, a5ss_event$three_prime) +
                pmax(a5ss_event$five_prime, a5ss_event$three_prime)) / 2
  a3ss_mid <- (pmin(a3ss_event$five_prime, a3ss_event$three_prime) +
                pmax(a3ss_event$five_prime, a3ss_event$three_prime)) / 2

  # Find nearest exon index for each event
  # A5SS affects the 3' end of an exon (+ strand) or 5' end (- strand)
  # A3SS affects the 5' start of an exon (+ strand) or 3' end (- strand)
  if (gene_strand == "+") {
    a5ss_exon_idx <- which.min(abs(exon_ends - a5ss_mid))
    a3ss_exon_idx <- which.min(abs(exon_starts - a3ss_mid))
  } else {
    a5ss_exon_idx <- which.min(abs(exon_starts - a5ss_mid))
    a3ss_exon_idx <- which.min(abs(exon_ends - a3ss_mid))
  }

  if (a5ss_exon_idx == a3ss_exon_idx) {
    "B2B"
  } else if (gene_strand == "+" && a5ss_exon_idx == a3ss_exon_idx - 1L) {
    "F2F"
  } else if (gene_strand == "-" && a5ss_exon_idx == a3ss_exon_idx + 1L) {
    "F2F"
  } else {
    "Distal"
  }
}
