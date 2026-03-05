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
#' For profiles containing both A5SS and A3SS events, classifies each
#' A5SS x A3SS pair into topology categories:
#' \itemize{
#'   \item \strong{F2F} (Face-to-Face): events flank the same intron
#'   \item \strong{B2B} (Back-to-Back): events are on the same exon
#'   \item \strong{Distal}: events are on non-adjacent structural elements
#' }
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
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
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
classifyTopology <- function(profiles, structures,
                              test = TRUE, n_perm = 1000L) {

  all_classifications <- list()
  perm_data <- list()

  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next

    a5ss <- events[events$event_type == "A5SS", , drop = FALSE]
    a3ss <- events[events$event_type == "A3SS", , drop = FALSE]

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
            event_a_type = "A5SS",
            event_b_type = "A3SS",
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


#' Classify topology of one A5SS x A3SS event pair
#'
#' @param a5ss_event One-row data frame for an A5SS event.
#' @param a3ss_event One-row data frame for an A3SS event.
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
