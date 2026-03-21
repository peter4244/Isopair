#' Compute Premature Termination Codon (PTC) Status
#'
#' Determines PTC status for coding isoforms using the 50-nucleotide rule:
#' a stop codon is premature if it is >50 nt upstream of the last
#' exon-exon junction (EJC) in mRNA space.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}
#'   with columns: \code{isoform_id}, \code{gene_id}, \code{strand},
#'   \code{n_exons}, \code{exon_starts} (list), \code{exon_ends} (list).
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with columns: \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{orf_length}.
#' @return A tibble with columns: \code{isoform_id}, \code{gene_id},
#'   \code{n_exons}, \code{orf_length}, \code{ptc_distance},
#'   \code{has_ptc}, \code{n_downstream_ejcs}, \code{stop_in_last_exon}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' ptc <- computePtcStatus(structures, cds)
#' head(ptc)
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr inner_join filter
#' @importFrom rlang .data
computePtcStatus <- function(structures, cds_metadata) {

  # Filter to coding isoforms present in both inputs
  coding <- dplyr::filter(cds_metadata, .data$coding_status == "coding")

  if (nrow(coding) == 0L) {
    return(tibble::tibble(
      isoform_id = character(0), gene_id = character(0),
      n_exons = integer(0), orf_length = integer(0),
      ptc_distance = numeric(0), has_ptc = logical(0),
      n_downstream_ejcs = integer(0), stop_in_last_exon = logical(0)
    ))
  }

  merged <- dplyr::inner_join(
    structures[, c("isoform_id", "gene_id", "strand", "n_exons",
                    "exon_starts", "exon_ends")],
    coding[, c("isoform_id", "cds_start", "cds_stop", "orf_length")],
    by = "isoform_id"
  )

  if (nrow(merged) == 0L) {
    return(tibble::tibble(
      isoform_id = character(0), gene_id = character(0),
      n_exons = integer(0), orf_length = integer(0),
      ptc_distance = numeric(0), has_ptc = logical(0),
      n_downstream_ejcs = integer(0), stop_in_last_exon = logical(0)
    ))
  }

  results <- vector("list", nrow(merged))

  for (i in seq_len(nrow(merged))) {
    row <- merged[i, ]
    feat <- .computePtcFeatures(
      cds_start = row$cds_start,
      cds_stop = row$cds_stop,
      strand = row$strand,
      exon_starts = row$exon_starts[[1]],
      exon_ends = row$exon_ends[[1]],
      n_exons = row$n_exons
    )

    results[[i]] <- tibble::tibble(
      isoform_id = row$isoform_id,
      gene_id = row$gene_id,
      n_exons = row$n_exons,
      orf_length = row$orf_length,
      ptc_distance = feat$ptc_distance,
      has_ptc = feat$has_ptc,
      n_downstream_ejcs = feat$n_downstream_ejcs,
      stop_in_last_exon = feat$stop_in_last_exon
    )
  }

  dplyr::bind_rows(results)
}


#' Compute PTC Features for a Single Isoform
#'
#' Internal function that calculates PTC distance, downstream EJC count,
#' and stop codon position relative to exon structure.
#'
#' @param cds_start Integer; min genomic CDS coordinate.
#' @param cds_stop Integer; max genomic CDS coordinate.
#' @param strand Character; "+" or "-".
#' @param exon_starts Integer vector of exon start coordinates.
#' @param exon_ends Integer vector of exon end coordinates.
#' @param n_exons Integer; number of exons.
#' @return A list with: \code{ptc_distance}, \code{has_ptc},
#'   \code{n_downstream_ejcs}, \code{stop_in_last_exon}.
#' @keywords internal
.computePtcFeatures <- function(cds_start, cds_stop, strand,
                                exon_starts, exon_ends, n_exons) {

  # Single-exon: no EJCs possible
  if (n_exons == 1L) {
    return(list(
      ptc_distance = NA_real_,
      has_ptc = FALSE,
      n_downstream_ejcs = 0L,
      stop_in_last_exon = TRUE
    ))
  }

  # Biological stop codon position
  stop_pos <- if (strand == "+") cds_stop else cds_start

  # Order exons 5'→3' biologically
  if (strand == "+") {
    ord <- order(exon_starts)
  } else {
    ord <- order(exon_starts, decreasing = TRUE)
  }
  e_starts <- exon_starts[ord]
  e_ends <- exon_ends[ord]

  # Exon lengths and cumulative mRNA positions
  exon_lengths <- e_ends - e_starts + 1L
  cum_starts <- c(0L, cumsum(exon_lengths[-length(exon_lengths)]))

  # Last EJC position in mRNA space
  last_ejc_mRNA <- sum(exon_lengths[-length(exon_lengths)])

  # Find which exon contains the stop codon
  stop_exon_idx <- NA_integer_
  for (j in seq_along(e_starts)) {
    if (stop_pos >= e_starts[j] && stop_pos <= e_ends[j]) {
      stop_exon_idx <- j
      break
    }
  }

  # Stop codon not found in any exon
  if (is.na(stop_exon_idx)) {
    return(list(
      ptc_distance = NA_real_,
      has_ptc = NA,
      n_downstream_ejcs = NA_integer_,
      stop_in_last_exon = NA
    ))
  }

  # Compute offset within the stop codon's exon (strand-aware)
  if (strand == "+") {
    offset <- stop_pos - e_starts[stop_exon_idx]
  } else {
    offset <- e_ends[stop_exon_idx] - stop_pos
  }

  # mRNA position of stop codon
  stop_mRNA <- cum_starts[stop_exon_idx] + offset

  # PTC distance: positive means stop is upstream of last EJC
  ptc_distance <- last_ejc_mRNA - stop_mRNA
  has_ptc <- ptc_distance > 50

  # Downstream EJCs and last-exon check
  n_downstream_ejcs <- as.integer(n_exons - stop_exon_idx)
  stop_in_last_exon <- (stop_exon_idx == n_exons)

  list(
    ptc_distance = as.numeric(ptc_distance),
    has_ptc = has_ptc,
    n_downstream_ejcs = n_downstream_ejcs,
    stop_in_last_exon = stop_in_last_exon
  )
}


#' Test PTC Association Between Pair Members
#'
#' Performs a paired McNemar test on PTC status discordance between
#' reference and comparator isoforms. Optionally stratifies by a grouping
#' variable using Fisher's exact test.
#'
#' @param ptc_status A tibble from \code{\link{computePtcStatus}()} with
#'   at least \code{isoform_id} and \code{has_ptc} columns.
#' @param pairs A data frame with \code{reference_isoform_id} and
#'   \code{comparator_isoform_id} columns.
#' @param group_col Optional character; column name in \code{pairs} to
#'   stratify Fisher's tests by.
#' @return A list with: \code{pair_summary} (tibble of per-pair PTC
#'   status), \code{ref_ptc_rate}, \code{comp_ptc_rate},
#'   \code{mcnemar_test} (htest object or NULL), and optionally
#'   \code{group_tests} (list of Fisher test results per group).
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' ptc <- computePtcStatus(structures, cds)
#' data(example_pairs)
#' result <- testPtcAssociation(ptc, example_pairs)
#' result$ref_ptc_rate
#' @export
#' @importFrom stats mcnemar.test fisher.test
testPtcAssociation <- function(ptc_status, pairs, group_col = NULL) {

  ptc_lookup <- ptc_status[, c("isoform_id", "has_ptc")]

  # Join PTC status for both pair members
  pair_ptc <- pairs
  ref_idx <- match(pair_ptc$reference_isoform_id, ptc_lookup$isoform_id)
  comp_idx <- match(pair_ptc$comparator_isoform_id, ptc_lookup$isoform_id)

  pair_ptc$ref_has_ptc <- ptc_lookup$has_ptc[ref_idx]
  pair_ptc$comp_has_ptc <- ptc_lookup$has_ptc[comp_idx]

  # Filter to pairs where both have non-NA PTC status
  pair_ptc <- pair_ptc[!is.na(pair_ptc$ref_has_ptc) &
                         !is.na(pair_ptc$comp_has_ptc), ]

  ref_ptc_rate <- if (nrow(pair_ptc) > 0L) {
    mean(pair_ptc$ref_has_ptc)
  } else NA_real_

  comp_ptc_rate <- if (nrow(pair_ptc) > 0L) {
    mean(pair_ptc$comp_has_ptc)
  } else NA_real_

  # McNemar test on discordant pairs
  mcnemar_result <- NULL
  if (nrow(pair_ptc) >= 2L) {
    # Build 2x2 contingency table
    tbl <- table(
      ref_ptc = factor(pair_ptc$ref_has_ptc, levels = c(FALSE, TRUE)),
      comp_ptc = factor(pair_ptc$comp_has_ptc, levels = c(FALSE, TRUE))
    )
    # McNemar requires discordant cells > 0
    if (tbl[1, 2] + tbl[2, 1] > 0L) {
      mcnemar_result <- tryCatch(
        stats::mcnemar.test(tbl),
        error = function(e) NULL
      )
    }
  }

  result <- list(
    pair_summary = tibble::as_tibble(pair_ptc),
    ref_ptc_rate = ref_ptc_rate,
    comp_ptc_rate = comp_ptc_rate,
    mcnemar_test = mcnemar_result
  )

  # Optional group-stratified Fisher tests
  if (!is.null(group_col) && group_col %in% names(pair_ptc) &&
      nrow(pair_ptc) > 0L) {
    groups <- unique(pair_ptc[[group_col]])
    group_tests <- lapply(groups, function(g) {
      subset_data <- pair_ptc[pair_ptc[[group_col]] == g, ]
      if (nrow(subset_data) < 2L) return(NULL)

      tbl <- table(
        ref_ptc = factor(subset_data$ref_has_ptc, levels = c(FALSE, TRUE)),
        comp_ptc = factor(subset_data$comp_has_ptc, levels = c(FALSE, TRUE))
      )
      tryCatch(
        stats::fisher.test(tbl),
        error = function(e) NULL
      )
    })
    names(group_tests) <- groups
    result$group_tests <- group_tests
  }

  result
}


#' Analyze Frame Disruption in Splicing Events
#'
#' @description
#' \strong{Deprecated.} Use \code{\link{analyzeFrameWalk}()} instead.
#'
#' \code{analyzeFrameDisruption()} is deprecated in favour of
#' \code{\link{analyzeFrameWalk}()}, which provides richer per-event output
#' (cumulative frame offsets, compensatory event detection, frameshifted CDS
#' extent). This function is now a thin wrapper that calls
#' \code{analyzeFrameWalk()} internally and transforms the output to match
#' the original return format.
#'
#' For each splicing event in a set of profiles, computes CDS overlap and
#' determines whether the event disrupts the reading frame (CDS overlap
#' not divisible by 3). Only frameshift-capable event types are evaluated:
#' SE, Missing_Internal, A5SS, A3SS, Partial_IR_5, Partial_IR_3.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} containing
#'   a \code{detailed_events} list column.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with columns: \code{isoform_id}, \code{coding_status},
#'   \code{cds_start}, \code{cds_stop}, \code{orf_length}, \code{strand}.
#' @return A list with two elements:
#'   \describe{
#'     \item{events}{Tibble of all events with added columns:
#'       \code{cds_overlap_bp}, \code{affects_cds}, \code{frame_disrupting}}
#'     \item{profile_summary}{Tibble with per-pair summary: \code{gene_id},
#'       \code{reference_isoform_id}, \code{comparator_isoform_id},
#'       \code{n_cds_events}, \code{n_frame_disrupting},
#'       \code{has_frame_disruption}}
#'   }
#' @seealso \code{\link{analyzeFrameWalk}} for the replacement function with
#'   richer output.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' cds <- extractCdsAnnotations(gtf_path, verbose = FALSE)
#' data(example_profiles)
#' fd <- analyzeFrameDisruption(example_profiles, cds)
#' fd$profile_summary
#' @export
#' @importFrom dplyr bind_rows left_join
#' @importFrom tibble tibble
analyzeFrameDisruption <- function(profiles, cds_metadata) {

  .Deprecated("analyzeFrameWalk")

  frameshift_types <- c("SE", "Missing_Internal", "A5SS", "A3SS",
                         "Partial_IR_5", "Partial_IR_3")

  # Run the full frame walk analysis
  fw <- analyzeFrameWalk(profiles, cds_metadata)

  # Build CDS lookup for event-level annotation (CDS overlap for ALL event
  # types, including non-frameshift types like IR that the walk ignores)
  cds_lookup <- cds_metadata[cds_metadata$coding_status == "coding",
                              c("isoform_id", "cds_start", "cds_stop")]

  all_events <- list()
  profile_summaries <- list()

  for (i in seq_len(nrow(profiles))) {
    gene <- profiles$gene_id[i]
    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]
    events <- profiles$detailed_events[[i]]

    if (is.null(events) || nrow(events) == 0L) {
      profile_summaries[[i]] <- tibble::tibble(
        gene_id = gene,
        reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        n_cds_events = 0L,
        n_frame_disrupting = 0L,
        has_frame_disruption = FALSE
      )
      next
    }

    # Look up CDS for reference isoform
    ref_cds <- cds_lookup[cds_lookup$isoform_id == ref_id, ]
    has_ref_cds <- nrow(ref_cds) > 0L

    # Compute CDS overlap for each event (all types, not just frameshift)
    evt_results <- events
    n_evt <- nrow(evt_results)

    cds_overlap_bp <- integer(n_evt)
    affects_cds <- logical(n_evt)
    frame_disrupting <- logical(n_evt)
    splice_group <- integer(n_evt)

    if (has_ref_cds) {
      ref_cds_start <- ref_cds$cds_start[1]
      ref_cds_stop <- ref_cds$cds_stop[1]

      for (j in seq_len(n_evt)) {
        event_start <- pmin(evt_results$five_prime[j],
                            evt_results$three_prime[j])
        event_end <- pmax(evt_results$five_prime[j],
                          evt_results$three_prime[j])

        overlap_start <- pmax(event_start, ref_cds_start)
        overlap_end <- pmin(event_end, ref_cds_stop)
        overlap_bp <- pmax(0L, as.integer(overlap_end - overlap_start + 1L))
        cds_overlap_bp[j] <- overlap_bp
        affects_cds[j] <- overlap_bp > 0L
      }

      # Derive splice_group and frame_disrupting from walk events
      fw_ev <- fw$events
      if (nrow(fw_ev) > 0L &&
          "reference_isoform_id" %in% names(fw_ev)) {
        fw_events_i <- fw_ev[
          fw_ev$reference_isoform_id == ref_id &
          fw_ev$comparator_isoform_id == comp_id, , drop = FALSE
        ]
      } else {
        fw_events_i <- fw_ev[integer(0), , drop = FALSE]
      }

      if (nrow(fw_events_i) > 0L) {
        # Determine which splice groups are frame-disrupting using the old
        # per-group semantics: a group is disrupting if its net signed CDS
        # change is not divisible by 3.
        for (sg in unique(fw_events_i$splice_group)) {
          sg_rows <- fw_events_i[fw_events_i$splice_group == sg, , drop = FALSE]
          group_net_change <- sum(sg_rows$signed_cds_change)
          group_disrupts <- (group_net_change %% 3L) != 0L

          # Map walk events back to original events by matching event_type,
          # direction, and genomic coordinates
          for (r in seq_len(nrow(sg_rows))) {
            match_idx <- which(
              evt_results$event_type == sg_rows$event_type[r] &
              evt_results$direction == sg_rows$direction[r] &
              pmin(evt_results$five_prime, evt_results$three_prime) ==
                sg_rows$genomic_start[r] &
              pmax(evt_results$five_prime, evt_results$three_prime) ==
                sg_rows$genomic_end[r]
            )
            if (length(match_idx) > 0L) {
              splice_group[match_idx[1]] <- sg
              if (group_disrupts) frame_disrupting[match_idx[1]] <- TRUE
            }
          }
        }
      }
    }

    evt_results$cds_overlap_bp <- cds_overlap_bp
    evt_results$affects_cds <- affects_cds
    evt_results$frame_disrupting <- frame_disrupting
    evt_results$splice_group <- splice_group

    all_events[[i]] <- evt_results

    # Profile-level summary: count disrupted splice groups
    fs_events <- evt_results[evt_results$event_type %in% frameshift_types, ]
    n_cds <- sum(fs_events$affects_cds)
    disrupted_groups <- unique(fs_events$splice_group[fs_events$frame_disrupting])
    disrupted_groups <- disrupted_groups[disrupted_groups > 0L]
    n_fd <- length(disrupted_groups)

    profile_summaries[[i]] <- tibble::tibble(
      gene_id = gene,
      reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      n_cds_events = n_cds,
      n_frame_disrupting = n_fd,
      has_frame_disruption = n_fd > 0L
    )
  }

  list(
    events = dplyr::bind_rows(all_events),
    profile_summary = dplyr::bind_rows(profile_summaries)
  )
}


#' Analyze Frame Walk Per Profile
#'
#' Walks events 5'->3' per profile, tracking cumulative reading frame offset.
#' Detects frameshifts, compensatory events that restore frame, and computes
#' the extent of CDS under a shifted reading frame.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with columns: \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{orf_length}, \code{strand}.
#' @param ptc_table Optional tibble with \code{isoform_id} and
#'   \code{has_ptc} columns (e.g., from \code{\link{computePtcStatus}()}).
#'   If provided, comparator PTC status is cross-referenced.
#' @return A list with two elements:
#'   \describe{
#'     \item{events}{Tibble — one row per CDS-affecting event with columns:
#'       \code{gene_id}, \code{reference_isoform_id},
#'       \code{comparator_isoform_id}, \code{event_type}, \code{direction},
#'       \code{genomic_start}, \code{genomic_end}, \code{cds_overlap_bp},
#'       \code{signed_cds_change}, \code{cumulative_frame_offset},
#'       \code{is_frameshift}, \code{is_compensatory}.}
#'     \item{profile_summary}{Tibble — one row per profile: \code{gene_id},
#'       \code{reference_isoform_id}, \code{comparator_isoform_id},
#'       \code{n_cds_events}, \code{n_splice_groups}, \code{n_frameshifts},
#'       \code{n_compensatory},
#'       \code{frame_resolved}, \code{total_frameshifted_cds_bp} (capped at
#'       \code{orf_length}; upper bound for multi-exon CDS),
#'       \code{pct_orf_frameshifted}, \code{comparator_has_ptc}.}
#'   }
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' cds <- extractCdsAnnotations(gtf_path, verbose = FALSE)
#' data(example_profiles)
#' fw <- analyzeFrameWalk(example_profiles, cds)
#' fw$profile_summary
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
analyzeFrameWalk <- function(profiles, cds_metadata, ptc_table = NULL) {

  frameshift_types <- c("SE", "Missing_Internal", "A5SS", "A3SS",
                         "Partial_IR_5", "Partial_IR_3")

  # CDS lookup
  cds_lookup <- cds_metadata[cds_metadata$coding_status == "coding",
                              c("isoform_id", "cds_start", "cds_stop",
                                "orf_length", "strand")]

  # PTC lookup
  ptc_lookup <- NULL
  if (!is.null(ptc_table) && nrow(ptc_table) > 0L) {
    ptc_lookup <- ptc_table[, c("isoform_id", "has_ptc")]
  }

  all_events <- list()
  profile_summaries <- list()

  for (i in seq_len(nrow(profiles))) {
    gene <- profiles$gene_id[i]
    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]
    events <- profiles$detailed_events[[i]]

    # Default summary for non-coding or no events
    ref_cds <- cds_lookup[cds_lookup$isoform_id == ref_id, ]

    # PTC cross-reference
    comp_ptc <- NA
    if (!is.null(ptc_lookup)) {
      ptc_match <- ptc_lookup[ptc_lookup$isoform_id == comp_id, ]
      if (nrow(ptc_match) > 0L) comp_ptc <- ptc_match$has_ptc[1L]
    }

    if (nrow(ref_cds) == 0L || is.null(events) || nrow(events) == 0L) {
      profile_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        n_cds_events = 0L, n_splice_groups = 0L,
        n_frameshifts = 0L, n_compensatory = 0L,
        frame_resolved = NA, total_frameshifted_cds_bp = 0L,
        pct_orf_frameshifted = NA_real_, comparator_has_ptc = comp_ptc
      )
      next
    }

    ref_cds_start <- ref_cds$cds_start[1L]
    ref_cds_stop <- ref_cds$cds_stop[1L]
    orf_length <- ref_cds$orf_length[1L]
    gene_strand <- ref_cds$strand[1L]

    # Filter to frameshift-capable events with CDS overlap
    fs_events <- events[events$event_type %in% frameshift_types, , drop = FALSE]
    if (nrow(fs_events) == 0L) {
      profile_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        n_cds_events = 0L, n_splice_groups = 0L,
        n_frameshifts = 0L, n_compensatory = 0L,
        frame_resolved = NA, total_frameshifted_cds_bp = 0L,
        pct_orf_frameshifted = NA_real_, comparator_has_ptc = comp_ptc
      )
      next
    }

    # Compute CDS overlap and genomic coords
    genomic_starts <- pmin(fs_events$five_prime, fs_events$three_prime)
    genomic_ends <- pmax(fs_events$five_prime, fs_events$three_prime)
    overlap_starts <- pmax(genomic_starts, ref_cds_start)
    overlap_ends <- pmin(genomic_ends, ref_cds_stop)
    cds_overlap_bp <- as.integer(pmax(0L, overlap_ends - overlap_starts + 1L))

    # Keep only events with CDS overlap
    has_cds <- cds_overlap_bp > 0L
    if (sum(has_cds) == 0L) {
      profile_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        n_cds_events = 0L, n_splice_groups = 0L,
        n_frameshifts = 0L, n_compensatory = 0L,
        frame_resolved = NA, total_frameshifted_cds_bp = 0L,
        pct_orf_frameshifted = NA_real_, comparator_has_ptc = comp_ptc
      )
      next
    }

    cds_events <- fs_events[has_cds, , drop = FALSE]
    cds_overlap <- cds_overlap_bp[has_cds]
    g_starts <- genomic_starts[has_cds]
    g_ends <- genomic_ends[has_cds]

    # Sort events 5'→3': ascending for + strand, descending for - strand
    if (!is.na(gene_strand) && gene_strand == "-") {
      sort_order <- order(g_starts, decreasing = TRUE)
    } else {
      sort_order <- order(g_starts)
    }
    cds_events <- cds_events[sort_order, , drop = FALSE]
    cds_overlap <- cds_overlap[sort_order]
    g_starts <- g_starts[sort_order]
    g_ends <- g_ends[sort_order]

    n_cds <- length(cds_overlap)

    # --- Splice group assignment ---
    # Events sharing a junction (ref or comp) arise from the same splicing
    # decision and should be evaluated as a unit for frameshift assessment.
    # Use union-find to identify connected components.
    sg_parent <- seq_len(n_cds)
    sg_find <- function(x) {
      while (sg_parent[x] != x) {
        sg_parent[x] <<- sg_parent[sg_parent[x]]
        x <- sg_parent[x]
      }
      x
    }
    sg_union <- function(a, b) {
      ra <- sg_find(a); rb <- sg_find(b)
      if (ra != rb) sg_parent[ra] <<- rb
    }

    # Parse junctions for each event (required for splice group assignment)
    if (is.null(cds_events$ref_junctions) || is.null(cds_events$comp_junctions)) {
      stop("Events must include 'ref_junctions' and 'comp_junctions' columns. ",
           "These are produced by detectEvents() and required for splice group ",
           "assignment in frameshift analysis.", call. = FALSE)
    }
    parse_jxns <- function(jxn_str) {
      if (is.na(jxn_str) || jxn_str == "") return(character(0))
      unlist(strsplit(jxn_str, ","))
    }
    ref_jxn_list <- lapply(cds_events$ref_junctions, parse_jxns)
    comp_jxn_list <- lapply(cds_events$comp_junctions, parse_jxns)

    # Build junction → event index map, then union events sharing junctions
    jxn_map <- list()
    for (k in seq_len(n_cds)) {
      all_jxns <- unique(c(ref_jxn_list[[k]], comp_jxn_list[[k]]))
      for (jxn in all_jxns) {
        if (jxn %in% names(jxn_map)) {
          sg_union(k, jxn_map[[jxn]])
        }
        jxn_map[[jxn]] <- k
      }
    }

    # Resolve group IDs
    splice_groups <- integer(n_cds)
    for (k in seq_len(n_cds)) splice_groups[k] <- sg_find(k)
    # Renumber to sequential 1..n_groups in 5'→3' order
    unique_groups <- unique(splice_groups)
    group_remap <- setNames(seq_along(unique_groups), unique_groups)
    splice_groups <- as.integer(group_remap[as.character(splice_groups)])

    # --- Frame walk at group level ---
    n_groups <- max(splice_groups)
    signed_changes <- integer(n_cds)
    for (k in seq_len(n_cds)) {
      signed_changes[k] <- if (cds_events$direction[k] == "GAIN") {
        cds_overlap[k]
      } else {
        -cds_overlap[k]
      }
    }

    # Compute per-group net signed change
    group_net_change <- integer(n_groups)
    group_last_idx <- integer(n_groups)  # last event index in each group
    group_first_idx <- integer(n_groups)
    group_first_idx[] <- n_cds + 1L
    for (k in seq_len(n_cds)) {
      g <- splice_groups[k]
      group_net_change[g] <- group_net_change[g] + signed_changes[k]
      if (k > group_last_idx[g]) group_last_idx[g] <- k
      if (k < group_first_idx[g]) group_first_idx[g] <- k
    }

    cum_offsets <- integer(n_cds)
    is_fs <- logical(n_cds)
    is_comp <- logical(n_cds)

    cumulative <- 0L
    n_frameshifts <- 0L
    n_compensatory <- 0L
    last_shift_pos <- NA_real_
    frameshifted_bp <- 0L
    processed_groups <- logical(n_groups)

    for (k in seq_len(n_cds)) {
      # Track per-event cumulative offset for transparency
      cumulative <- cumulative + signed_changes[k]
      cum_offsets[k] <- ((cumulative %% 3L) + 3L) %% 3L

      g <- splice_groups[k]
      if (k < group_last_idx[g]) next  # not last event in group yet

      # Evaluate frameshift at group boundary
      net_change <- group_net_change[g]
      prev_group_offset <- (((cumulative - net_change) %% 3L) + 3L) %% 3L
      new_group_offset <- cum_offsets[k]

      if (prev_group_offset == 0L && new_group_offset != 0L) {
        is_fs[k] <- TRUE
        n_frameshifts <- n_frameshifts + 1L
        # 5' boundary of the group's genomic span
        first_k <- group_first_idx[g]
        last_shift_pos <- if (!is.na(gene_strand) && gene_strand == "-") {
          g_starts[first_k]
        } else {
          g_ends[k]
        }
      }

      if (prev_group_offset != 0L && new_group_offset == 0L) {
        is_comp[k] <- TRUE
        n_compensatory <- n_compensatory + 1L
        if (!is.na(last_shift_pos)) {
          if (!is.na(gene_strand) && gene_strand == "-") {
            shift_start <- min(last_shift_pos, g_ends[k])
            shift_end <- max(last_shift_pos, g_ends[k])
          } else {
            shift_start <- min(last_shift_pos, g_starts[k])
            shift_end <- max(last_shift_pos, g_starts[k])
          }
          bp_start <- max(shift_start, ref_cds_start)
          bp_end <- min(shift_end, ref_cds_stop)
          frameshifted_bp <- frameshifted_bp +
            as.integer(max(0L, bp_end - bp_start + 1L))
        }
        last_shift_pos <- NA_real_
      }
    }

    # Unresolved frameshift: remaining CDS from last_shift_pos to CDS end
    final_offset <- ((cumulative %% 3L) + 3L) %% 3L
    frame_resolved <- (final_offset == 0L)

    if (!frame_resolved && !is.na(last_shift_pos)) {
      if (!is.na(gene_strand) && gene_strand == "-") {
        bp_start <- ref_cds_start
        bp_end <- min(last_shift_pos, ref_cds_stop)
      } else {
        bp_start <- max(last_shift_pos, ref_cds_start)
        bp_end <- ref_cds_stop
      }
      frameshifted_bp <- frameshifted_bp +
        as.integer(max(0L, bp_end - bp_start + 1L))
    }

    # Cap frameshifted_bp at orf_length: the genomic interval calculation
    # may overcount for multi-exon CDS (includes intervening introns).
    # The capped value is a conservative upper bound.
    frameshifted_bp <- min(frameshifted_bp, orf_length)

    pct_orf <- if (orf_length > 0L) {
      100 * frameshifted_bp / orf_length
    } else NA_real_

    # Build event-level output
    all_events[[i]] <- tibble::tibble(
      gene_id = gene, reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      event_type = cds_events$event_type,
      direction = cds_events$direction,
      genomic_start = as.integer(g_starts),
      genomic_end = as.integer(g_ends),
      cds_overlap_bp = cds_overlap,
      signed_cds_change = signed_changes,
      cumulative_frame_offset = cum_offsets,
      splice_group = splice_groups,
      is_frameshift = is_fs,
      is_compensatory = is_comp
    )

    profile_summaries[[i]] <- tibble::tibble(
      gene_id = gene, reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      n_cds_events = n_cds,
      n_splice_groups = n_groups,
      n_frameshifts = n_frameshifts,
      n_compensatory = n_compensatory,
      frame_resolved = frame_resolved,
      total_frameshifted_cds_bp = frameshifted_bp,
      pct_orf_frameshifted = pct_orf,
      comparator_has_ptc = comp_ptc
    )
  }

  list(
    events = dplyr::bind_rows(all_events),
    profile_summary = dplyr::bind_rows(profile_summaries)
  )
}


# =============================================================================
# Cross-Isoform Reading Frame Comparison
# =============================================================================

#' Get Start Codon Position
#'
#' Returns the genomic position of the 5'-most CDS nucleotide for an isoform.
#' @param exons Data frame with cds_exon_start, cds_exon_end for one isoform.
#' @param strand "+" or "-".
#' @return Integer genomic position.
#' @keywords internal
.getStartCodonPos <- function(exons, strand) {
  if (strand == "-") {
    max(exons$cds_exon_end)
  } else {
    min(exons$cds_exon_start)
  }
}


#' Compute CDS Phase at a Genomic Position
#'
#' Walks through an isoform's CDS exons in 5'-to-3' order, accumulating
#' coding bp, and returns the reading frame phase (0/1/2) at a specific
#' genomic position.
#'
#' @param exons Data frame with cds_exon_start, cds_exon_end for one isoform,
#'   already sorted genomically (ascending start).
#' @param strand "+" or "-".
#' @param pos Integer genomic position (must be within a CDS exon).
#' @return Integer phase (0, 1, or 2).
#' @keywords internal
.computePhaseAtPosition <- function(exons, strand, pos) {
  n <- nrow(exons)
  if (n == 0L) return(NA_integer_)

  # Walk 5'→3'
  if (strand == "-") {
    # Descending genomic order (highest end first)
    ord <- order(exons$cds_exon_end, decreasing = TRUE)
  } else {
    ord <- order(exons$cds_exon_start)
  }

  cumulative <- 0L
  for (idx in ord) {
    es <- exons$cds_exon_start[idx]
    ee <- exons$cds_exon_end[idx]

    if (pos >= es && pos <= ee) {
      # Position is in this exon — add partial contribution
      if (strand == "-") {
        cumulative <- cumulative + as.integer(ee - pos)
      } else {
        cumulative <- cumulative + as.integer(pos - es)
      }
      return(as.integer(cumulative %% 3L))
    }

    # Full exon — accumulate all bp
    cumulative <- cumulative + as.integer(ee - es + 1L)
  }

  # Position not found in any CDS exon
  NA_integer_
}


#' Intersect Two Sets of Genomic Intervals
#'
#' Computes the genomic intersection of two sets of genomic intervals.
#' Both inputs are assumed to be non-overlapping within each set.
#'
#' @param a Data frame with columns \code{start} and \code{end} (integer).
#' @param b Data frame with columns \code{start} and \code{end} (integer).
#' @return Data frame with columns \code{start} and \code{end} for the
#'   intersection regions.
#' @keywords internal
.intersectCdsIntervals <- function(a, b) {
  if (nrow(a) == 0L || nrow(b) == 0L) {
    return(data.frame(start = integer(0), end = integer(0)))
  }

  # Sort both by start
  a <- a[order(a$start), , drop = FALSE]
  b <- b[order(b$start), , drop = FALSE]

  result_start <- integer(0)
  result_end <- integer(0)
  i <- 1L
  j <- 1L

  while (i <= nrow(a) && j <= nrow(b)) {
    # Compute overlap
    ov_start <- max(a$start[i], b$start[j])
    ov_end <- min(a$end[i], b$end[j])

    if (ov_start <= ov_end) {
      result_start <- c(result_start, ov_start)
      result_end <- c(result_end, ov_end)
    }

    # Advance the interval that ends first (both if tied)
    a_end <- a$end[i]
    b_end <- b$end[j]
    if (a_end <= b_end) i <- i + 1L
    if (b_end <= a_end) j <- j + 1L
  }

  data.frame(start = result_start, end = result_end)
}


#' Compare Reading Frames Between Isoform Pairs
#'
#' For each isoform pair, determines whether the reference and comparator
#' translate their shared CDS regions in the same reading frame. Classifies
#' pairs by start codon relationship (same vs different) and frame agreement
#' at shared CDS positions.
#'
#' This function complements \code{\link{analyzeFrameWalk}()}, which tracks
#' cumulative frame offsets through events using only the reference CDS.
#' \code{compareIsoformFrames()} uses both isoforms' CDS annotations
#' independently, correctly handling pairs with different start codons.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with columns
#'   \code{gene_id}, \code{reference_isoform_id}, \code{comparator_isoform_id}.
#' @param cds_exons A tibble from \code{\link{extractCdsExons}()} with
#'   per-exon CDS intervals. Must come from the same GTF as \code{cds_metadata}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()} with
#'   \code{isoform_id} and \code{coding_status} columns.
#' @return A list with two elements:
#'   \describe{
#'     \item{pair_summary}{Tibble with one row per pair: \code{gene_id},
#'       \code{reference_isoform_id}, \code{comparator_isoform_id},
#'       \code{ref_start_codon}, \code{comp_start_codon},
#'       \code{same_start_codon}, \code{n_shared_cds_regions},
#'       \code{shared_cds_bp}, \code{frame_category},
#'       \code{pct_shared_cds_in_frame}, \code{pct_shared_cds_frameshifted}.}
#'     \item{region_detail}{Tibble with one row per shared CDS region per pair:
#'       \code{gene_id}, \code{reference_isoform_id},
#'       \code{comparator_isoform_id}, \code{region_start}, \code{region_end},
#'       \code{region_bp}, \code{ref_phase_at_start},
#'       \code{comp_phase_at_start}, \code{phase_match}.}
#'   }
#' @examples
#' gtf_path <- system.file("extdata", "test_cds_cases.gtf", package = "Isopair")
#' cds_exons <- extractCdsExons(gtf_path, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, verbose = FALSE)
#' pairs <- data.frame(
#'   gene_id = "GENE_CDS_A1",
#'   reference_isoform_id = "TX_A1_1",
#'   comparator_isoform_id = "TX_A1_1",
#'   stringsAsFactors = FALSE
#' )
#' result <- compareIsoformFrames(pairs, cds_exons, cds)
#' result$pair_summary$frame_category
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
compareIsoformFrames <- function(profiles, cds_exons, cds_metadata) {

  # Coding status lookup
  coding_ids <- cds_metadata$isoform_id[cds_metadata$coding_status == "coding"]

  # Pre-split cds_exons by isoform for fast lookup
  exon_list <- split(
    cds_exons[, c("cds_exon_start", "cds_exon_end", "strand")],
    cds_exons$isoform_id
  )

  pair_summaries <- list()
  region_details <- list()

  for (i in seq_len(nrow(profiles))) {
    gene <- profiles$gene_id[i]
    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]

    # Check coding status
    ref_coding <- ref_id %in% coding_ids
    comp_coding <- comp_id %in% coding_ids

    if (!ref_coding || !comp_coding) {
      pair_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        ref_start_codon = NA_integer_, comp_start_codon = NA_integer_,
        same_start_codon = NA,
        n_shared_cds_regions = 0L, shared_cds_bp = 0L,
        frame_category = "non_coding",
        pct_shared_cds_in_frame = NA_real_,
        pct_shared_cds_frameshifted = NA_real_
      )
      next
    }

    # Get per-exon CDS intervals
    ref_exons <- exon_list[[ref_id]]
    comp_exons <- exon_list[[comp_id]]

    if (is.null(ref_exons) || nrow(ref_exons) == 0L ||
        is.null(comp_exons) || nrow(comp_exons) == 0L) {
      pair_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        ref_start_codon = NA_integer_, comp_start_codon = NA_integer_,
        same_start_codon = NA,
        n_shared_cds_regions = 0L, shared_cds_bp = 0L,
        frame_category = "non_coding",
        pct_shared_cds_in_frame = NA_real_,
        pct_shared_cds_frameshifted = NA_real_
      )
      next
    }

    strand <- ref_exons$strand[1L]

    # Start codon positions
    ref_start <- .getStartCodonPos(ref_exons, strand)
    comp_start <- .getStartCodonPos(comp_exons, strand)
    same_start <- (ref_start == comp_start)

    # Intersect CDS intervals
    ref_intervals <- data.frame(
      start = ref_exons$cds_exon_start,
      end = ref_exons$cds_exon_end
    )
    comp_intervals <- data.frame(
      start = comp_exons$cds_exon_start,
      end = comp_exons$cds_exon_end
    )
    shared <- .intersectCdsIntervals(ref_intervals, comp_intervals)

    if (nrow(shared) == 0L) {
      pair_summaries[[i]] <- tibble::tibble(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        ref_start_codon = ref_start, comp_start_codon = comp_start,
        same_start_codon = same_start,
        n_shared_cds_regions = 0L, shared_cds_bp = 0L,
        frame_category = "no_shared_cds",
        pct_shared_cds_in_frame = NA_real_,
        pct_shared_cds_frameshifted = NA_real_
      )
      next
    }

    # Compute phase at start of each shared region for both isoforms
    n_regions <- nrow(shared)
    ref_phases <- integer(n_regions)
    comp_phases <- integer(n_regions)
    region_bp <- as.integer(shared$end - shared$start + 1L)

    for (r in seq_len(n_regions)) {
      # Phase at 5' end of shared region
      if (strand == "-") {
        phase_pos <- shared$end[r]
      } else {
        phase_pos <- shared$start[r]
      }
      ref_phases[r] <- .computePhaseAtPosition(ref_exons, strand, phase_pos)
      comp_phases[r] <- .computePhaseAtPosition(comp_exons, strand, phase_pos)
    }

    phase_match <- (ref_phases == comp_phases)
    bp_in_frame <- sum(region_bp[phase_match])
    bp_frameshifted <- sum(region_bp[!phase_match])
    total_shared_bp <- sum(region_bp)

    # Classify
    all_match <- all(phase_match)
    if (same_start) {
      frame_cat <- if (all_match) "same_start_in_frame" else "same_start_frameshift"
    } else {
      if (all_match) {
        frame_cat <- "diff_start_same_frame"
      } else {
        # Check if the phase offset is constant (start-codon-only) or changes
        # (splicing-induced frameshift on top of start offset)
        offsets <- (ref_phases - comp_phases) %% 3L
        if (length(unique(offsets)) == 1L) {
          frame_cat <- "diff_start_diff_frame"
        } else {
          frame_cat <- "diff_start_diff_frame_with_frameshift"
        }
      }
    }

    pair_summaries[[i]] <- tibble::tibble(
      gene_id = gene, reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      ref_start_codon = ref_start, comp_start_codon = comp_start,
      same_start_codon = same_start,
      n_shared_cds_regions = n_regions,
      shared_cds_bp = total_shared_bp,
      frame_category = frame_cat,
      pct_shared_cds_in_frame = round(100 * bp_in_frame / total_shared_bp, 1),
      pct_shared_cds_frameshifted = round(100 * bp_frameshifted / total_shared_bp, 1)
    )

    region_details[[i]] <- tibble::tibble(
      gene_id = gene, reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      region_start = shared$start,
      region_end = shared$end,
      region_bp = region_bp,
      ref_phase_at_start = ref_phases,
      comp_phase_at_start = comp_phases,
      phase_match = phase_match
    )
  }

  list(
    pair_summary = dplyr::bind_rows(pair_summaries),
    region_detail = dplyr::bind_rows(region_details)
  )
}
