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
#' For each splicing event in a set of profiles, computes CDS overlap and
#' determines whether the event disrupts the reading frame (CDS overlap
#' not divisible by 3). Only frameshift-capable event types are evaluated:
#' SE, Missing_Internal, A5SS, A3SS, Partial_IR_5, Partial_IR_3.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} containing
#'   a \code{detailed_events} list column.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with columns: \code{isoform_id}, \code{coding_status},
#'   \code{cds_start}, \code{cds_stop}.
#' @return A list with two elements:
#'   \describe{
#'     \item{events}{Tibble of all events with added columns:
#'       \code{cds_overlap_bp}, \code{affects_cds}, \code{frame_disrupting}}
#'     \item{profile_summary}{Tibble with per-pair summary: \code{gene_id},
#'       \code{reference_isoform_id}, \code{comparator_isoform_id},
#'       \code{n_cds_events}, \code{n_frame_disrupting},
#'       \code{has_frame_disruption}}
#'   }
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

  frameshift_types <- c("SE", "Missing_Internal", "A5SS", "A3SS",
                         "Partial_IR_5", "Partial_IR_3")

  # Build CDS lookup
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

    # Compute CDS overlap for each event
    evt_results <- events
    n_evt <- nrow(evt_results)

    cds_overlap_bp <- integer(n_evt)
    affects_cds <- logical(n_evt)
    frame_disrupting <- logical(n_evt)

    if (has_ref_cds) {
      ref_cds_start <- ref_cds$cds_start[1]
      ref_cds_stop <- ref_cds$cds_stop[1]

      for (j in seq_len(n_evt)) {
        # Get genomic coordinates of event
        event_start <- pmin(evt_results$five_prime[j],
                            evt_results$three_prime[j])
        event_end <- pmax(evt_results$five_prime[j],
                          evt_results$three_prime[j])

        # Compute overlap with CDS
        overlap_start <- pmax(event_start, ref_cds_start)
        overlap_end <- pmin(event_end, ref_cds_stop)
        overlap_bp <- pmax(0L, as.integer(overlap_end - overlap_start + 1L))
        cds_overlap_bp[j] <- overlap_bp
        affects_cds[j] <- overlap_bp > 0L

        # Frame disruption: only for frameshift-capable events
        if (affects_cds[j] &&
            evt_results$event_type[j] %in% frameshift_types) {
          frame_disrupting[j] <- (overlap_bp %% 3L) != 0L
        }
      }
    }

    evt_results$cds_overlap_bp <- cds_overlap_bp
    evt_results$affects_cds <- affects_cds
    evt_results$frame_disrupting <- frame_disrupting

    all_events[[i]] <- evt_results

    # Profile-level summary
    fs_events <- evt_results[evt_results$event_type %in% frameshift_types, ]
    n_cds <- sum(fs_events$affects_cds)
    n_fd <- sum(fs_events$frame_disrupting)

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
