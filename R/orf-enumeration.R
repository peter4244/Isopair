# Per-ORF enumeration + classification
#
# Isopair's per-transcript functions (traceReferenceAtg, computePtcStatus,
# attributePtcEvents) evaluate a single ORF per isoform-or-pair. In biology,
# a transcript can host multiple ORFs (main CDS, uORFs, alt-start internal
# ORFs, overlapping frames) with independent translation and NMD fates.
# enumerateOrfs() provides the per-ORF view: for each transcript in the
# input, every *plausibly translatable* ATG (Kozak-passing by default) is
# enumerated and its ORF classified.
#
# Kozak-first design: a transcript contains many ATGs but only those with
# reasonable initiation context are plausible translation starts. The
# function scores every ATG with the PWM first, filters to the threshold,
# and only then walks the reading frame to classify NMD fate. This is the
# biologically principled unit for downstream per-ORF analyses. Pass
# kozak_filter = FALSE to recover the unfiltered scan (all ATGs above
# min_orf_nt), e.g. for benchmarking or exhaustive enumeration.


#' Derive an empirical Kozak-score threshold from annotated CDS starts
#'
#' Scores the Kozak context of every annotated CDS start in the input
#' (one per coding isoform in \code{cds_metadata}, where both the
#' transcript sequence and the mapped start position are available) and
#' returns the requested lower-tail quantile of those scores.
#'
#' The returned value is the data-driven "plausibly translatable"
#' threshold — by default the 5th percentile of real annotated starts,
#' so any ATG scoring at or above this threshold is "at least as
#' initiation-competent as 95% of annotated CDS starts" in the training
#' set. Pass this directly as \code{kozak_threshold} in
#' \code{\link{enumerateOrfs}()}.
#'
#' This is more principled than picking a fixed number by hand because
#' the PWM's expected score range depends on the weights in use and the
#' composition of the input transcriptome.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @param sequences Named character vector of transcript sequences.
#' @param quantile Numeric in \code{(0, 1)}; default \code{0.05} (5th
#'   percentile). Use \code{0.10} for a stricter threshold,
#'   \code{0.01} for more inclusive.
#' @param weights Optional PWM; defaults to \code{.defaultKozakPWM()}.
#' @return A named list with \code{threshold} (the quantile cutoff),
#'   \code{n_used} (annotated starts successfully scored),
#'   \code{scores} (numeric vector of all annotated-start scores).
#' @seealso \code{\link{enumerateOrfs}}, \code{\link{scoreKozakPWM}}
#' @examples
#' \dontrun{
#' thr <- empiricalKozakThreshold(structures, cds, sequences,
#'                                 quantile = 0.05)
#' orfs <- enumerateOrfs(structures, cds, sequences,
#'                       kozak_threshold = thr$threshold)
#' }
#' @export
empiricalKozakThreshold <- function(structures, cds_metadata, sequences,
                                    quantile = 0.05,
                                    weights  = NULL) {
  if (is.null(weights)) weights <- .defaultKozakPWM()
  stopifnot(is.numeric(quantile), length(quantile) == 1L,
            quantile > 0, quantile < 1)

  coding <- cds_metadata[cds_metadata$coding_status == "coding", , drop = FALSE]
  annot_atg_g <- stats::setNames(
    ifelse(coding$strand == "+", coding$cds_start, coding$cds_stop),
    coding$isoform_id
  )

  scores <- numeric(0)
  for (i in seq_len(nrow(structures))) {
    iso_id <- structures$isoform_id[i]
    if (!iso_id %in% names(sequences))     next
    if (!iso_id %in% names(annot_atg_g))   next
    atg_g <- annot_atg_g[iso_id]
    if (is.na(atg_g))                      next

    atg_tx <- suppressWarnings(genomicToTranscript(
      atg_g,
      structures$exon_starts[[i]], structures$exon_ends[[i]],
      structures$strand[i]
    ))
    if (is.na(atg_tx)) next

    sc_df <- scoreKozakPWM(sequences[iso_id], atg_positions = atg_tx,
                           weights = weights)
    sc <- if (is.data.frame(sc_df)) sc_df$score[1] else sc_df[1]
    if (!is.na(sc)) scores <- c(scores, sc)
  }

  thr <- if (length(scores) > 0L)
    as.numeric(stats::quantile(scores, probs = quantile, na.rm = TRUE)) else NA_real_

  list(threshold = thr, n_used = length(scores), scores = scores)
}


#' Enumerate and Classify Plausibly-Translated ORFs in a Set of Transcripts
#'
#' For each isoform in \code{structures}, scores the Kozak context of every
#' ATG in the transcript, keeps only ATGs above the threshold, and for those
#' traces the reading frame to the first in-frame stop, classifying the ORF
#' as NMD-triggering (\code{effectively_ptc}), terminating cleanly
#' (\code{no_downstream_ejc}), or running off the transcript without a stop
#' (\code{no_stop_in_frame}).
#'
#' Kozak filtering is the default because only ATGs with reasonable
#' initiation context are plausibly translated; enumerating every ATG
#' produces mostly noise. The default PWM comes from Isopair's internal
#' \code{.defaultKozakPWM()} (log-odds at positions -6..-1, +4, +5), and
#' the default score threshold \code{kozak_threshold = 0} admits ATGs
#' above the null (random-sequence) baseline. A strong-canonical Kozak
#' typically scores ~+2 or higher; annotated CDS starts are usually >= +1.
#' Tighten the threshold to make the per-ORF output more biology-weighted;
#' loosen (or set \code{kozak_filter = FALSE}) to recover every ATG-indexed
#' ORF above \code{min_orf_nt}.
#'
#' Returns one row per (isoform, plausibly-translated ORF) — NOT one row
#' per transcript. Transcript-level questions ("is this isoform an NMD
#' substrate?") are group-by roll-ups over this per-ORF output.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()},
#'   used to flag whether each ATG matches the isoform's annotated CDS
#'   start (\code{is_annotated_cds}).
#' @param sequences A named character vector of transcript sequences
#'   (names = isoform IDs; uppercase DNA in transcript orientation).
#' @param min_orf_nt Integer; minimum ORF length (nt, ATG-inclusive to
#'   stop-exclusive) to emit a row. Default 30 (10 aa).
#' @param ejc_threshold Integer; minimum nt distance between stop codon 3'
#'   end and next downstream junction for NMD triggering. Default 50.
#' @param include_no_stop Logical; emit rows for ATGs whose reading frame
#'   runs off the transcript (category \code{no_stop_in_frame}). Default
#'   \code{TRUE}.
#' @param kozak_filter Logical; when \code{TRUE} (default), score every
#'   ATG's Kozak context and retain only those at or above
#'   \code{kozak_threshold} before ORF tracing. When \code{FALSE}, every
#'   ATG above \code{min_orf_nt} is emitted.
#' @param kozak_threshold Numeric; log-odds Kozak score threshold. Default
#'   \code{0} (above random). Real CDS starts generally score >= 1.
#' @param kozak_weights Optional PWM passed to
#'   \code{\link{scoreKozakPWM}()}. When \code{NULL} (default), Isopair's
#'   internal default PWM is used.
#' @return A tibble with columns: \code{isoform_id}, \code{atg_tx_pos},
#'   \code{kozak_score}, \code{stop_tx_pos} (NA when no_stop),
#'   \code{orf_length} (NA when no_stop), \code{n_downstream_ejc},
#'   \code{is_annotated_cds}, \code{category} (one of
#'   \code{effectively_ptc}, \code{no_downstream_ejc},
#'   \code{no_stop_in_frame}).
#' @seealso \code{\link{traceReferenceAtg}}, \code{\link{computePtcStatus}},
#'   \code{\link{scoreKozakPWM}}
#' @examples
#' \dontrun{
#' orfs <- enumerateOrfs(structures, cds, sequences,
#'                       kozak_threshold = 0.5)
#' # Transcript is NMD substrate if any plausibly-translated ORF is PTC
#' dplyr::group_by(orfs, isoform_id) |>
#'   dplyr::summarise(any_ptc = any(category == "effectively_ptc"))
#' }
#' @export
enumerateOrfs <- function(structures, cds_metadata, sequences,
                          min_orf_nt      = 30L,
                          ejc_threshold   = 50L,
                          include_no_stop = TRUE,
                          kozak_filter    = TRUE,
                          kozak_threshold = 0,
                          kozak_weights   = NULL) {
  stop_codons <- c("TAA", "TAG", "TGA")

  if (kozak_filter && is.null(kozak_weights))
    kozak_weights <- .defaultKozakPWM()

  # Per-isoform annotated-CDS-start lookup (strand-aware ATG).
  coding <- cds_metadata[cds_metadata$coding_status == "coding", ]
  annot_atg_g <- stats::setNames(
    ifelse(coding$strand == "+", coding$cds_start, coding$cds_stop),
    coding$isoform_id
  )

  out_parts <- vector("list", nrow(structures))
  for (i in seq_len(nrow(structures))) {
    iso_id <- structures$isoform_id[i]
    if (!iso_id %in% names(sequences)) next
    starts <- structures$exon_starts[[i]]
    ends   <- structures$exon_ends[[i]]
    strand <- structures$strand[i]
    seq_i  <- sequences[iso_id]

    # Find every ATG (transcript-space, 1-based positions).
    m <- gregexpr("ATG", seq_i, fixed = TRUE)[[1]]
    if (length(m) == 1L && m[1] == -1L) next
    atg_positions <- as.integer(m)

    # Score Kozak for every ATG, then gate by threshold.
    kozak_scores <- rep(NA_real_, length(atg_positions))
    if (!is.null(kozak_weights)) {
      sc_df <- scoreKozakPWM(rep(seq_i, length(atg_positions)),
                             atg_positions = atg_positions,
                             weights = kozak_weights)
      kozak_scores <- if (is.data.frame(sc_df)) sc_df$score else sc_df
    }
    if (kozak_filter) {
      keep <- !is.na(kozak_scores) & kozak_scores >= kozak_threshold
      atg_positions <- atg_positions[keep]
      kozak_scores  <- kozak_scores[keep]
      if (length(atg_positions) == 0L) next
    }

    # Translate annotated-CDS ATG (genomic) to transcript pos for this iso.
    annot_tx_pos <- NA_integer_
    if (iso_id %in% names(annot_atg_g) && !is.na(annot_atg_g[iso_id])) {
      annot_tx_pos <- suppressWarnings(
        genomicToTranscript(annot_atg_g[iso_id], starts, ends, strand)
      )
    }

    junctions <- .getJunctionPositions(starts, ends, strand)

    rows <- vector("list", length(atg_positions))
    for (k in seq_along(atg_positions)) {
      atg_tx <- atg_positions[k]
      stop_p <- .findFirstStop(seq_i, atg_tx, stop_codons)

      if (is.na(stop_p)) {
        if (!include_no_stop) next
        cat_val <- "no_stop_in_frame"
        orf_len <- NA_integer_
        n_dejc  <- NA_integer_
      } else {
        orf_len <- stop_p - atg_tx
        if (orf_len < min_orf_nt) next
        stop_end <- stop_p + 2L
        n_dejc   <- sum(junctions > (stop_end + ejc_threshold - 1L))
        cat_val  <- if (n_dejc > 0L) "effectively_ptc" else "no_downstream_ejc"
      }

      rows[[k]] <- tibble::tibble(
        isoform_id       = iso_id,
        atg_tx_pos       = atg_tx,
        kozak_score      = kozak_scores[k],
        stop_tx_pos      = if (is.na(stop_p)) NA_integer_ else as.integer(stop_p),
        orf_length       = if (is.na(orf_len)) NA_integer_ else as.integer(orf_len),
        n_downstream_ejc = if (is.na(n_dejc))  NA_integer_ else as.integer(n_dejc),
        is_annotated_cds = !is.na(annot_tx_pos) && atg_tx == annot_tx_pos,
        category         = cat_val
      )
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) out_parts[[i]] <- dplyr::bind_rows(rows)
  }

  dplyr::bind_rows(out_parts)
}
