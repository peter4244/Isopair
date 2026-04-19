# Per-ORF enumeration + classification
#
# Isopair's per-transcript functions (traceReferenceAtg, computePtcStatus,
# attributePtcEvents) evaluate a single ORF per isoform-or-pair. In biology,
# a transcript can host multiple ORFs (main CDS, uORFs, alt-start internal
# ORFs, overlapping frames) with independent translation and NMD fates.
# enumerateOrfs() provides the per-ORF view: for each transcript in the
# input, every viable ATG is enumerated and its ORF classified.
#
# This is the data shape downstream per-ORF analyses should start from. The
# existing per-transcript functions are convenience aggregations over this
# per-ORF space.


#' Enumerate and Classify Every Viable ORF in a Set of Transcripts
#'
#' For each isoform in \code{structures}, scans the transcript sequence for
#' every ATG codon, traces each one to its first in-frame stop, and — for
#' ORFs meeting the length threshold — classifies the ORF as NMD-triggering
#' (\code{effectively_ptc}), terminating cleanly
#' (\code{no_downstream_ejc}), or running off the transcript without a stop
#' (\code{no_stop_in_frame}).
#'
#' Returns one row per (isoform, viable ORF) — NOT one row per transcript.
#' Transcript-level questions ("is this isoform an NMD substrate?") are
#' group-by roll-ups over this per-ORF output.
#'
#' Optionally scores the Kozak context of each ATG via
#' \code{\link{scoreKozakPWM}()} — when \code{kozak_weights} is supplied, a
#' \code{kozak_score} column is added; when \code{kozak_score_min} is also
#' set, ATGs below that threshold are excluded.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()},
#'   used only to flag whether a given ATG coincides with the isoform's own
#'   annotated CDS start (adds an \code{is_annotated_cds} column).
#' @param sequences A named character vector of transcript sequences (names =
#'   isoform IDs, values = uppercase DNA in transcript orientation).
#' @param min_orf_nt Integer; minimum ORF length (nt, ATG-inclusive to stop-
#'   exclusive) to emit a row. Default 30 (10 aa).
#' @param ejc_threshold Integer; minimum nt distance between stop codon 3'
#'   end and next downstream junction for NMD triggering. Default 50.
#' @param include_no_stop Logical; when \code{TRUE}, emit rows for ATGs whose
#'   reading frame runs off the transcript without hitting a stop, categorized
#'   as \code{"no_stop_in_frame"}. Default \code{TRUE}.
#' @param kozak_weights Optional PWM passed to \code{\link{scoreKozakPWM}()}.
#'   When supplied, adds a \code{kozak_score} column.
#' @param kozak_score_min Optional numeric; when supplied with
#'   \code{kozak_weights}, ATGs with scores below this threshold are filtered
#'   out.
#' @return A tibble with columns:
#'   \code{isoform_id}, \code{atg_tx_pos} (1-based transcript position),
#'   \code{stop_tx_pos} (NA when the ORF runs off the transcript),
#'   \code{orf_length} (nt; NA when no_stop),
#'   \code{n_downstream_ejc}, \code{is_annotated_cds} (TRUE when this ATG
#'   matches the isoform's own CDS start from \code{cds_metadata}),
#'   \code{kozak_score} (when \code{kozak_weights} is supplied),
#'   \code{category} — one of \code{"effectively_ptc"},
#'   \code{"no_downstream_ejc"}, \code{"no_stop_in_frame"}.
#' @seealso \code{\link{traceReferenceAtg}}, \code{\link{computePtcStatus}}
#' @examples
#' \dontrun{
#' orfs <- enumerateOrfs(structures, cds, sequences)
#' # Transcript is NMD substrate if ANY ORF is effectively_ptc
#' is_nmd <- orfs |>
#'   dplyr::group_by(isoform_id) |>
#'   dplyr::summarise(any_ptc = any(category == "effectively_ptc"))
#' }
#' @export
enumerateOrfs <- function(structures, cds_metadata, sequences,
                          min_orf_nt      = 30L,
                          ejc_threshold   = 50L,
                          include_no_stop = TRUE,
                          kozak_weights   = NULL,
                          kozak_score_min = NULL) {
  stop_codons <- c("TAA", "TAG", "TGA")

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

  out <- dplyr::bind_rows(out_parts)
  if (nrow(out) == 0L) return(out)

  # Optional Kozak scoring / filtering.
  if (!is.null(kozak_weights)) {
    seqs_vec <- sequences[out$isoform_id]
    out$kozak_score <- scoreKozakPWM(seqs_vec, out$atg_tx_pos,
                                      weights = kozak_weights)
    if (!is.null(kozak_score_min)) {
      out <- out[!is.na(out$kozak_score) & out$kozak_score >= kozak_score_min, ,
                 drop = FALSE]
    }
  }

  out
}
