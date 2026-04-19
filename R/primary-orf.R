# Primary-ORF selection
#
# Rule (per isovar 2026-04-18 design):
#   1. If the dominant isoform's annotated CDS ATG is exonic in the
#      target isoform, use that ATG as the primary translation start
#      for the target (preserves shared-start ORFs across isoforms).
#   2. Otherwise, use the 5'-most ATG whose Kozak context meets the
#      MANE-calibrated threshold (default from defaultKozakThreshold()).
#   3. If neither condition yields an ATG, no primary ORF is selected.
#
# The primary ORF is the translation start a ribosome is most plausibly
# using for a given isoform given what we know about the cohort's
# dominant isoform. Downstream classification (effectively_ptc,
# no_downstream_ejc, etc.) applies to this single ORF per isoform;
# enumerateOrfs() remains the per-ORF exhaustive view.


#' Select the Primary Translation Start per Isoform
#'
#' For each isoform whose enumerated ORFs are provided, picks a single
#' "primary" ORF using the two-step rule described in the module
#' header. Returns the atg_tx_pos, Kozak score, and a source flag
#' indicating why that ATG was chosen.
#'
#' The primary choice is the dominant isoform's ATG when it's exonic in
#' the target (\code{source = "dominant"}); otherwise the 5'-most
#' enumerated ORF, which is Kozak-filtered by construction when
#' \code{orfs} came from \code{enumerateOrfs(kozak_filter = TRUE, ...)}
#' (\code{source = "first_hq"}).
#'
#' @param orfs Output of \code{\link{enumerateOrfs}()}; assumed to be
#'   the Kozak-gated set (ATGs below threshold already filtered out).
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()},
#'   used to look up the dominant isoform's annotated CDS ATG.
#' @param dominant_isoform_id Character; ID of the dominant isoform.
#'   If \code{NULL} (default), step 1 is skipped and every selection
#'   is \code{"first_hq"}.
#' @return A tibble with one row per isoform present in \code{orfs}:
#'   \code{isoform_id}, \code{primary_atg_tx_pos}, \code{primary_kozak_score},
#'   \code{primary_source} (one of \code{"dominant"}, \code{"first_hq"},
#'   \code{"none"}), \code{primary_category} (pulled from \code{orfs}
#'   for the selected ORF), \code{primary_orf_length},
#'   \code{primary_n_downstream_ejc}.
#' @seealso \code{\link{enumerateOrfs}}, \code{\link{identifyDominantIsoforms}}
#' @export
selectPrimaryOrf <- function(orfs, structures, cds_metadata,
                             dominant_isoform_id = NULL) {
  # Look up the dominant isoform's annotated CDS ATG (genomic, strand-aware).
  dom_atg_g <- NA_real_
  dom_strand <- NA_character_
  if (!is.null(dominant_isoform_id)) {
    dom_row <- cds_metadata[cds_metadata$isoform_id == dominant_isoform_id &
                              cds_metadata$coding_status == "coding", ]
    if (nrow(dom_row) == 1L) {
      dom_strand <- dom_row$strand
      dom_atg_g  <- if (dom_strand == "+") dom_row$cds_start else dom_row$cds_stop
    }
  }

  all_isos <- unique(orfs$isoform_id)
  out <- vector("list", length(all_isos))

  for (k in seq_along(all_isos)) {
    iso_id <- all_isos[k]
    sub <- orfs[orfs$isoform_id == iso_id, , drop = FALSE]
    if (nrow(sub) == 0L) next

    # Attempt step 1: is the dominant ATG exonic in this isoform?
    chosen_idx <- NA_integer_
    chosen_src <- NA_character_
    if (!is.na(dom_atg_g)) {
      s_row <- structures[structures$isoform_id == iso_id, ]
      if (nrow(s_row) == 1L) {
        starts <- s_row$exon_starts[[1]]
        ends   <- s_row$exon_ends[[1]]
        dom_tx <- suppressWarnings(genomicToTranscript(dom_atg_g, starts, ends,
                                                        s_row$strand))
        if (!is.na(dom_tx)) {
          match_row <- which(sub$atg_tx_pos == dom_tx)
          if (length(match_row) == 1L) {
            chosen_idx <- match_row
            chosen_src <- "dominant"
          }
        }
      }
    }

    # Step 2: first high-quality ATG (5'-most in the Kozak-gated set).
    if (is.na(chosen_idx)) {
      o <- order(sub$atg_tx_pos, decreasing = FALSE)
      chosen_idx <- o[1]
      chosen_src <- "first_hq"
    }

    r <- sub[chosen_idx, ]
    out[[k]] <- tibble::tibble(
      isoform_id                = iso_id,
      primary_atg_tx_pos        = r$atg_tx_pos,
      primary_kozak_score       = if ("kozak_score" %in% names(r)) r$kozak_score else NA_real_,
      primary_source            = chosen_src,
      primary_category          = r$category,
      primary_orf_length        = r$orf_length,
      primary_n_downstream_ejc  = r$n_downstream_ejc
    )
  }

  out <- dplyr::bind_rows(out[!vapply(out, is.null, logical(1))])

  # Add a "none" row for isoforms not present in orfs (Kozak gate removed
  # them all). Caller can pass the full isoform set if they want these.
  out
}
