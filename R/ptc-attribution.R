# PTC-Causing Event Attribution
#
# Functions for attributing premature termination codons to specific splice
# events, and for tracing reference ATG reading frames through comparator
# isoforms to identify hidden PTCs.


#' Attribute PTC-Causing Splice Events
#'
#' For each PTC-positive isoform pair, identifies the specific splice event
#' responsible for the premature stop codon using a three-mechanism framework:
#' \enumerate{
#'   \item \strong{Frameshift}: a splice event shifts the reading frame,
#'     creating a premature stop codon in the new frame (identified by
#'     \code{\link{analyzeFrameWalk}()}).
#'   \item \strong{In-frame stop}: a splice event introduces exonic sequence
#'     containing an in-frame stop codon (event coordinates contain the PTC
#'     position).
#'   \item \strong{3'UTR splice}: the comparator and reference share the same
#'     stop codon, but a downstream splice event repositions an EJC, making the
#'     shared stop premature (attributed by \code{\link{attribute3UtrSplice}()}).
#' }
#'
#' Split-codon handling: a stop codon can span a splice junction (2 nt from
#' one exon, 1 from the next). Coordinate containment checks use a +/-2 bp
#' buffer, and ATG-to-stop region filtering uses overlap semantics with a 3 bp
#' extension.
#'
#' @param pairs A data frame with columns \code{comparator_isoform_id} and
#'   \code{reference_isoform_id}. Each row is one PTC-positive pair.
#' @param fw_events A data frame from \code{\link{analyzeFrameWalk}()$events}
#'   with columns: \code{reference_isoform_id}, \code{comparator_isoform_id},
#'   \code{event_type}, \code{genomic_start}, \code{genomic_end},
#'   \code{is_frameshift}.
#' @param profiles A data frame from \code{\link{buildProfiles}()} with a
#'   \code{detailed_events} list-column and \code{comparator_isoform_id}.
#' @param ptc_genomic_pos Named numeric vector: comparator_isoform_id to
#'   strand-aware genomic position of the premature stop codon. Use
#'   \code{\link{getStrandAwareStop}()} to compute this.
#' @param atg_genomic_pos Optional named numeric vector: comparator_isoform_id
#'   to genomic position of the CDS start ATG. If provided, frameshift events
#'   are filtered to those overlapping the ATG-to-stop region. Use \code{NULL}
#'   (default) to consider all frameshift events for a pair.
#' @param strand_vec Named character vector: comparator_isoform_id to strand
#'   ("+" or "-").
#' @param is_frameshift_vec Optional named logical vector:
#'   comparator_isoform_id to whether the pair is in a frameshift frame
#'   category. If \code{NULL} (default), determined from \code{fw_events}
#'   (any frameshift event for the pair implies TRUE).
#' @return A tibble with columns: \code{comparator_isoform_id},
#'   \code{mechanism} ("Frameshift", "In-frame stop", or "unresolved"),
#'   \code{ptc_causing_event} (event type name or "Unresolved"),
#'   \code{attribution} ("direct" or "unresolved").
#' @seealso \code{\link{analyzeFrameWalk}}, \code{\link{computePtcStatus}},
#'   \code{\link{attribute3UtrSplice}}
#' @examples
#' # See vignette("Isopair") for a complete workflow
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
attributePtcEvents <- function(pairs,
                                fw_events,
                                profiles,
                                ptc_genomic_pos,
                                atg_genomic_pos = NULL,
                                strand_vec,
                                is_frameshift_vec = NULL) {

  # Key on ref::comp to handle cases where the same comparator appears
  # in multiple pairs with different references (and different events)
  prof_keys <- paste0(profiles$reference_isoform_id, "::",
                       profiles$comparator_isoform_id)
  prof_idx_lookup <- stats::setNames(seq_len(nrow(profiles)), prof_keys)

  result_rows <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    comp_id <- pairs$comparator_isoform_id[i]
    ref_id  <- pairs$reference_isoform_id[i]
    stop_g  <- ptc_genomic_pos[comp_id]
    strand  <- strand_vec[comp_id]

    if (is.na(stop_g) || is.na(strand)) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "unresolved",
        ptc_causing_event = "no_stop_pos", attribution = "unresolved")
      next
    }

    # Determine if this pair is in a frameshift category
    if (!is.null(is_frameshift_vec)) {
      is_fs <- isTRUE(is_frameshift_vec[comp_id])
    } else {
      pair_fw <- fw_events[fw_events$reference_isoform_id == ref_id &
                             fw_events$comparator_isoform_id == comp_id, ]
      is_fs <- any(pair_fw$is_frameshift, na.rm = TRUE)
    }

    # --- Mechanism 1: Frameshift ---
    if (is_fs) {
      pair_fw <- fw_events[fw_events$reference_isoform_id == ref_id &
                             fw_events$comparator_isoform_id == comp_id, ]

      if (nrow(pair_fw) > 0L) {
        # If ATG position provided, filter to events overlapping ATG-to-stop
        # region (+3 bp buffer for split-codon splice junctions)
        if (!is.null(atg_genomic_pos)) {
          atg_g <- atg_genomic_pos[comp_id]
          if (!is.na(atg_g)) {
            region_lo <- min(atg_g, stop_g) - 3L
            region_hi <- max(atg_g, stop_g) + 3L
            pair_fw <- pair_fw[pair_fw$genomic_start <= region_hi &
                                 pair_fw$genomic_end >= region_lo, ]
          }
        }

        fs_events <- pair_fw[pair_fw$is_frameshift == TRUE, ]
        if (nrow(fs_events) > 0L) {
          # First frameshift event in transcription direction
          if (!is.na(strand) && strand == "+") {
            first_fs <- fs_events[which.min(fs_events$genomic_start), ]
          } else {
            first_fs <- fs_events[which.max(fs_events$genomic_end), ]
          }
          result_rows[[i]] <- tibble::tibble(
            comparator_isoform_id = comp_id, mechanism = "Frameshift",
            ptc_causing_event = first_fs$event_type[1],
            attribution = "direct")
          next
        }
      }

      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "Frameshift",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    # --- Mechanism 2: In-frame stop ---
    pair_key <- paste0(ref_id, "::", comp_id)
    prof_idx <- prof_idx_lookup[pair_key]
    if (is.na(prof_idx)) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    de <- profiles$detailed_events[[prof_idx]]
    if (is.null(de) || nrow(de) == 0L) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    # ±2 bp buffer for split-codon splice junctions
    ev_min <- pmin(de$five_prime, de$three_prime)
    ev_max <- pmax(de$five_prime, de$three_prime)
    containing <- de[stop_g >= (ev_min - 2L) & stop_g <= (ev_max + 2L), ]

    if (nrow(containing) > 0L) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = containing$event_type[1],
        attribution = "direct")
    } else {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "In-frame stop",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
    }
  }

  dplyr::bind_rows(result_rows)
}


#' Attribute 3'UTR Splice PTC-Causing Events
#'
#' For same-stop PTC-positive pairs (or pairs with downstream EJCs but no
#' truncated ORF), the PTC arises from a splice event downstream of the stop
#' codon that repositions an EJC. This function finds the nearest downstream
#' event.
#'
#' @param pairs A data frame with \code{comparator_isoform_id} and
#'   \code{reference_isoform_id}.
#' @param profiles A data frame from \code{\link{buildProfiles}()} with
#'   \code{detailed_events} list-column.
#' @param stop_genomic_pos Named numeric vector: comparator_isoform_id to stop
#'   codon genomic position.
#' @param strand_vec Named character vector: comparator_isoform_id to strand.
#' @return A tibble with columns: \code{comparator_isoform_id},
#'   \code{mechanism} (always "3'UTR splice"), \code{ptc_causing_event},
#'   \code{attribution}.
#' @seealso \code{\link{attributePtcEvents}}
#' @examples
#' # See vignette("Isopair") for a complete workflow
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
attribute3UtrSplice <- function(pairs, profiles, stop_genomic_pos, strand_vec) {

  # Key on ref::comp (same as attributePtcEvents) to handle duplicate comparators
  prof_keys <- paste0(profiles$reference_isoform_id, "::",
                       profiles$comparator_isoform_id)
  prof_idx_lookup <- stats::setNames(seq_len(nrow(profiles)), prof_keys)

  result_rows <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    comp_id <- pairs$comparator_isoform_id[i]
    ref_id  <- pairs$reference_isoform_id[i]
    stop_g  <- stop_genomic_pos[comp_id]
    strand  <- strand_vec[comp_id]

    if (is.na(stop_g) || is.na(strand)) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    pair_key <- paste0(ref_id, "::", comp_id)
    prof_idx <- prof_idx_lookup[pair_key]
    if (is.na(prof_idx)) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    de <- profiles$detailed_events[[prof_idx]]
    if (is.null(de) || nrow(de) == 0L) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    # Find events downstream of stop in transcription direction.
    # Use minimum absolute distance to either event boundary to handle
    # events that straddle the stop codon position.
    if (strand == "+") {
      utr3_ev <- de[de$five_prime >= stop_g | de$three_prime >= stop_g, ]
    } else {
      utr3_ev <- de[de$five_prime <= stop_g | de$three_prime <= stop_g, ]
    }

    if (nrow(utr3_ev) == 0L) {
      result_rows[[i]] <- tibble::tibble(
        comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
        ptc_causing_event = "Unresolved", attribution = "unresolved")
      next
    }

    dist_to_stop <- pmin(abs(utr3_ev$five_prime - stop_g),
                          abs(utr3_ev$three_prime - stop_g))
    nearest_idx <- which.min(dist_to_stop)

    result_rows[[i]] <- tibble::tibble(
      comparator_isoform_id = comp_id, mechanism = "3'UTR splice",
      ptc_causing_event = utr3_ev$event_type[nearest_idx],
      attribution = "direct")
  }

  dplyr::bind_rows(result_rows)
}


#' Trace Reference ATG Through Comparator Isoform
#'
#' For each isoform pair, checks whether the reference isoform's CDS start ATG
#' is exonic in the comparator, traces the reading frame to the first in-frame
#' stop codon, and counts downstream exon-exon junctions. This identifies
#' "hidden PTCs" — comparator isoforms that appear PTC-negative under their
#' own CDS prediction but are effectively PTC-positive when analyzed from the
#' reference reading frame.
#'
#' @param pairs A data frame with \code{reference_isoform_id} and
#'   \code{comparator_isoform_id}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @param sequences A named character vector of transcript sequences (names =
#'   isoform IDs, values = uppercase DNA sequences in transcript orientation).
#' @param ejc_threshold Integer; minimum distance (nt) between stop codon and
#'   downstream junction for NMD triggering (default 50).
#' @param resolve_alt_start Logical; if TRUE, when the reference ATG is not
#'   exonic in the comparator (the former \code{ref_atg_lost} terminal case),
#'   locate the first viable ATG in the comparator and continue the trace
#'   from that alternative start. Default FALSE preserves legacy behavior.
#' @param min_alt_orf_nt Integer; minimum alt-start ORF length (nt) to count
#'   as viable. Default 30 (10 aa).
#' @return A tibble with one row per pair and columns:
#'   \code{reference_isoform_id}, \code{comparator_isoform_id},
#'   \code{ref_atg_genomic} (strand-aware ATG position),
#'   \code{ref_atg_exonic_in_comp} (logical),
#'   \code{ref_orf_length} (nt in reference),
#'   \code{comp_orf_length} (nt in comparator from ref ATG),
#'   \code{comp_stop_tx_pos} (transcript position of stop in comparator),
#'   \code{n_downstream_ejc} (junctions downstream of stop),
#'   \code{alt_start_tx_pos} (transcript position of alt ATG when used; NA
#'   otherwise), \code{alt_start_orf_length} (nt; same condition),
#'   \code{category} (one of: "effectively_ptc", "truncated_no_ejc",
#'   "ref_atg_lost", "no_downstream_ejc", "no_ref_cds", "mapping_failed",
#'   and when \code{resolve_alt_start = TRUE}:
#'   "alt_start_effectively_ptc", "alt_start_no_downstream_ejc",
#'   "ref_atg_lost_no_viable_start").
#' @seealso \code{\link{computePtcStatus}}, \code{\link{attributePtcEvents}}
#' @examples
#' # See vignette("Isopair") for a complete workflow
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
traceReferenceAtg <- function(pairs, structures, cds_metadata, sequences,
                               ejc_threshold = 50L,
                               resolve_alt_start = FALSE,
                               min_alt_orf_nt = 30L) {

  # Build lookups
  structs_lookup <- stats::setNames(
    lapply(seq_len(nrow(structures)), function(i) {
      list(starts = structures$exon_starts[[i]],
           ends = structures$exon_ends[[i]],
           strand = structures$strand[i])
    }), structures$isoform_id)

  coding <- cds_metadata[cds_metadata$coding_status == "coding", ]
  cds_5prime <- stats::setNames(
    ifelse(coding$strand == "+", coding$cds_start, coding$cds_stop),
    coding$isoform_id)
  cds_strand <- stats::setNames(coding$strand, coding$isoform_id)

  stop_codons <- c("TAA", "TAG", "TGA")

  results <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    ref_id  <- pairs$reference_isoform_id[i]
    comp_id <- pairs$comparator_isoform_id[i]

    ref_atg_g <- cds_5prime[ref_id]
    strand    <- cds_strand[ref_id]

    # No reference CDS
    if (is.na(ref_atg_g) || is.na(strand)) {
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = NA_real_, ref_atg_exonic_in_comp = NA,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "no_ref_cds")
      next
    }

    # Check if ATG codon is exonic in comparator
    comp_s <- structs_lookup[[comp_id]]
    if (is.null(comp_s)) {
      # Missing comparator structures — can't resolve alt start.
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = FALSE,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "ref_atg_lost")
      next
    }

    atg_exonic <- .isCodonExonic(ref_atg_g, comp_s$starts, comp_s$ends,
                                  strand)
    if (!atg_exonic) {
      # Reference ATG not in any comparator exon. Classical "ref_atg_lost"
      # behavior unless resolve_alt_start = TRUE, in which case we scan
      # the comparator for the first viable ATG and trace from it.
      if (resolve_alt_start && comp_id %in% names(sequences)) {
        comp_seq_local <- sequences[comp_id]
        alt <- .findFirstViableATG(comp_seq_local, stop_codons,
                                   min_orf_nt = min_alt_orf_nt)
        if (!is.na(alt$atg_tx)) {
          junctions <- .getJunctionPositions(comp_s$starts, comp_s$ends, strand)
          stop_end  <- alt$stop_tx + 2L
          n_dejc    <- sum(junctions > (stop_end + ejc_threshold - 1L))
          cat_val <- if (n_dejc > 0L) "alt_start_effectively_ptc"
                     else "alt_start_no_downstream_ejc"
          results[[i]] <- tibble::tibble(
            reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
            ref_atg_genomic = as.numeric(ref_atg_g),
            ref_atg_exonic_in_comp = FALSE,
            ref_orf_length = NA_integer_,
            comp_orf_length = NA_integer_,
            comp_stop_tx_pos = NA_integer_,
            n_downstream_ejc = n_dejc,
            alt_start_tx_pos = alt$atg_tx,
            alt_start_orf_length = alt$orf_length,
            category = cat_val)
          next
        }
        # No viable alt start found
        results[[i]] <- tibble::tibble(
          reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
          ref_atg_genomic = as.numeric(ref_atg_g),
          ref_atg_exonic_in_comp = FALSE,
          ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
          comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
          alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
          category = "ref_atg_lost_no_viable_start")
        next
      }
      # Default (resolve_alt_start = FALSE): original terminal ref_atg_lost
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = FALSE,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "ref_atg_lost")
      next
    }

    # Map ATG to transcript positions
    ref_s <- structs_lookup[[ref_id]]
    if (is.null(ref_s) || !ref_id %in% names(sequences) ||
        !comp_id %in% names(sequences)) {
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = TRUE,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "mapping_failed")
      next
    }

    ref_tx  <- genomicToTranscript(ref_atg_g, ref_s$starts, ref_s$ends, strand)
    comp_tx <- genomicToTranscript(ref_atg_g, comp_s$starts, comp_s$ends, strand)

    if (is.na(ref_tx) || is.na(comp_tx)) {
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = TRUE,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "mapping_failed") ## dup above — second mapping_failed path
      next
    }

    ref_seq  <- sequences[ref_id]
    comp_seq <- sequences[comp_id]

    # Verify ATG
    if (substr(ref_seq, ref_tx, ref_tx + 2L) != "ATG" ||
        substr(comp_seq, comp_tx, comp_tx + 2L) != "ATG") {
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = TRUE,
        ref_orf_length = NA_integer_, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
        category = "mapping_failed")
      next
    }

    # Walk ORF in reference
    ref_stop <- .findFirstStop(ref_seq, ref_tx, stop_codons)
    ref_orf_len <- if (!is.na(ref_stop)) ref_stop - ref_tx else NA_integer_

    # Walk ORF in comparator
    comp_stop <- .findFirstStop(comp_seq, comp_tx, stop_codons)
    if (is.na(comp_stop)) {
      results[[i]] <- tibble::tibble(
        reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
        ref_atg_genomic = as.numeric(ref_atg_g),
        ref_atg_exonic_in_comp = TRUE,
        ref_orf_length = ref_orf_len, comp_orf_length = NA_integer_,
        comp_stop_tx_pos = NA_integer_, n_downstream_ejc = NA_integer_,
        category = "mapping_failed")
      next
    }

    comp_orf_len <- comp_stop - comp_tx

    # Count downstream EJCs
    junctions <- .getJunctionPositions(comp_s$starts, comp_s$ends, strand)
    stop_end  <- comp_stop + 2L
    n_dejc    <- sum(junctions > (stop_end + ejc_threshold - 1L))

    # Classify
    if (n_dejc > 0L) {
      cat_val <- "effectively_ptc"
    } else if (!is.na(ref_orf_len) && comp_orf_len < ref_orf_len) {
      cat_val <- "truncated_no_ejc"
    } else {
      cat_val <- "no_downstream_ejc"
    }

    results[[i]] <- tibble::tibble(
      reference_isoform_id = ref_id, comparator_isoform_id = comp_id,
      ref_atg_genomic = as.numeric(ref_atg_g),
      ref_atg_exonic_in_comp = TRUE,
      ref_orf_length = ref_orf_len, comp_orf_length = comp_orf_len,
      comp_stop_tx_pos = comp_stop, n_downstream_ejc = n_dejc,
      alt_start_tx_pos = NA_integer_, alt_start_orf_length = NA_integer_,
      category = cat_val)
  }

  dplyr::bind_rows(results)
}


# ---------- Internal helpers ----------

#' Check if a 3-nucleotide codon is exonic
#' @keywords internal
.isCodonExonic <- function(pos, exon_starts, exon_ends, strand) {
  if (is.na(pos)) return(FALSE)
  if (strand == "+") {
    positions <- c(pos, pos + 1L, pos + 2L)
  } else {
    positions <- c(pos, pos - 1L, pos - 2L)
  }
  all(vapply(positions, function(p) {
    any(p >= exon_starts & p <= exon_ends)
  }, logical(1)))
}

#' Find first in-frame stop codon in a sequence
#' @keywords internal
.findFirstStop <- function(seq_str, atg_pos, stop_codons) {
  pos <- atg_pos + 3L
  while (pos + 2L <= nchar(seq_str)) {
    if (substr(seq_str, pos, pos + 2L) %in% stop_codons) return(pos)
    pos <- pos + 3L
  }
  NA_integer_
}

#' Find the first viable ATG in a transcript sequence
#'
#' Scans \code{seq_str} for every ATG codon, follows each to the first
#' in-frame stop, and returns the first ATG whose ORF length (in nt,
#' ATG inclusive -> first stop exclusive) meets or exceeds
#' \code{min_orf_nt}. When no such ATG is found, returns a list with
#' \code{atg_tx = NA}.
#'
#' @param seq_str Uppercase DNA sequence in transcript orientation.
#' @param stop_codons Character vector of stop codons
#'   (usually \code{c("TAA","TAG","TGA")}).
#' @param min_orf_nt Integer; minimum ORF length (nt) to treat an ATG
#'   as viable. Default 30 (10 aa).
#' @return list(atg_tx, stop_tx, orf_length) with NA values when no ATG
#'   passes the threshold.
#' @keywords internal
.findFirstViableATG <- function(seq_str, stop_codons, min_orf_nt = 30L) {
  n <- nchar(seq_str)
  if (n < 3L) return(list(atg_tx = NA_integer_, stop_tx = NA_integer_,
                          orf_length = NA_integer_))
  # Find all ATG positions in the transcript (1-based). gregexpr returns
  # a vector; positions are ordered 5' -> 3'.
  m <- gregexpr("ATG", seq_str, fixed = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L)
    return(list(atg_tx = NA_integer_, stop_tx = NA_integer_,
                orf_length = NA_integer_))
  for (atg in m) {
    stop_p <- .findFirstStop(seq_str, as.integer(atg), stop_codons)
    if (!is.na(stop_p)) {
      orf_len <- stop_p - as.integer(atg)
      if (orf_len >= min_orf_nt)
        return(list(atg_tx = as.integer(atg), stop_tx = as.integer(stop_p),
                    orf_length = as.integer(orf_len)))
    }
  }
  list(atg_tx = NA_integer_, stop_tx = NA_integer_, orf_length = NA_integer_)
}

#' Get junction positions in transcript space
#' @keywords internal
.getJunctionPositions <- function(exon_starts, exon_ends, strand) {
  n <- length(exon_starts)
  if (n <= 1L) return(integer(0))
  exon_lengths <- exon_ends - exon_starts + 1L
  if (strand == "-") exon_lengths <- rev(exon_lengths)
  cumsum(exon_lengths)[-n]
}
