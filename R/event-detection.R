# Event Detection Engine
#
# CRITICAL: Event Direction Semantics (GAIN vs LOSS)
#
# GAIN and LOSS are defined from the COMPARATOR'S PERSPECTIVE
#
# LOSS = Comparator LOST sequence (reference has MORE)
#   Reconstruction: ADD regions to comparator to build reference
#
# GAIN = Comparator GAINED sequence (reference has LESS)
#   Reconstruction: REMOVE regions from comparator to build reference
#
# 12 event types in 4 categories:
#   Terminal: Alt_TSS, Alt_TES
#   Boundary: A5SS, A3SS, Partial_IR_5, Partial_IR_3
#   Exon-level: SE, Missing_Internal
#   Intron retention: IR, IR_diff_5, IR_diff_3, IR_diff_5_3

#' Detect Structural Events Between Two Isoforms
#'
#' Runs the full hierarchical event detection pipeline for a single
#' reference/comparator pair:
#'
#' 1. Boundary determination (TSS/TES positions, genomic spans)
#' 2. Within-boundary events (only when isoforms overlap):
#'    IR -> Boundary (A5SS/A3SS/Partial_IR) -> SE/Missing_Internal
#' 3. Terminal events (Alt_TSS/Alt_TES, always runs)
#'
#' @param reference_exons Data frame of reference isoform exons with columns:
#'   chr, exon_start, exon_end, strand, gene_id, transcript_id.
#' @param comparator_exons Data frame of comparator isoform exons (same schema).
#' @param gene_id Gene identifier string.
#' @param reference_id Reference isoform ID.
#' @param comparator_id Comparator isoform ID.
#' @param strand Gene strand ("+" or "-").
#' @param tss_tolerance TSS change detection tolerance in bp (default 20).
#' @param tes_tolerance TES change detection tolerance in bp (default 20).
#' @param splice_site_threshold Boundary difference below this is a splice site
#'   event (A5SS/A3SS); at or above is partial intron retention (default 100).
#' @return A tibble of detected events with columns: gene_id,
#'   reference_isoform_id, comparator_isoform_id, event_type, direction, chr,
#'   five_prime, three_prime, strand, bp_diff, missing_terminal_exons,
#'   orphan_terminal_exons, ir_split_exons, ref_junctions, comp_junctions.
#'   Returns a zero-row tibble if no events detected.
#' @examples
#' # Two isoforms: reference has 3 exons, comparator skips the middle one
#' ref_exons <- data.frame(
#'   chr = "chr1", exon_start = c(100L, 300L, 500L),
#'   exon_end = c(200L, 400L, 600L), strand = "+",
#'   gene_id = "gene1", transcript_id = "tx_ref"
#' )
#' comp_exons <- data.frame(
#'   chr = "chr1", exon_start = c(100L, 500L),
#'   exon_end = c(200L, 600L), strand = "+",
#'   gene_id = "gene1", transcript_id = "tx_comp"
#' )
#' events <- detectEvents(ref_exons, comp_exons,
#'   gene_id = "gene1", reference_id = "tx_ref",
#'   comparator_id = "tx_comp", strand = "+"
#' )
#' events$event_type
#' @export
#' @importFrom dplyr filter bind_rows
#' @importFrom tibble tibble
#' @importFrom rlang .data
detectEvents <- function(reference_exons, comparator_exons,
                         gene_id, reference_id, comparator_id, strand,
                         tss_tolerance = 20L, tes_tolerance = 20L,
                         splice_site_threshold = 100L) {

  # Order exons biologically (TSS -> TES)
  ref_ordered <- .orderExonsBiological(reference_exons, strand)
  comp_ordered <- .orderExonsBiological(comparator_exons, strand)

  # Pre-compute splice junctions
  ref_junctions_all <- .computeJunctions(ref_ordered)
  comp_junctions_all <- .computeJunctions(comp_ordered)

  # =========================================================================
  # STEP 1: Boundary Determination
  # =========================================================================
  first_ref <- ref_ordered[1L, ]
  first_comp <- comp_ordered[1L, ]
  last_ref <- ref_ordered[nrow(ref_ordered), ]
  last_comp <- comp_ordered[nrow(comp_ordered), ]

  comp_genomic_start <- min(comp_ordered$exon_start)
  comp_genomic_end <- max(comp_ordered$exon_end)
  ref_genomic_start <- min(ref_ordered$exon_start)
  ref_genomic_end <- max(ref_ordered$exon_end)

  all_events <- list()

  # =========================================================================
  # STEP 2: Within-Boundary Events (skipped if isoforms do not overlap)
  # =========================================================================
  overlap_check <- .checkIsoformOverlap(comp_ordered, ref_ordered)

  if (overlap_check$has_overlap) {

    # -----------------------------------------------------------------------
    # STEP 2a: IR Detection
    # -----------------------------------------------------------------------
    ir_exon_pairs <- list()

    # IR GAIN: comparator exon spans multiple reference exons
    for (i in seq_len(nrow(comp_ordered))) {
      if (.detectIrSimple(comp_ordered[i, ], ref_ordered)) {
        comp_exon <- comp_ordered[i, ]

        overlapping_regions <- list()
        overlapping_ref_indices <- integer()
        first_ref_overlap <- NULL
        last_ref_overlap <- NULL

        for (j in seq_len(nrow(ref_ordered))) {
          ref_exon <- ref_ordered[j, ]
          overlaps <- (ref_exon$exon_start <= comp_exon$exon_end) &&
                      (ref_exon$exon_end >= comp_exon$exon_start)
          if (overlaps) {
            overlapping_regions[[length(overlapping_regions) + 1L]] <-
              sprintf("%d-%d", ref_exon$exon_start, ref_exon$exon_end)
            overlapping_ref_indices <- c(overlapping_ref_indices, j)
            if (is.null(first_ref_overlap)) first_ref_overlap <- ref_exon
            last_ref_overlap <- ref_exon
          }
        }

        # Classify IR_diff type
        ir_event_type <- "IR"
        if (!is.null(first_ref_overlap)) {
          if (strand == "+") {
            diff_5 <- comp_exon$exon_start > first_ref_overlap$exon_start
            diff_3 <- comp_exon$exon_end < last_ref_overlap$exon_end
          } else {
            diff_5 <- comp_exon$exon_end < first_ref_overlap$exon_end
            diff_3 <- comp_exon$exon_start > last_ref_overlap$exon_start
          }
          ir_event_type <- if (diff_5 && diff_3) "IR_diff_5_3" else
                           if (diff_5) "IR_diff_5" else
                           if (diff_3) "IR_diff_3" else "IR"
        }

        ir_range_start <- min(vapply(overlapping_ref_indices,
                                     function(k) ref_ordered$exon_start[k],
                                     numeric(1)))
        ir_range_end <- max(vapply(overlapping_ref_indices,
                                   function(k) ref_ordered$exon_end[k],
                                   numeric(1)))
        ref_jxns_ir <- .junctionsWithinRange(ref_junctions_all,
                                             ir_range_start, ir_range_end)

        all_events[[length(all_events) + 1L]] <- tibble::tibble(
          gene_id = gene_id,
          reference_isoform_id = reference_id,
          comparator_isoform_id = comparator_id,
          event_type = ir_event_type,
          direction = "GAIN",
          chr = comp_exon$chr,
          five_prime = if (strand == "+") comp_exon$exon_start else comp_exon$exon_end,
          three_prime = if (strand == "+") comp_exon$exon_end else comp_exon$exon_start,
          strand = strand,
          bp_diff = NA_integer_,
          missing_terminal_exons = "",
          orphan_terminal_exons = "",
          ir_split_exons = paste(overlapping_regions, collapse = ", "),
          ref_junctions = .formatJunctions(ref_jxns_ir),
          comp_junctions = ""
        )

        for (ref_idx in overlapping_ref_indices)
          ir_exon_pairs[[length(ir_exon_pairs) + 1L]] <-
            list(comp = i, dom = ref_idx)
      }
    }

    # IR LOSS: reference exon spans multiple comparator exons
    for (i in seq_len(nrow(ref_ordered))) {
      if (.detectIrSimple(ref_ordered[i, ], comp_ordered)) {
        ref_exon <- ref_ordered[i, ]

        overlapping_regions <- list()
        overlapping_comp_indices <- integer()
        first_comp_overlap <- NULL
        last_comp_overlap <- NULL

        for (j in seq_len(nrow(comp_ordered))) {
          comp_exon <- comp_ordered[j, ]
          overlaps <- (comp_exon$exon_start <= ref_exon$exon_end) &&
                      (comp_exon$exon_end >= ref_exon$exon_start)
          if (overlaps) {
            overlapping_regions[[length(overlapping_regions) + 1L]] <-
              sprintf("%d-%d",
                      max(comp_exon$exon_start, ref_exon$exon_start),
                      min(comp_exon$exon_end, ref_exon$exon_end))
            overlapping_comp_indices <- c(overlapping_comp_indices, j)
            if (is.null(first_comp_overlap)) first_comp_overlap <- comp_exon
            last_comp_overlap <- comp_exon
          }
        }

        ir_event_type <- "IR"
        if (!is.null(first_comp_overlap)) {
          if (strand == "+") {
            diff_5 <- ref_exon$exon_start < first_comp_overlap$exon_start
            diff_3 <- ref_exon$exon_end > last_comp_overlap$exon_end
          } else {
            diff_5 <- ref_exon$exon_end > first_comp_overlap$exon_end
            diff_3 <- ref_exon$exon_start < last_comp_overlap$exon_start
          }
          ir_event_type <- if (diff_5 && diff_3) "IR_diff_5_3" else
                           if (diff_5) "IR_diff_5" else
                           if (diff_3) "IR_diff_3" else "IR"
        }

        comp_jxns_ir <- .junctionsWithinRange(comp_junctions_all,
                                              ref_exon$exon_start,
                                              ref_exon$exon_end)

        all_events[[length(all_events) + 1L]] <- tibble::tibble(
          gene_id = gene_id,
          reference_isoform_id = reference_id,
          comparator_isoform_id = comparator_id,
          event_type = ir_event_type,
          direction = "LOSS",
          chr = ref_exon$chr,
          five_prime = if (strand == "+") ref_exon$exon_start else ref_exon$exon_end,
          three_prime = if (strand == "+") ref_exon$exon_end else ref_exon$exon_start,
          strand = strand,
          bp_diff = NA_integer_,
          missing_terminal_exons = "",
          orphan_terminal_exons = "",
          ir_split_exons = paste(overlapping_regions, collapse = ", "),
          ref_junctions = "",
          comp_junctions = .formatJunctions(comp_jxns_ir)
        )

        for (comp_idx in overlapping_comp_indices)
          ir_exon_pairs[[length(ir_exon_pairs) + 1L]] <-
            list(comp = comp_idx, dom = i)
      }
    }

    # -----------------------------------------------------------------------
    # STEP 2b: Boundary Events (non-IR overlapping exon pairs)
    # -----------------------------------------------------------------------
    for (i in seq_len(nrow(comp_ordered))) {
      comp_exon <- comp_ordered[i, ]
      is_first_comp <- (i == 1L)
      is_last_comp <- (i == nrow(comp_ordered))

      for (j in seq_len(nrow(ref_ordered))) {
        ref_exon <- ref_ordered[j, ]
        is_first_ref <- (j == 1L)
        is_last_ref <- (j == nrow(ref_ordered))

        overlaps <- (ref_exon$exon_start <= comp_exon$exon_end) &&
                    (ref_exon$exon_end >= comp_exon$exon_start)
        if (!overlaps) next

        has_ir <- any(vapply(ir_exon_pairs, function(pair) {
          pair$comp == i && pair$dom == j
        }, logical(1)))
        if (has_ir) next

        event_result <- .detectSharedBoundaryEvent(
          ref_exon, comp_exon, strand,
          is_first_exon = is_first_ref && is_first_comp,
          is_last_exon = is_last_ref && is_last_comp,
          is_first_exon_ref = is_first_ref,
          is_last_exon_ref = is_last_ref,
          is_first_exon_comp = is_first_comp,
          is_last_exon_comp = is_last_comp,
          splice_site_threshold = splice_site_threshold
        )

        if (event_result$event_type != "none") {
          event_type <- event_result$event_type
          direction <- if (!is.null(event_result$direction)) {
            event_result$direction
          } else {
            "-"
          }

          coords <- .computeBoundaryCoords(event_type, direction, strand,
                                           ref_exon, comp_exon)

          evt_start <- min(coords$five_prime, coords$three_prime)
          evt_end <- max(coords$five_prime, coords$three_prime)
          ref_jxns_evt <- .junctionsTouchingRange(ref_junctions_all,
                                                  evt_start, evt_end)
          comp_jxns_evt <- .junctionsTouchingRange(comp_junctions_all,
                                                   evt_start, evt_end)

          all_events[[length(all_events) + 1L]] <- tibble::tibble(
            gene_id = gene_id,
            reference_isoform_id = reference_id,
            comparator_isoform_id = comparator_id,
            event_type = event_type,
            direction = direction,
            chr = ref_exon$chr,
            five_prime = coords$five_prime,
            three_prime = coords$three_prime,
            strand = strand,
            bp_diff = event_result$bp_diff,
            missing_terminal_exons = "",
            orphan_terminal_exons = "",
            ir_split_exons = "",
            ref_junctions = .formatJunctions(ref_jxns_evt),
            comp_junctions = .formatJunctions(comp_jxns_evt)
          )

          # Second event for asymmetric dual-boundary pairs
          if (!is.null(event_result$second_event)) {
            event_type2 <- event_result$second_event$event_type
            direction2 <- event_result$second_event$direction

            coords2 <- .computeBoundaryCoords(event_type2, direction2, strand,
                                              ref_exon, comp_exon)

            evt_start2 <- min(coords2$five_prime, coords2$three_prime)
            evt_end2 <- max(coords2$five_prime, coords2$three_prime)
            ref_jxns_evt2 <- .junctionsTouchingRange(ref_junctions_all,
                                                     evt_start2, evt_end2)
            comp_jxns_evt2 <- .junctionsTouchingRange(comp_junctions_all,
                                                      evt_start2, evt_end2)

            all_events[[length(all_events) + 1L]] <- tibble::tibble(
              gene_id = gene_id,
              reference_isoform_id = reference_id,
              comparator_isoform_id = comparator_id,
              event_type = event_type2,
              direction = direction2,
              chr = ref_exon$chr,
              five_prime = coords2$five_prime,
              three_prime = coords2$three_prime,
              strand = strand,
              bp_diff = event_result$second_event$bp_diff,
              missing_terminal_exons = "",
              orphan_terminal_exons = "",
              ir_split_exons = "",
              ref_junctions = .formatJunctions(ref_jxns_evt2),
              comp_junctions = .formatJunctions(comp_jxns_evt2)
            )
          }
        }
      }
    }

    # -----------------------------------------------------------------------
    # STEP 2c: SE / Missing_Internal Detection
    # -----------------------------------------------------------------------

    # LOSS: ref exon absent in comp
    for (i in seq_len(nrow(ref_ordered))) {
      ref_exon <- ref_ordered[i, ]

      if (ref_exon$exon_start < comp_genomic_start ||
          ref_exon$exon_end > comp_genomic_end) next
      if (.exonOverlapsAny(ref_exon, comp_ordered)) next

      has_prev <- i > 1L
      has_next <- i < nrow(ref_ordered)
      flanking_ok <- has_prev && has_next &&
                     .exonOverlapsAny(ref_ordered[i - 1L, ], comp_ordered) &&
                     .exonOverlapsAny(ref_ordered[i + 1L, ], comp_ordered)
      etype <- if (flanking_ok) "SE" else "Missing_Internal"

      ref_jxns_se <- .junctionsTouchingRange(ref_junctions_all,
                                             ref_exon$exon_start,
                                             ref_exon$exon_end)
      comp_jxns_se <- .junctionsSpanningRange(comp_junctions_all,
                                              ref_exon$exon_start,
                                              ref_exon$exon_end)

      all_events[[length(all_events) + 1L]] <- tibble::tibble(
        gene_id = gene_id,
        reference_isoform_id = reference_id,
        comparator_isoform_id = comparator_id,
        event_type = etype,
        direction = "LOSS",
        chr = ref_exon$chr,
        five_prime = if (strand == "+") ref_exon$exon_start else ref_exon$exon_end,
        three_prime = if (strand == "+") ref_exon$exon_end else ref_exon$exon_start,
        strand = strand,
        bp_diff = ref_exon$exon_end - ref_exon$exon_start + 1L,
        missing_terminal_exons = "",
        orphan_terminal_exons = "",
        ir_split_exons = "",
        ref_junctions = .formatJunctions(ref_jxns_se),
        comp_junctions = .formatJunctions(comp_jxns_se)
      )
    }

    # GAIN: comp exon absent in ref
    for (i in seq_len(nrow(comp_ordered))) {
      comp_exon <- comp_ordered[i, ]

      if (comp_exon$exon_start < ref_genomic_start ||
          comp_exon$exon_end > ref_genomic_end) next
      if (.exonOverlapsAny(comp_exon, ref_ordered)) next
      if (i == 1L || i == nrow(comp_ordered)) next

      flanking_ok <- .exonOverlapsAny(comp_ordered[i - 1L, ], ref_ordered) &&
                     .exonOverlapsAny(comp_ordered[i + 1L, ], ref_ordered)
      etype <- if (flanking_ok) "SE" else "Missing_Internal"

      ref_jxns_se <- .junctionsSpanningRange(ref_junctions_all,
                                             comp_exon$exon_start,
                                             comp_exon$exon_end)
      comp_jxns_se <- .junctionsTouchingRange(comp_junctions_all,
                                              comp_exon$exon_start,
                                              comp_exon$exon_end)

      all_events[[length(all_events) + 1L]] <- tibble::tibble(
        gene_id = gene_id,
        reference_isoform_id = reference_id,
        comparator_isoform_id = comparator_id,
        event_type = etype,
        direction = "GAIN",
        chr = comp_exon$chr,
        five_prime = if (strand == "+") comp_exon$exon_start else comp_exon$exon_end,
        three_prime = if (strand == "+") comp_exon$exon_end else comp_exon$exon_start,
        strand = strand,
        bp_diff = comp_exon$exon_end - comp_exon$exon_start + 1L,
        missing_terminal_exons = "",
        orphan_terminal_exons = "",
        ir_split_exons = "",
        ref_junctions = .formatJunctions(ref_jxns_se),
        comp_junctions = .formatJunctions(comp_jxns_se)
      )
    }

  } else {
    # Non-overlapping: gap zone Missing_Internal LOSS events
    for (i in seq_len(nrow(ref_ordered))) {
      ref_exon <- ref_ordered[i, ]
      if (ref_exon$exon_start < comp_genomic_start ||
          ref_exon$exon_end > comp_genomic_end) next
      if (.exonOverlapsAny(ref_exon, comp_ordered)) next

      ref_jxns_gap <- .junctionsTouchingRange(ref_junctions_all,
                                              ref_exon$exon_start,
                                              ref_exon$exon_end)
      comp_jxns_gap <- .junctionsSpanningRange(comp_junctions_all,
                                               ref_exon$exon_start,
                                               ref_exon$exon_end)

      all_events[[length(all_events) + 1L]] <- tibble::tibble(
        gene_id = gene_id,
        reference_isoform_id = reference_id,
        comparator_isoform_id = comparator_id,
        event_type = "Missing_Internal",
        direction = "LOSS",
        chr = ref_exon$chr,
        five_prime = if (strand == "+") ref_exon$exon_start else ref_exon$exon_end,
        three_prime = if (strand == "+") ref_exon$exon_end else ref_exon$exon_start,
        strand = strand,
        bp_diff = ref_exon$exon_end - ref_exon$exon_start + 1L,
        missing_terminal_exons = "",
        orphan_terminal_exons = "",
        ir_split_exons = "",
        ref_junctions = .formatJunctions(ref_jxns_gap),
        comp_junctions = .formatJunctions(comp_jxns_gap)
      )
    }
  }

  # =========================================================================
  # STEP 3: Terminal Events (always runs)
  # =========================================================================

  # Alt_TSS
  if (.detectTssChange(first_ref, first_comp, strand, tss_tolerance)) {
    if (strand == "+") {
      ref_tss <- first_ref$exon_start; comp_tss <- first_comp$exon_start
    } else {
      ref_tss <- first_ref$exon_end; comp_tss <- first_comp$exon_end
    }

    direction <- if ((strand == "+" && ref_tss < comp_tss) ||
                     (strand == "-" && ref_tss > comp_tss)) "LOSS" else "GAIN"

    missing_exons <- if (direction == "LOSS") {
      .computeMissingTerminalExonsTss(ref_ordered, comp_ordered, strand)
    } else {
      .computeMissingTerminalExonsTss(comp_ordered, ref_ordered, strand)
    }

    orphan_exons <- .computeOrphanTerminalExonsTss(comp_ordered, ref_ordered)

    ref_jxns_tss <- if (length(ref_junctions_all) > 0L) {
      ref_junctions_all[1L]
    } else {
      character(0)
    }
    comp_jxns_tss <- if (length(comp_junctions_all) > 0L) {
      comp_junctions_all[1L]
    } else {
      character(0)
    }

    all_events[[length(all_events) + 1L]] <- tibble::tibble(
      gene_id = gene_id,
      reference_isoform_id = reference_id,
      comparator_isoform_id = comparator_id,
      event_type = "Alt_TSS",
      direction = direction,
      chr = first_ref$chr,
      five_prime = ref_tss,
      three_prime = comp_tss,
      strand = strand,
      bp_diff = abs(ref_tss - comp_tss),
      missing_terminal_exons = missing_exons,
      orphan_terminal_exons = orphan_exons,
      ir_split_exons = "",
      ref_junctions = .formatJunctions(ref_jxns_tss),
      comp_junctions = .formatJunctions(comp_jxns_tss)
    )
  }

  # Alt_TES
  if (.detectTesChange(last_ref, last_comp, strand, tes_tolerance)) {
    if (strand == "+") {
      ref_tes <- last_ref$exon_end; comp_tes <- last_comp$exon_end
    } else {
      ref_tes <- last_ref$exon_start; comp_tes <- last_comp$exon_start
    }

    direction <- if ((strand == "+" && ref_tes > comp_tes) ||
                     (strand == "-" && ref_tes < comp_tes)) "LOSS" else "GAIN"

    missing_exons <- if (direction == "LOSS") {
      .computeMissingTerminalExonsTes(ref_ordered, comp_ordered, strand)
    } else {
      .computeMissingTerminalExonsTes(comp_ordered, ref_ordered, strand)
    }

    orphan_exons <- .computeOrphanTerminalExonsTes(comp_ordered, ref_ordered)

    ref_jxns_tes <- if (length(ref_junctions_all) > 0L) {
      ref_junctions_all[length(ref_junctions_all)]
    } else {
      character(0)
    }
    comp_jxns_tes <- if (length(comp_junctions_all) > 0L) {
      comp_junctions_all[length(comp_junctions_all)]
    } else {
      character(0)
    }

    all_events[[length(all_events) + 1L]] <- tibble::tibble(
      gene_id = gene_id,
      reference_isoform_id = reference_id,
      comparator_isoform_id = comparator_id,
      event_type = "Alt_TES",
      direction = direction,
      chr = last_ref$chr,
      five_prime = ref_tes,
      three_prime = comp_tes,
      strand = strand,
      bp_diff = abs(ref_tes - comp_tes),
      missing_terminal_exons = missing_exons,
      orphan_terminal_exons = orphan_exons,
      ir_split_exons = "",
      ref_junctions = .formatJunctions(ref_jxns_tes),
      comp_junctions = .formatJunctions(comp_jxns_tes)
    )
  }

  if (length(all_events) > 0L) return(dplyr::bind_rows(all_events))
  tibble::tibble()
}


# ============================================================================
# Internal Helper Functions
# ============================================================================

#' @keywords internal
.detectTssChange <- function(exon_ref, exon_comp, strand, tolerance) {
  if (strand == "+") {
    tss_ref <- exon_ref$exon_start; tss_comp <- exon_comp$exon_start
  } else {
    tss_ref <- exon_ref$exon_end; tss_comp <- exon_comp$exon_end
  }
  if (abs(tss_ref - tss_comp) > tolerance) return(TRUE)
  first_overlap <- exon_ref$exon_start <= exon_comp$exon_end &&
                   exon_ref$exon_end >= exon_comp$exon_start
  !first_overlap
}

#' @keywords internal
.detectTesChange <- function(exon_ref, exon_comp, strand, tolerance) {
  if (strand == "+") {
    tes_ref <- exon_ref$exon_end; tes_comp <- exon_comp$exon_end
  } else {
    tes_ref <- exon_ref$exon_start; tes_comp <- exon_comp$exon_start
  }
  if (abs(tes_ref - tes_comp) > tolerance) return(TRUE)
  last_overlap <- exon_ref$exon_start <= exon_comp$exon_end &&
                  exon_ref$exon_end >= exon_comp$exon_start
  !last_overlap
}

#' @keywords internal
.computeMissingTerminalExonsTss <- function(ref_exons, comp_exons, strand) {
  if (strand == "+") {
    ref_tss <- ref_exons$exon_start[1L]
    comp_tss <- comp_exons$exon_start[1L]
    if (ref_tss >= comp_tss) return("")
    cutoff <- comp_tss
  } else {
    ref_tss <- ref_exons$exon_end[1L]
    comp_tss <- comp_exons$exon_end[1L]
    if (ref_tss <= comp_tss) return("")
    cutoff <- comp_tss
  }

  missing_ranges <- list()
  for (i in seq_len(nrow(ref_exons))) {
    exon <- ref_exons[i, ]
    if (strand == "+") {
      if (exon$exon_end < cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start < cutoff && exon$exon_end >= cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, cutoff - 1L)
        break
      } else {
        break
      }
    } else {
      if (exon$exon_start > cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start <= cutoff && exon$exon_end > cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", cutoff + 1L, exon$exon_end)
        break
      } else {
        break
      }
    }
  }
  if (length(missing_ranges) > 0L) paste(missing_ranges, collapse = ",") else ""
}

#' @keywords internal
.computeMissingTerminalExonsTes <- function(ref_exons, comp_exons, strand) {
  last_ref_idx <- nrow(ref_exons)
  last_comp_idx <- nrow(comp_exons)

  if (strand == "+") {
    ref_tes <- ref_exons$exon_end[last_ref_idx]
    comp_tes <- comp_exons$exon_end[last_comp_idx]
    if (ref_tes <= comp_tes) return("")
    cutoff <- comp_tes
  } else {
    ref_tes <- ref_exons$exon_start[last_ref_idx]
    comp_tes <- comp_exons$exon_start[last_comp_idx]
    if (ref_tes >= comp_tes) return("")
    cutoff <- comp_tes
  }

  missing_ranges <- list()
  for (i in last_ref_idx:1L) {
    exon <- ref_exons[i, ]
    if (strand == "+") {
      if (exon$exon_start > cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_start <= cutoff && exon$exon_end > cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", cutoff + 1L, exon$exon_end)
        break
      } else {
        break
      }
    } else {
      if (exon$exon_end < cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, exon$exon_end)
      } else if (exon$exon_end >= cutoff && exon$exon_start < cutoff) {
        missing_ranges[[length(missing_ranges) + 1L]] <-
          sprintf("%d-%d", exon$exon_start, cutoff - 1L)
        break
      } else {
        break
      }
    }
  }

  if (length(missing_ranges) > 0L) {
    paste(rev(missing_ranges), collapse = ",")
  } else {
    ""
  }
}

#' @keywords internal
.computeOrphanTerminalExonsTss <- function(walk_ordered, ref_ordered) {
  orphan_ranges <- list()
  for (i in seq_len(nrow(walk_ordered))) {
    exon <- walk_ordered[i, ]
    has_overlap <- any(vapply(seq_len(nrow(ref_ordered)), function(j) {
      r <- ref_ordered[j, ]
      exon$exon_start <= r$exon_end && exon$exon_end >= r$exon_start
    }, logical(1)))
    if (has_overlap) break
    orphan_ranges[[length(orphan_ranges) + 1L]] <-
      sprintf("%d-%d", exon$exon_start, exon$exon_end)
  }
  if (length(orphan_ranges) > 0L) paste(orphan_ranges, collapse = ",") else ""
}

#' @keywords internal
.computeOrphanTerminalExonsTes <- function(walk_ordered, ref_ordered) {
  orphan_ranges <- list()
  for (i in rev(seq_len(nrow(walk_ordered)))) {
    exon <- walk_ordered[i, ]
    has_overlap <- any(vapply(seq_len(nrow(ref_ordered)), function(j) {
      r <- ref_ordered[j, ]
      exon$exon_start <= r$exon_end && exon$exon_end >= r$exon_start
    }, logical(1)))
    if (has_overlap) break
    orphan_ranges[[length(orphan_ranges) + 1L]] <-
      sprintf("%d-%d", exon$exon_start, exon$exon_end)
  }
  if (length(orphan_ranges) > 0L) {
    paste(rev(orphan_ranges), collapse = ",")
  } else {
    ""
  }
}

#' Detect shared-boundary and dual-boundary splicing events
#' @keywords internal
.detectSharedBoundaryEvent <- function(exon_ref, exon_comp, strand,
                                       is_first_exon = FALSE,
                                       is_last_exon = FALSE,
                                       is_first_exon_ref = FALSE,
                                       is_last_exon_ref = FALSE,
                                       is_first_exon_comp = FALSE,
                                       is_last_exon_comp = FALSE,
                                       splice_site_threshold = 100L) {
  if (strand == "+") {
    shares_acceptor <- (exon_ref$exon_start == exon_comp$exon_start)
    differs_donor <- (exon_ref$exon_end != exon_comp$exon_end)
    donor_diff <- if (differs_donor) abs(exon_ref$exon_end - exon_comp$exon_end) else 0L
    shares_donor <- (exon_ref$exon_end == exon_comp$exon_end)
    differs_acceptor <- (exon_ref$exon_start != exon_comp$exon_start)
    acceptor_diff <- if (differs_acceptor) abs(exon_ref$exon_start - exon_comp$exon_start) else 0L
  } else {
    shares_acceptor <- (exon_ref$exon_end == exon_comp$exon_end)
    differs_donor <- (exon_ref$exon_start != exon_comp$exon_start)
    donor_diff <- if (differs_donor) abs(exon_ref$exon_start - exon_comp$exon_start) else 0L
    shares_donor <- (exon_ref$exon_start == exon_comp$exon_start)
    differs_acceptor <- (exon_ref$exon_end != exon_comp$exon_end)
    acceptor_diff <- if (differs_acceptor) abs(exon_ref$exon_end - exon_comp$exon_end) else 0L
  }

  event_type <- "none"
  bp_diff <- 0L
  direction <- NULL
  second_event <- NULL

  if (shares_acceptor && differs_donor) {
    bp_diff <- donor_diff
    if (strand == "+") {
      direction <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
    } else {
      direction <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
    }
    donor_is_terminal <- is_last_exon
    if (donor_is_terminal) {
      event_type <- "none"
    } else if (bp_diff < splice_site_threshold) {
      event_type <- "A5SS"
    } else {
      event_type <- "Partial_IR_5"
    }

  } else if (shares_donor && differs_acceptor) {
    bp_diff <- acceptor_diff
    if (strand == "+") {
      direction <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
    } else {
      direction <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
    }
    acceptor_is_terminal <- is_first_exon
    if (acceptor_is_terminal) {
      event_type <- "none"
    } else if (bp_diff < splice_site_threshold) {
      event_type <- "A3SS"
    } else {
      event_type <- "Partial_IR_3"
    }

  } else if (differs_acceptor && differs_donor) {
    is_first_any <- is_first_exon || is_first_exon_ref || is_first_exon_comp
    is_last_any <- is_last_exon || is_last_exon_ref || is_last_exon_comp

    if (is_first_any || is_last_any) {
      if (is_first_any) {
        bp_diff <- donor_diff
        event_type <- if (donor_diff < splice_site_threshold) "A5SS" else "Partial_IR_5"
        if (strand == "+") {
          direction <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
        } else {
          direction <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
        }
        is_asymmetric_first <- xor(is_first_exon_ref, is_first_exon_comp)
        if (is_asymmetric_first || is_last_any) {
          acc_type <- if (acceptor_diff < splice_site_threshold) "A3SS" else "Partial_IR_3"
          if (strand == "+") {
            acc_dir <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
          } else {
            acc_dir <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
          }
          second_event <- list(event_type = acc_type, bp_diff = acceptor_diff,
                               direction = acc_dir)
        }
      } else if (is_last_any) {
        bp_diff <- acceptor_diff
        event_type <- if (acceptor_diff < splice_site_threshold) "A3SS" else "Partial_IR_3"
        if (strand == "+") {
          direction <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
        } else {
          direction <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
        }
        is_asymmetric_last <- xor(is_last_exon_ref, is_last_exon_comp)
        if (is_asymmetric_last) {
          don_type <- if (donor_diff < splice_site_threshold) "A5SS" else "Partial_IR_5"
          if (strand == "+") {
            don_dir <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
          } else {
            don_dir <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
          }
          second_event <- list(event_type = don_type, bp_diff = donor_diff,
                               direction = don_dir)
        }
      }
    } else {
      bp_diff <- donor_diff
      event_type <- if (donor_diff < splice_site_threshold) "A5SS" else "Partial_IR_5"
      if (strand == "+") {
        direction <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
      } else {
        direction <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
      }
      acc_type <- if (acceptor_diff < splice_site_threshold) "A3SS" else "Partial_IR_3"
      if (strand == "+") {
        acc_dir <- if (exon_ref$exon_start < exon_comp$exon_start) "LOSS" else "GAIN"
      } else {
        acc_dir <- if (exon_ref$exon_end > exon_comp$exon_end) "LOSS" else "GAIN"
      }
      second_event <- list(event_type = acc_type, bp_diff = acceptor_diff,
                           direction = acc_dir)
    }
  }

  list(event_type = event_type, bp_diff = bp_diff, direction = direction,
       second_event = second_event)
}

#' Compute five_prime/three_prime coordinates for boundary events
#' @keywords internal
.computeBoundaryCoords <- function(event_type, direction, strand,
                                   ref_exon, comp_exon) {
  if (event_type %in% c("A3SS", "Partial_IR_3")) {
    if (strand == "+") {
      if (direction == "LOSS") {
        five_prime <- ref_exon$exon_start; three_prime <- comp_exon$exon_start - 1L
      } else {
        five_prime <- comp_exon$exon_start; three_prime <- ref_exon$exon_start - 1L
      }
    } else {
      if (direction == "LOSS") {
        five_prime <- ref_exon$exon_end; three_prime <- comp_exon$exon_end + 1L
      } else {
        five_prime <- comp_exon$exon_end; three_prime <- ref_exon$exon_end + 1L
      }
    }
  } else if (event_type %in% c("A5SS", "Partial_IR_5")) {
    if (strand == "+") {
      if (direction == "LOSS") {
        five_prime <- comp_exon$exon_end + 1L; three_prime <- ref_exon$exon_end
      } else {
        five_prime <- ref_exon$exon_end + 1L; three_prime <- comp_exon$exon_end
      }
    } else {
      if (direction == "LOSS") {
        five_prime <- comp_exon$exon_start; three_prime <- ref_exon$exon_start - 1L
      } else {
        five_prime <- ref_exon$exon_start - 1L; three_prime <- comp_exon$exon_start
      }
    }
  } else {
    if (strand == "+") {
      five_prime <- min(ref_exon$exon_start, comp_exon$exon_start)
      three_prime <- max(ref_exon$exon_end, comp_exon$exon_end)
    } else {
      five_prime <- max(ref_exon$exon_end, comp_exon$exon_end)
      three_prime <- min(ref_exon$exon_start, comp_exon$exon_start)
    }
  }
  list(five_prime = five_prime, three_prime = three_prime)
}

#' @keywords internal
.detectIrSimple <- function(exon, other_exons) {
  overlapping <- other_exons[
    (other_exons$exon_start >= exon$exon_start &
     other_exons$exon_start <= exon$exon_end) |
    (other_exons$exon_end >= exon$exon_start &
     other_exons$exon_end <= exon$exon_end) |
    (other_exons$exon_start <= exon$exon_start &
     other_exons$exon_end >= exon$exon_end), ]
  nrow(overlapping) >= 2L
}

#' @keywords internal
.exonOverlapsAny <- function(exon, exon_set) {
  any(vapply(seq_len(nrow(exon_set)), function(k) {
    e <- exon_set[k, ]
    exon$exon_start <= e$exon_end && exon$exon_end >= e$exon_start
  }, logical(1)))
}

#' @keywords internal
.checkIsoformOverlap <- function(comp_ordered, ref_ordered) {
  has_overlap <- FALSE
  for (i in seq_len(nrow(comp_ordered))) {
    for (j in seq_len(nrow(ref_ordered))) {
      if (comp_ordered$exon_start[i] <= ref_ordered$exon_end[j] &&
          comp_ordered$exon_end[i] >= ref_ordered$exon_start[j]) {
        has_overlap <- TRUE
        break
      }
    }
    if (has_overlap) break
  }
  ref_exon_coords <- paste(
    sprintf("%d-%d", ref_ordered$exon_start, ref_ordered$exon_end),
    collapse = ", "
  )
  list(has_overlap = has_overlap, dom_exon_coords = ref_exon_coords)
}

#' @keywords internal
.checkSpansFlanking <- function(exon_ref, exon_comp,
                                flanking_exons_ref, flanking_exons_comp) {
  if (nrow(flanking_exons_ref) > 0L) {
    for (i in seq_len(nrow(flanking_exons_ref))) {
      flank <- flanking_exons_ref[i, ]
      if ((exon_comp$exon_start <= flank$exon_end) &&
          (exon_comp$exon_end >= flank$exon_start)) return(TRUE)
    }
  }
  if (nrow(flanking_exons_comp) > 0L) {
    for (i in seq_len(nrow(flanking_exons_comp))) {
      flank <- flanking_exons_comp[i, ]
      if ((exon_ref$exon_start <= flank$exon_end) &&
          (exon_ref$exon_end >= flank$exon_start)) return(TRUE)
    }
  }
  FALSE
}

# Junction helpers
#' @keywords internal
.junctionsTouchingRange <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0L) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    (cc[1] >= range_start && cc[1] <= range_end) ||
    (cc[2] >= range_start && cc[2] <= range_end)
  }, logical(1))]
}

#' @keywords internal
.junctionsWithinRange <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0L) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    cc[1] >= range_start && cc[2] <= range_end
  }, logical(1))]
}

#' @keywords internal
.junctionsSpanningRange <- function(jxn_vec, range_start, range_end) {
  if (length(jxn_vec) == 0L) return(character(0))
  jxn_vec[vapply(jxn_vec, function(jxn) {
    cc <- as.integer(strsplit(jxn, ":")[[1]])
    cc[1] <= range_start && cc[2] >= range_end
  }, logical(1))]
}

#' @keywords internal
.formatJunctions <- function(jxn_vec) {
  if (length(jxn_vec) == 0L) return("")
  paste(jxn_vec, collapse = ",")
}
