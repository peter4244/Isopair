# Reconstruction Engine
#
# Reconstructs the reference isoform from comparator exons + detected events.
# Two-phase approach:
#   Phase 1: Internal events (IR, SE, splice sites) -> merge after
#   Phase 2: Terminal events (Alt_TSS, Alt_TES) -> NO merge between events

#' Reconstruct Reference Isoform from Comparator and Events
#'
#' Starting from comparator exon coordinates, applies detected events to
#' reconstruct the reference isoform structure. Events are applied in two
#' phases: internal events first (with merging), then terminal events
#' (without merging between events).
#'
#' @param comparator_exons Data frame of comparator exons with columns:
#'   chr, exon_start, exon_end, strand, gene_id, transcript_id.
#' @param events Tibble of detected events from [detectEvents()].
#' @return Data frame of reconstructed reference exons.
#' @examples
#' # Detect events, then reconstruct reference from comparator
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
#'   "gene1", "tx_ref", "tx_comp", "+"
#' )
#' reconstructed <- reconstructDominant(comp_exons, events)
#' nrow(reconstructed)
#' @export
#' @importFrom dplyr filter arrange distinct bind_rows mutate case_when
#' @importFrom rlang .data
reconstructDominant <- function(comparator_exons, events) {
  reconstructed <- comparator_exons

  # Handle empty events (zero-row tibble or data frame)
  if (nrow(events) == 0L) return(reconstructed)

  # Separate events into internal vs terminal
  internal_events_raw <- events[events$event_type %in%
    c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3",
      "SE", "Missing_Internal", "A5SS", "A3SS",
      "Partial_IR_5", "Partial_IR_3"), ]

  # Filter out Partial_IR events overlapping with IR events
  ir_events <- internal_events_raw[internal_events_raw$event_type %in%
    c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3"), ]
  partial_ir_events <- internal_events_raw[internal_events_raw$event_type %in%
    c("Partial_IR_5", "Partial_IR_3"), ]
  other_events <- internal_events_raw[!internal_events_raw$event_type %in%
    c("IR", "Partial_IR_5", "Partial_IR_3"), ]

  if (nrow(ir_events) > 0L && nrow(partial_ir_events) > 0L) {
    keep_partial <- logical(nrow(partial_ir_events))
    for (i in seq_len(nrow(partial_ir_events))) {
      pe <- partial_ir_events[i, ]
      pe_start <- min(pe$five_prime, pe$three_prime)
      pe_end <- max(pe$five_prime, pe$three_prime)

      overlaps <- FALSE
      for (j in seq_len(nrow(ir_events))) {
        ie <- ir_events[j, ]
        ie_start <- min(ie$five_prime, ie$three_prime)
        ie_end <- max(ie$five_prime, ie$three_prime)
        if (ie_start <= pe_end && ie_end >= pe_start) {
          overlaps <- TRUE
          break
        }
      }
      keep_partial[i] <- !overlaps
    }
    partial_ir_events <- partial_ir_events[keep_partial, ]
  }

  # Combine and order: IR first, then SE, then others
  internal_events <- dplyr::bind_rows(ir_events, partial_ir_events,
                                      other_events)
  if (nrow(internal_events) > 0L) {
    order_priority <- dplyr::case_when(
      internal_events$event_type %in%
        c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3") ~ 1L,
      internal_events$event_type == "SE" ~ 2L,
      TRUE ~ 3L
    )
    internal_events <- internal_events[order(order_priority), ]
  }

  terminal_events <- events[events$event_type %in% c("Alt_TSS", "Alt_TES"), ]
  if (nrow(terminal_events) > 0L) {
    # LOSS before GAIN, then by event type
    terminal_events <- terminal_events[
      order(terminal_events$direction != "LOSS", terminal_events$event_type), ]
  }

  # PHASE 1: Apply internal events
  for (i in seq_len(nrow(internal_events))) {
    event <- internal_events[i, ]
    reconstructed <- tryCatch(
      .applyEvent(reconstructed, event),
      error = function(e) {
        warning(sprintf("Error applying %s event: %s",
                        event$event_type, e$message))
        reconstructed
      }
    )
  }

  # Pre-merge orphan removal
  if (nrow(terminal_events) > 0L) {
    for (i in seq_len(nrow(terminal_events))) {
      te <- terminal_events[i, ]
      if (!is.na(te$orphan_terminal_exons) &&
          te$orphan_terminal_exons != "") {
        ranges_list <- strsplit(te$orphan_terminal_exons, ",")[[1]]
        for (range_str in ranges_list) {
          coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
          reconstructed <- reconstructed[
            !(reconstructed$exon_start >= coords[1] &
              reconstructed$exon_end <= coords[2]), ]
        }
      }
    }
  }

  # Merge after internal events
  reconstructed <- dplyr::arrange(reconstructed, .data$exon_start,
                                  .data$exon_end)
  reconstructed <- dplyr::distinct(reconstructed, .data$exon_start,
                                   .data$exon_end, .keep_all = TRUE)

  if (nrow(reconstructed) > 0L) {
    strand_val <- reconstructed$strand[1L]
    reconstructed <- .mergeAdjacentExons(reconstructed, strand_val)
  }

  # PHASE 2: Apply terminal events WITHOUT merging between events
  for (i in seq_len(nrow(terminal_events))) {
    event <- terminal_events[i, ]
    reconstructed <- tryCatch(
      .applyEvent(reconstructed, event),
      error = function(e) {
        warning(sprintf("Error applying %s event: %s",
                        event$event_type, e$message))
        reconstructed
      }
    )
  }

  # Final cleanup
  reconstructed <- dplyr::arrange(reconstructed, .data$exon_start,
                                  .data$exon_end)
  reconstructed <- dplyr::distinct(reconstructed, .data$exon_start,
                                   .data$exon_end, .keep_all = TRUE)

  if (nrow(reconstructed) > 0L) {
    strand_val <- reconstructed$strand[1L]
    reconstructed <- .mergeAdjacentExons(reconstructed, strand_val)
  }

  reconstructed
}


#' Verify Reconstruction Accuracy
#'
#' Checks whether reconstructed exons match the original reference exons.
#' Internal boundaries require exact match; terminal boundaries (TSS/TES)
#' allow tolerance.
#'
#' @param original_exons Exons from the true reference isoform.
#' @param reconstructed_exons Exons from reconstruction.
#' @param strand Gene strand ("+" or "-").
#' @param tss_tolerance TSS tolerance in bp (default 20).
#' @param tes_tolerance TES tolerance in bp (default 20).
#' @return List with: status ("PASS" or "FAIL"), reason (character),
#'   n_exons_original, n_exons_reconstructed.
#' @examples
#' # Verify that two identical exon sets match
#' exons_a <- data.frame(exon_start = c(100L, 300L), exon_end = c(200L, 400L))
#' exons_b <- data.frame(exon_start = c(100L, 300L), exon_end = c(200L, 400L))
#' result <- verifyReconstruction(exons_a, exons_b, strand = "+")
#' result$status
#' @export
#' @importFrom dplyr arrange
#' @importFrom rlang .data
verifyReconstruction <- function(original_exons, reconstructed_exons, strand,
                                 tss_tolerance = 20L, tes_tolerance = 20L) {
  if (nrow(original_exons) != nrow(reconstructed_exons)) {
    return(list(
      status = "FAIL",
      reason = sprintf("Exon count mismatch: %d original vs %d reconstructed",
                        nrow(original_exons), nrow(reconstructed_exons)),
      n_exons_original = nrow(original_exons),
      n_exons_reconstructed = nrow(reconstructed_exons)
    ))
  }

  original_sorted <- dplyr::arrange(original_exons, .data$exon_start,
                                    .data$exon_end)
  reconstructed_sorted <- dplyr::arrange(reconstructed_exons, .data$exon_start,
                                         .data$exon_end)

  n_exons <- nrow(original_sorted)

  for (i in seq_len(n_exons)) {
    orig <- original_sorted[i, ]
    recon <- reconstructed_sorted[i, ]

    is_first_genomic <- (i == 1L)
    is_last_genomic <- (i == n_exons)

    if (strand == "+") {
      start_tol <- if (is_first_genomic) tss_tolerance else 0L
      end_tol <- if (is_last_genomic) tes_tolerance else 0L
    } else {
      start_tol <- if (is_first_genomic) tes_tolerance else 0L
      end_tol <- if (is_last_genomic) tss_tolerance else 0L
    }

    start_diff <- abs(orig$exon_start - recon$exon_start)
    end_diff <- abs(orig$exon_end - recon$exon_end)

    if (start_diff > start_tol || end_diff > end_tol) {
      return(list(
        status = "FAIL",
        reason = sprintf("Exon %d mismatch: orig [%d-%d] vs recon [%d-%d]",
                          i, orig$exon_start, orig$exon_end,
                          recon$exon_start, recon$exon_end),
        n_exons_original = nrow(original_exons),
        n_exons_reconstructed = nrow(reconstructed_exons)
      ))
    }
  }

  list(
    status = "PASS",
    reason = "Match",
    n_exons_original = nrow(original_exons),
    n_exons_reconstructed = nrow(reconstructed_exons)
  )
}


# ============================================================================
# Internal Event Application Functions
# ============================================================================

#' Apply a single event to exon structure
#' @keywords internal
#' @importFrom dplyr filter arrange bind_rows mutate distinct
#' @importFrom tibble tibble
#' @importFrom rlang .data
.applyEvent <- function(exons, event) {
  event_type <- event$event_type
  direction <- event$direction

  # SE / Missing_Internal
  if (event_type %in% c("SE", "Missing_Internal")) {
    event_start <- min(event$five_prime, event$three_prime)
    event_end <- max(event$five_prime, event$three_prime)

    if (direction == "GAIN") {
      exons <- exons[!(exons$exon_start >= event_start &
                       exons$exon_end <= event_end), ]
    } else {
      new_exon <- tibble::tibble(
        chr = exons$chr[1L], exon_start = event_start, exon_end = event_end,
        strand = exons$strand[1L], gene_id = exons$gene_id[1L],
        transcript_id = exons$transcript_id[1L]
      )
      exons <- dplyr::bind_rows(exons, new_exon) |>
        dplyr::arrange(.data$exon_start)
    }
    return(exons)
  }

  # A5SS / A3SS / Partial_IR
  if (event_type %in% c("A5SS", "A3SS", "Partial_IR_5", "Partial_IR_3")) {
    exons <- .modifyExonBoundary(exons, event)
    return(exons)
  }

  # Alt_TSS / Alt_TES
  if (event_type %in% c("Alt_TSS", "Alt_TES")) {
    exons <- .modifyTerminalExon(exons, event)

    if (!is.na(event$orphan_terminal_exons) &&
        event$orphan_terminal_exons != "") {
      ranges_list <- strsplit(event$orphan_terminal_exons, ",")[[1]]
      for (range_str in ranges_list) {
        coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
        exons <- exons[!(exons$exon_start >= coords[1] &
                         exons$exon_end <= coords[2]), ]
      }
    }
    return(exons)
  }

  # IR events
  if (event_type %in% c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3")) {
    ir_start <- min(event$five_prime, event$three_prime)
    ir_end <- max(event$five_prime, event$three_prime)

    meta_chr <- if (nrow(exons) > 0L) exons$chr[1L] else event$chr
    meta_strand <- if (nrow(exons) > 0L) exons$strand[1L] else event$strand
    meta_gene <- if (nrow(exons) > 0L) exons$gene_id[1L] else event$gene_id
    meta_tx <- if (nrow(exons) > 0L) exons$transcript_id[1L] else NA_character_

    if (direction == "GAIN") {
      for (i in seq_len(nrow(exons))) {
        if (exons$exon_start[i] <= ir_start && exons$exon_end[i] >= ir_end) {
          exons <- exons[-i, ]
          break
        }
      }

      if (!is.na(event$ir_split_exons) && event$ir_split_exons != "") {
        split_ranges <- strsplit(event$ir_split_exons, ",")[[1]]
        for (range_str in split_ranges) {
          coords <- as.integer(strsplit(trimws(range_str), "-")[[1]])
          new_exon <- tibble::tibble(
            chr = meta_chr, exon_start = coords[1], exon_end = coords[2],
            strand = meta_strand, gene_id = meta_gene,
            transcript_id = meta_tx
          )
          exons <- dplyr::bind_rows(exons, new_exon)
        }
        exons <- dplyr::arrange(exons, .data$exon_start)
      }
    } else {
      exons <- exons[!(exons$exon_start <= ir_end &
                       exons$exon_end >= ir_start), ]
      merged_exon <- tibble::tibble(
        chr = meta_chr, exon_start = ir_start, exon_end = ir_end,
        strand = meta_strand, gene_id = meta_gene,
        transcript_id = meta_tx
      )
      exons <- dplyr::bind_rows(exons, merged_exon) |>
        dplyr::arrange(.data$exon_start)
    }
    return(exons)
  }

  exons
}

#' Modify exon boundary for splice site events
#' @keywords internal
.modifyExonBoundary <- function(exons, event) {
  strand_val <- event$strand
  event_type <- event$event_type
  direction <- event$direction
  fp <- event$five_prime
  tp <- event$three_prime

  if (event_type %in% c("A5SS", "Partial_IR_5")) {
    if (strand_val == "+") {
      if (direction == "LOSS") {
        boundary_coord <- min(fp, tp) - 1L
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_end[exon_idx] <- max(fp, tp)
      } else {
        boundary_coord <- max(fp, tp)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_end[exon_idx] <- min(fp, tp) - 1L
      }
    } else {
      if (direction == "LOSS") {
        boundary_coord <- max(fp, tp) + 1L
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_start[exon_idx] <- min(fp, tp) + 1L
      } else {
        boundary_coord <- min(fp, tp)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_start[exon_idx] <- max(fp, tp) + 1L
      }
    }
  } else if (event_type %in% c("A3SS", "Partial_IR_3")) {
    if (strand_val == "+") {
      if (direction == "LOSS") {
        boundary_coord <- max(fp, tp) + 1L
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_start[exon_idx] <- min(fp, tp)
      } else {
        boundary_coord <- min(fp, tp)
        exon_idx <- which(abs(exons$exon_start - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_start[exon_idx] <- max(fp, tp) + 1L
      }
    } else {
      if (direction == "LOSS") {
        boundary_coord <- min(fp, tp) - 1L
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_end[exon_idx] <- max(fp, tp)
      } else {
        boundary_coord <- max(fp, tp)
        exon_idx <- which(abs(exons$exon_end - boundary_coord) <= 1L)[1L]
        if (!is.na(exon_idx)) exons$exon_end[exon_idx] <- min(fp, tp) - 1L
      }
    }
  }

  exons
}

#' Modify terminal exon for Alt_TSS/Alt_TES events
#' @keywords internal
#' @importFrom dplyr arrange distinct filter mutate bind_rows
#' @importFrom tibble tibble
#' @importFrom rlang .data
.modifyTerminalExon <- function(exons, event) {
  five_prime <- event$five_prime

  if (event$direction == "LOSS") {
    if (!is.na(event$missing_terminal_exons) &&
        event$missing_terminal_exons != "") {
      ranges_str <- strsplit(event$missing_terminal_exons, ",")[[1]]
      for (rs in ranges_str) {
        coords <- as.integer(strsplit(trimws(rs), "-")[[1]])
        new_exon <- tibble::tibble(
          chr = exons$chr[1L], exon_start = coords[1], exon_end = coords[2],
          strand = exons$strand[1L], gene_id = exons$gene_id[1L],
          transcript_id = exons$transcript_id[1L]
        )
        exons <- dplyr::bind_rows(exons, new_exon)
      }
      exons <- exons |>
        dplyr::arrange(.data$exon_start) |>
        dplyr::distinct(.data$exon_start, .data$exon_end, .keep_all = TRUE)
    }
  } else {
    remove_high <- (event$event_type == "Alt_TES" && event$strand == "+") ||
                   (event$event_type == "Alt_TSS" && event$strand == "-")

    if (remove_high) {
      exons <- exons[exons$exon_start <= five_prime, ]
      exons$exon_end <- pmin(exons$exon_end, five_prime)
    } else {
      exons <- exons[exons$exon_end >= five_prime, ]
      exons$exon_start <- pmax(exons$exon_start, five_prime)
    }

    exons <- exons[exons$exon_start <= exons$exon_end, ]
  }

  exons
}
