#' Plot an Isoform Pair Comparison
#'
#' Visualizes the exon structures of a reference and comparator isoform pair,
#' optionally showing the reconstructed reference and detected splicing events.
#' When CDS metadata is provided, exons are colored by region type (5'UTR,
#' CDS, 3'UTR).
#'
#' @param reference_exons A data frame with `exon_start` and `exon_end`
#'   columns for the reference isoform.
#' @param comparator_exons A data frame with `exon_start` and `exon_end`
#'   columns for the comparator isoform.
#' @param events A data frame from [detectEvents()] (can be 0-row or NULL).
#' @param gene_id Character; gene identifier for the title.
#' @param reference_id Character; reference isoform ID.
#' @param comparator_id Character; comparator isoform ID.
#' @param strand Character; "+" or "-".
#' @param reconstructed_exons Optional data frame with `exon_start` and
#'   `exon_end` columns for the reconstructed reference isoform.
#' @param cds_metadata Optional data frame with columns `isoform_id`,
#'   `cds_start`, `cds_end`. When provided, exons are segmented and
#'   colored by region type.
#' @param show_events Logical; whether to show event annotation brackets
#'   below the comparator track. Default TRUE.
#' @param pair_label Optional character string for the plot title prefix.
#' @return A ggplot object.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ref <- data.frame(exon_start = c(100L, 300L, 500L),
#'     exon_end = c(200L, 400L, 600L))
#'   comp <- data.frame(exon_start = c(100L, 500L),
#'     exon_end = c(200L, 600L))
#'   events <- detectEvents(
#'     cbind(ref, chr = "chr1", strand = "+", gene_id = "g1",
#'       transcript_id = "ref"),
#'     cbind(comp, chr = "chr1", strand = "+", gene_id = "g1",
#'       transcript_id = "comp"),
#'     "g1", "ref", "comp", "+"
#'   )
#'   p <- plotIsoformPair(ref, comp, events, "g1", "ref", "comp", "+",
#'     show_events = FALSE)
#'   class(p) # "gg" "ggplot"
#' }
#' }
#' @export
#' @importFrom dplyr arrange mutate filter bind_rows select lead
plotIsoformPair <- function(reference_exons, comparator_exons, events,
                            gene_id, reference_id, comparator_id, strand,
                            reconstructed_exons = NULL,
                            cds_metadata = NULL,
                            show_events = TRUE,
                            pair_label = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotIsoformPair(). ",
         "Install it with: install.packages('ggplot2')")
  }

  show_recon <- !is.null(reconstructed_exons) && nrow(reconstructed_exons) > 0
  has_cds <- !is.null(cds_metadata) && nrow(cds_metadata) > 0

  # Compatibility: extractCdsAnnotations() returns cds_stop, but

  # .buildTrackData() expects cds_end
  if (has_cds &&
      "cds_stop" %in% names(cds_metadata) &&
      !"cds_end" %in% names(cds_metadata)) {
    cds_metadata$cds_end <- cds_metadata$cds_stop
  }

  # Track y-positions
  ref_y <- if (show_recon) 3 else 2
  recon_y <- 2
  comp_y <- 1

  # --- Build exon segments (with optional CDS coloring) ---
  ref_plot <- .buildTrackData(reference_exons, "Reference", ref_y,
                              reference_id, cds_metadata, has_cds)
  comp_plot <- .buildTrackData(comparator_exons, "Comparator", comp_y,
                               comparator_id, cds_metadata, has_cds)

  if (show_recon) {
    # Reconstructed track: always uniform color (CDS unknown)
    recon_plot <- reconstructed_exons |>
      dplyr::arrange(.data$exon_start) |>
      dplyr::mutate(
        track = "Reconstructed",
        y_pos = recon_y,
        region = "reconstructed"
      )
    all_exons <- dplyr::bind_rows(ref_plot, recon_plot, comp_plot)
  } else {
    all_exons <- dplyr::bind_rows(ref_plot, comp_plot)
  }

  x_min <- min(all_exons$exon_start)
  x_max <- max(all_exons$exon_end)
  x_range <- x_max - x_min
  x_buffer <- max(x_range * 0.05, 50)

  # --- Intron connectors ---
  intron_data <- .makeIntrons(all_exons)

  p <- ggplot2::ggplot()

  # Intron lines with strand direction chevrons
  if (nrow(intron_data) > 0) {
    p <- p + ggplot2::geom_segment(
      data = intron_data,
      ggplot2::aes(x = .data$exon_end, xend = .data$next_start,
                   y = .data$y_pos, yend = .data$y_pos),
      color = "gray50", linewidth = 0.4
    )

    # Add chevrons indicating transcription direction
    chevron_data <- .makeChevrons(intron_data, strand)
    if (nrow(chevron_data) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data = chevron_data,
          ggplot2::aes(x = .data$x, xend = .data$xend,
                       y = .data$y, yend = .data$yend),
          color = "gray50", linewidth = 0.3
        )
    }
  }

  # --- Exon rectangles ---
  if (has_cds) {
    # Color by region
    region_colors <- c(
      "5UTR" = "#5DADE2", "CDS" = "#1E8449", "3UTR" = "#E67E22",
      "non_coding" = "#BDC3C7", "reconstructed" = "#E67E22"
    )
    p <- p +
      ggplot2::geom_rect(
        data = all_exons,
        ggplot2::aes(xmin = .data$exon_start, xmax = .data$exon_end,
                     ymin = .data$y_pos - 0.25, ymax = .data$y_pos + 0.25,
                     fill = .data$region),
        color = "black", linewidth = 0.3
      ) +
      ggplot2::scale_fill_manual(
        values = region_colors,
        name = "Region",
        breaks = c("5UTR", "CDS", "3UTR", "non_coding"),
        labels = c("5'UTR", "CDS", "3'UTR", "Non-coding")
      )
  } else {
    # Color by track
    track_colors <- c(
      "Reference" = "#27AE60",
      "Reconstructed" = "#E67E22",
      "Comparator" = "#3498DB"
    )
    p <- p +
      ggplot2::geom_rect(
        data = all_exons,
        ggplot2::aes(xmin = .data$exon_start, xmax = .data$exon_end,
                     ymin = .data$y_pos - 0.25, ymax = .data$y_pos + 0.25,
                     fill = .data$track),
        color = "black", linewidth = 0.3
      ) +
      ggplot2::scale_fill_manual(values = track_colors, name = NULL)
  }

  # --- Mismatch highlights on reconstructed track ---
  if (show_recon) {
    mismatch_data <- .identifyMismatches(reference_exons, reconstructed_exons,
                                         strand)
    mismatched <- dplyr::filter(mismatch_data, .data$mismatch)
    if (nrow(mismatched) > 0) {
      p <- p +
        ggplot2::geom_rect(
          data = mismatched,
          ggplot2::aes(xmin = .data$exon_start, xmax = .data$exon_end,
                       ymin = recon_y - 0.35, ymax = recon_y + 0.35),
          fill = NA, color = "red", linewidth = 1.2
        )
    }
  }

  # --- Event annotations ---
  if (show_events && !is.null(events) && nrow(events) > 0) {
    event_annots <- lapply(seq_len(nrow(events)), function(i) {
      .eventToAnnotation(events[i, ])
    })
    event_annots <- dplyr::bind_rows(event_annots)
    event_annots <- dplyr::filter(event_annots,
                                  !is.na(.data$x_start), !is.na(.data$x_end))

    if (nrow(event_annots) > 0) {
      event_annots <- event_annots |>
        dplyr::arrange(.data$x_start) |>
        dplyr::mutate(row = dplyr::row_number())

      y_base <- 0.4
      y_step <- 0.3

      for (i in seq_len(nrow(event_annots))) {
        ea <- event_annots[i, ]
        y_level <- y_base - (i - 1) * y_step

        dir_symbol <- if (ea$direction == "LOSS") "\u2193" else "\u2191"
        event_label <- sprintf("%s %s", ea$label, dir_symbol)
        bracket_color <- if (ea$direction == "LOSS") "#C0392B" else "#2980B9"

        p <- p +
          ggplot2::geom_segment(
            ggplot2::aes(x = !!ea$x_start, xend = !!ea$x_end,
                         y = !!y_level, yend = !!y_level),
            color = bracket_color, linewidth = 1.2
          ) +
          ggplot2::geom_segment(
            ggplot2::aes(x = !!ea$x_start, xend = !!ea$x_start,
                         y = !!y_level, yend = !!(y_level + 0.12)),
            color = bracket_color, linewidth = 0.8
          ) +
          ggplot2::geom_segment(
            ggplot2::aes(x = !!ea$x_end, xend = !!ea$x_end,
                         y = !!y_level, yend = !!(y_level + 0.12)),
            color = bracket_color, linewidth = 0.8
          ) +
          ggplot2::annotate(
            "text",
            x = (ea$x_start + ea$x_end) / 2, y = y_level - 0.12,
            label = event_label, size = 2.7, color = bracket_color,
            fontface = "bold"
          )
      }
    }
  }

  # --- Track labels ---
  if (show_recon) {
    label_df <- tibble::tibble(
      y_pos = c(3, 2, 1),
      label = c(sprintf("Reference\n(%s)", reference_id),
                "Reconstructed",
                sprintf("Comparator\n(%s)", comparator_id))
    )
    y_limits <- c(-0.5, 3.7)
  } else {
    label_df <- tibble::tibble(
      y_pos = c(2, 1),
      label = c(sprintf("Reference\n(%s)", reference_id),
                sprintf("Comparator\n(%s)", comparator_id))
    )
    y_limits <- c(-0.5, 2.7)
  }

  # Adjust y limits for many events
  if (show_events && !is.null(events) && nrow(events) > 3) {
    y_limits[1] <- min(y_limits[1], 0.4 - nrow(events) * 0.3 - 0.3)
  }

  p <- p +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x_min - x_buffer * 1.5, y = .data$y_pos,
                   label = .data$label),
      hjust = 1, size = 3.5, fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(breaks = NULL, limits = y_limits) +
    ggplot2::coord_cartesian(
      xlim = c(x_min - x_buffer * 3, x_max + x_buffer), clip = "off"
    ) +
    ggplot2::labs(
      x = sprintf("Genomic position (%s strand)", strand),
      y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(5, 10, 5, 100)
    )

  # --- Title ---
  event_str <- if (!is.null(events) && nrow(events) > 0) {
    paste(events$event_type, collapse = ", ")
  } else {
    "none"
  }

  title_prefix <- if (!is.null(pair_label) && nchar(pair_label) > 0) {
    paste0(pair_label, " ")
  } else {
    ""
  }
  title_text <- sprintf("%s%s", title_prefix, gene_id)
  subtitle_text <- sprintf("Events: %s", event_str)

  p <- p +
    ggplot2::ggtitle(title_text, subtitle = subtitle_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray30")
    )

  p
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Build track data with optional CDS segmentation
#' @keywords internal
.buildTrackData <- function(exons, track_name, y_pos, isoform_id,
                            cds_metadata, has_cds) {
  exons <- dplyr::arrange(exons, .data$exon_start)

  if (has_cds) {
    cds_row <- cds_metadata[cds_metadata$isoform_id == isoform_id, ]
    if (nrow(cds_row) > 0) {
      cds_start <- cds_row$cds_start[1]
      cds_end <- cds_row$cds_end[1]
      strand <- if ("strand" %in% names(cds_row)) cds_row$strand[1] else "+"
      segments <- lapply(seq_len(nrow(exons)), function(i) {
        .segmentExonByCds(exons$exon_start[i], exons$exon_end[i],
                          cds_start, cds_end, strand)
      })
      result <- dplyr::bind_rows(segments)
      result$track <- track_name
      result$y_pos <- y_pos
      return(result)
    }
  }

  # No CDS: uniform region
  exons$track <- track_name
  exons$y_pos <- y_pos
  exons$region <- "non_coding"
  exons
}


#' Segment a single exon by CDS boundaries
#'
#' Returns up to 3 segments (5'UTR, CDS, 3'UTR) with region labels.
#' Strand-aware labeling: on + strand, before CDS = 5'UTR, after = 3'UTR;
#' on - strand, before CDS (lower coords) = 3'UTR, after = 5'UTR.
#'
#' @param exon_start Integer; exon start coordinate.
#' @param exon_end Integer; exon end coordinate.
#' @param cds_start Integer; CDS start coordinate.
#' @param cds_end Integer; CDS end coordinate.
#' @param strand Character; "+" or "-".
#' @return A tibble with exon_start, exon_end, and region columns.
#' @keywords internal
.segmentExonByCds <- function(exon_start, exon_end, cds_start, cds_end,
                              strand) {
  segments <- list()

  # Determine labels based on strand
  if (strand == "+") {
    before_cds <- "5UTR"
    after_cds <- "3UTR"
  } else {
    before_cds <- "3UTR"
    after_cds <- "5UTR"
  }

  # Exon entirely before CDS
  if (exon_end < cds_start) {
    return(tibble::tibble(
      exon_start = exon_start, exon_end = exon_end, region = before_cds))
  }

  # Exon entirely after CDS
  if (exon_start > cds_end) {
    return(tibble::tibble(
      exon_start = exon_start, exon_end = exon_end, region = after_cds))
  }

  # Portion before CDS
  if (exon_start < cds_start) {
    segments[[length(segments) + 1L]] <- tibble::tibble(
      exon_start = exon_start, exon_end = cds_start - 1L,
      region = before_cds)
  }

  # CDS portion
  overlap_start <- max(exon_start, cds_start)
  overlap_end <- min(exon_end, cds_end)
  if (overlap_start <= overlap_end) {
    segments[[length(segments) + 1L]] <- tibble::tibble(
      exon_start = overlap_start, exon_end = overlap_end, region = "CDS")
  }

  # Portion after CDS
  if (exon_end > cds_end) {
    segments[[length(segments) + 1L]] <- tibble::tibble(
      exon_start = cds_end + 1L, exon_end = exon_end, region = after_cds)
  }

  dplyr::bind_rows(segments)
}


#' Build intron connector segments for all tracks
#' @keywords internal
.makeIntrons <- function(all_exons) {
  tracks <- unique(all_exons$track)
  introns <- lapply(tracks, function(tr) {
    track_data <- all_exons[all_exons$track == tr, ]
    track_data <- dplyr::arrange(track_data, .data$exon_start)
    if (nrow(track_data) < 2L) return(NULL)
    # Get unique exon boundaries for intron connectors
    # (CDS segmentation may split one exon into multiple rows)
    exon_bounds <- .collapseSegments(track_data)
    if (nrow(exon_bounds) < 2L) return(NULL)
    tibble::tibble(
      exon_end = exon_bounds$exon_end[-nrow(exon_bounds)],
      next_start = exon_bounds$exon_start[-1L],
      y_pos = exon_bounds$y_pos[1L]
    )
  })
  dplyr::bind_rows(introns)
}


#' Collapse adjacent CDS segments back to original exon boundaries
#' @keywords internal
.collapseSegments <- function(track_data) {
  track_data <- dplyr::arrange(track_data, .data$exon_start)
  if (nrow(track_data) == 0L) return(track_data)

  merged <- list()
  current_start <- track_data$exon_start[1]
  current_end <- track_data$exon_end[1]
  current_y <- track_data$y_pos[1]

  for (i in seq_len(nrow(track_data))[-1L]) {
    if (track_data$exon_start[i] <= current_end + 1L) {
      current_end <- max(current_end, track_data$exon_end[i])
    } else {
      merged[[length(merged) + 1L]] <- tibble::tibble(
        exon_start = current_start, exon_end = current_end,
        y_pos = current_y)
      current_start <- track_data$exon_start[i]
      current_end <- track_data$exon_end[i]
    }
  }
  merged[[length(merged) + 1L]] <- tibble::tibble(
    exon_start = current_start, exon_end = current_end,
    y_pos = current_y)

  dplyr::bind_rows(merged)
}


#' Build chevron (arrowhead) segments indicating strand direction on introns
#'
#' Places small arrowhead marks at regular intervals along each intron line.
#' For "+" strand, arrowheads point right (>); for "-" strand, they point
#' left (<).
#'
#' @param intron_data A tibble with columns `exon_end`, `next_start`, `y_pos`
#'   from [.makeIntrons()].
#' @param strand Character; "+" or "-".
#' @param spacing Numeric; minimum genomic distance between chevrons. Default
#'   is NULL, which auto-calculates as ~5% of the total plotted range.
#' @param size Numeric; half-height of the arrowhead in y-axis units.
#'   Default 0.08.
#' @return A tibble with columns `x`, `xend`, `y`, `yend` for each arrowhead
#'   arm segment.
#' @keywords internal
.makeChevrons <- function(intron_data, strand, spacing = NULL, size = 0.08) {
  if (nrow(intron_data) == 0L) {
    return(data.frame(x = numeric(0), xend = numeric(0),
                      y = numeric(0), yend = numeric(0)))
  }

  # Auto-calculate spacing if not provided
  if (is.null(spacing)) {
    total_range <- max(intron_data$next_start) - min(intron_data$exon_end)
    spacing <- max(total_range * 0.03, 200)
  }

  # Width of each arrowhead arm in genomic coordinates
  arm_width <- spacing * 0.3

  # Collect all chevron segments across all introns
  all_x <- numeric(0)
  all_xend <- numeric(0)
  all_y <- numeric(0)
  all_yend <- numeric(0)

  for (i in seq_len(nrow(intron_data))) {
    intron_start <- intron_data$exon_end[i]
    intron_end <- intron_data$next_start[i]
    yv <- intron_data$y_pos[i]
    intron_len <- intron_end - intron_start

    # Skip very short introns
    if (intron_len < arm_width * 3) next

    # Place chevrons at even intervals within the intron
    n_chev <- max(1L, floor(intron_len / spacing))
    step <- intron_len / (n_chev + 1)
    positions <- intron_start + step * seq_len(n_chev)

    for (cx in positions) {
      if (strand == "+") {
        # Right-pointing arrowhead: > shape
        all_x    <- c(all_x,    cx - arm_width, cx - arm_width)
        all_xend <- c(all_xend, cx,             cx)
      } else {
        # Left-pointing arrowhead: < shape
        all_x    <- c(all_x,    cx + arm_width, cx + arm_width)
        all_xend <- c(all_xend, cx,             cx)
      }
      all_y    <- c(all_y,    yv + size, yv - size)
      all_yend <- c(all_yend, yv,        yv)
    }
  }

  data.frame(x = all_x, xend = all_xend, y = all_y, yend = all_yend)
}


#' Compare reference and reconstructed exon sets, identify mismatches
#'
#' @param reference_exons Data frame with exon_start, exon_end (sorted).
#' @param reconstructed_exons Data frame with exon_start, exon_end (sorted).
#' @param strand "+" or "-".
#' @param tss_tolerance Integer; tolerance for TSS boundary comparison.
#' @param tes_tolerance Integer; tolerance for TES boundary comparison.
#' @return A tibble with exon_start, exon_end, mismatch (logical), detail.
#' @keywords internal
.identifyMismatches <- function(reference_exons, reconstructed_exons, strand,
                                tss_tolerance = 20L, tes_tolerance = 20L) {
  ref_sorted <- dplyr::arrange(reference_exons, .data$exon_start)
  recon_sorted <- dplyr::arrange(reconstructed_exons, .data$exon_start)

  n_ref <- nrow(ref_sorted)
  n_recon <- nrow(recon_sorted)

  if (n_ref != n_recon) {
    # Exon count mismatch: flag all reconstructed exons
    mismatches <- lapply(seq_len(n_recon), function(i) {
      tibble::tibble(
        exon_start = recon_sorted$exon_start[i],
        exon_end = recon_sorted$exon_end[i],
        mismatch = TRUE,
        detail = sprintf("Exon count: %d ref vs %d recon", n_ref, n_recon)
      )
    })
    return(dplyr::bind_rows(mismatches))
  }

  mismatches <- lapply(seq_len(n_ref), function(i) {
    is_first <- (i == 1L)
    is_last <- (i == n_ref)

    if (strand == "+") {
      start_tol <- if (is_first) tss_tolerance else 0L
      end_tol <- if (is_last) tes_tolerance else 0L
    } else {
      start_tol <- if (is_first) tes_tolerance else 0L
      end_tol <- if (is_last) tss_tolerance else 0L
    }

    start_diff <- abs(ref_sorted$exon_start[i] - recon_sorted$exon_start[i])
    end_diff <- abs(ref_sorted$exon_end[i] - recon_sorted$exon_end[i])

    has_mismatch <- start_diff > start_tol || end_diff > end_tol
    detail <- if (has_mismatch) {
      sprintf("Exon %d: [%d-%d] vs [%d-%d]", i,
              ref_sorted$exon_start[i], ref_sorted$exon_end[i],
              recon_sorted$exon_start[i], recon_sorted$exon_end[i])
    } else {
      ""
    }

    tibble::tibble(
      exon_start = recon_sorted$exon_start[i],
      exon_end = recon_sorted$exon_end[i],
      mismatch = has_mismatch,
      detail = detail
    )
  })

  dplyr::bind_rows(mismatches)
}


#' Extract genomic coordinates from an event row for annotation
#' @param event One-row data frame from [detectEvents()].
#' @return A tibble with x_start, x_end, label, direction.
#' @keywords internal
.eventToAnnotation <- function(event) {
  tibble::tibble(
    x_start = as.numeric(event$five_prime),
    x_end = as.numeric(event$three_prime),
    label = as.character(event$event_type),
    direction = as.character(event$direction)
  )
}
