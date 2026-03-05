# Internal utility functions for Isopair
#
# These are shared helper functions used across the package.
# None are exported.

#' Order exons in biological transcription order (5' to 3')
#'
#' Plus strand: ascending by exon_start (TSS has lowest coordinate).
#' Minus strand: descending by exon_start (TSS has highest coordinate).
#'
#' @param exons Data frame with exon_start and exon_end columns.
#' @param strand Character: "+" or "-".
#' @return Exons sorted in biological order with biological_exon_number column.
#' @keywords internal
#' @importFrom dplyr arrange mutate row_number desc
#' @importFrom rlang .data
.orderExonsBiological <- function(exons, strand) {
  if (strand == "+") {
    exons <- exons |>
      dplyr::arrange(.data$exon_start) |>
      dplyr::mutate(biological_exon_number = dplyr::row_number())
  } else if (strand == "-") {
    exons <- exons |>
      dplyr::arrange(dplyr::desc(.data$exon_start)) |>
      dplyr::mutate(biological_exon_number = dplyr::row_number())
  } else {
    stop("Invalid strand: must be '+' or '-'")
  }
  exons
}

#' Compute splice junctions from ordered exons
#'
#' Returns junction strings as "left:right" where left < right (genomic
#' coordinates). Works correctly for both strands because biological ordering
#' means consecutive exons alternate genomic position.
#'
#' @param exons_ordered Exons in biological order (TSS to TES).
#' @return Character vector of "left:right" junction strings.
#' @keywords internal
.computeJunctions <- function(exons_ordered) {
  if (nrow(exons_ordered) < 2L) return(character(0))
  vapply(seq_len(nrow(exons_ordered) - 1L), function(i) {
    e1 <- exons_ordered[i, ]
    e2 <- exons_ordered[i + 1L, ]
    sprintf("%d:%d", min(e1$exon_end, e2$exon_end),
            max(e1$exon_start, e2$exon_start))
  }, character(1))
}

#' Merge adjacent or overlapping exons into continuous segments
#'
#' Sorts by exon_start, then iteratively merges where current exon_end
#' is at least next exon_start - 1 (adjacent or overlapping).
#'
#' @param exons Data frame with exon_start and exon_end columns.
#' @param strand Strand (not used currently, kept for interface consistency).
#' @return Merged exon data frame.
#' @keywords internal
#' @importFrom dplyr arrange bind_rows
#' @importFrom rlang .data
.mergeAdjacentExons <- function(exons, strand) {
  if (nrow(exons) == 0L) return(exons)

  exons <- dplyr::arrange(exons, .data$exon_start)

  merged <- list()
  current <- exons[1L, ]

  for (i in seq_len(nrow(exons))[-1L]) {
    next_exon <- exons[i, ]
    if (current$exon_end >= next_exon$exon_start - 1L) {
      current$exon_end <- max(current$exon_end, next_exon$exon_end)
    } else {
      merged[[length(merged) + 1L]] <- current
      current <- next_exon
    }
  }

  merged[[length(merged) + 1L]] <- current
  dplyr::bind_rows(merged)
}
