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


#' Convert Genomic Position to Transcript Position
#'
#' Maps a genomic coordinate to the corresponding 1-based position in
#' transcript (mRNA) space, accounting for exon structure and strand.
#'
#' @param genomic_pos Integer; the genomic coordinate to map.
#' @param exon_starts Integer vector of exon start coordinates (1-based,
#'   inclusive).
#' @param exon_ends Integer vector of exon end coordinates (1-based,
#'   inclusive).
#' @param strand Character; "+" or "-".
#' @return Integer transcript position (1-based), or \code{NA_integer_} if the
#'   position is not exonic.
#' @examples
#' # Plus-strand gene with two exons: [100,200] and [300,400]
#' genomicToTranscript(150, c(100, 300), c(200, 400), "+")  # returns 51
#' genomicToTranscript(350, c(100, 300), c(200, 400), "+")  # returns 152
#'
#' # Minus-strand: transcript reads from high to low genomic coordinate
#' genomicToTranscript(350, c(100, 300), c(200, 400), "-")  # returns 51
#' @export
genomicToTranscript <- function(genomic_pos, exon_starts, exon_ends, strand) {
  if (is.na(genomic_pos)) return(NA_integer_)
  if (strand == "-") {
    exon_starts <- rev(exon_starts)
    exon_ends   <- rev(exon_ends)
  }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    if (strand == "+") {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (genomic_pos - es) + 1L)
    } else {
      if (genomic_pos >= es && genomic_pos <= ee)
        return(cum + (ee - genomic_pos) + 1L)
    }
    cum <- cum + (ee - es + 1L)
  }
  NA_integer_
}


#' Convert Transcript Position to Genomic Position
#'
#' Maps a 1-based transcript (mRNA) position to the corresponding genomic
#' coordinate, accounting for exon structure and strand.
#'
#' @param tx_pos Integer; the 1-based transcript position.
#' @param exon_starts Integer vector of exon start coordinates.
#' @param exon_ends Integer vector of exon end coordinates.
#' @param strand Character; "+" or "-".
#' @return Integer genomic position, or \code{NA_integer_} if the position
#'   exceeds the transcript length.
#' @examples
#' # Plus-strand gene with two exons: [100,200] and [300,400]
#' transcriptToGenomic(51, c(100, 300), c(200, 400), "+")  # returns 150
#' transcriptToGenomic(152, c(100, 300), c(200, 400), "+")  # returns 350
#' @export
transcriptToGenomic <- function(tx_pos, exon_starts, exon_ends, strand) {
  if (is.na(tx_pos)) return(NA_integer_)
  if (strand == "-") {
    exon_starts <- rev(exon_starts)
    exon_ends   <- rev(exon_ends)
  }
  cum <- 0L
  for (j in seq_along(exon_starts)) {
    es <- exon_starts[j]; ee <- exon_ends[j]
    exon_len <- ee - es + 1L
    if (tx_pos <= cum + exon_len) {
      offset <- tx_pos - cum - 1L
      if (strand == "+") return(es + offset)
      else return(ee - offset)
    }
    cum <- cum + exon_len
  }
  NA_integer_
}


#' Get Strand-Aware Stop Codon Position
#'
#' Returns the strand-correct genomic position of the stop codon for each
#' isoform. On + strand the stop is at \code{cds_stop} (larger coordinate);
#' on - strand it is at \code{cds_start} (smaller coordinate).
#'
#' @param isoform_ids Character vector of isoform IDs.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with columns: \code{isoform_id}, \code{cds_start}, \code{cds_stop},
#'   \code{strand}.
#' @return Named numeric vector of stop codon genomic positions.
#' @examples
#' # Create minimal CDS metadata
#' cds <- data.frame(
#'   isoform_id = c("tx1", "tx2"),
#'   cds_start = c(100, 500), cds_stop = c(400, 800),
#'   strand = c("+", "-"), coding_status = "coding")
#' getStrandAwareStop(c("tx1", "tx2"), cds)  # returns c(tx1=400, tx2=500)
#' @export
getStrandAwareStop <- function(isoform_ids, cds_metadata) {
  idx <- match(isoform_ids, cds_metadata$isoform_id)
  stops <- ifelse(cds_metadata$strand[idx] == "-",
                   cds_metadata$cds_start[idx],
                   cds_metadata$cds_stop[idx])
  stats::setNames(stops, isoform_ids)
}


#' Compute Exonic UTR Lengths
#'
#' Computes the exonic 5'UTR and 3'UTR lengths for coding isoforms from their
#' CDS boundaries and exon structure.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @param isoform_ids Optional character vector of isoform IDs to process.
#'   If \code{NULL} (default), processes all coding isoforms.
#' @return A tibble with columns: \code{isoform_id}, \code{utr5_bp},
#'   \code{utr3_bp}, \code{cds_bp}, \code{tx_length}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' utrs <- computeUtrLengths(structures, cds)
#' head(utrs)
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
computeUtrLengths <- function(structures, cds_metadata, isoform_ids = NULL) {

  coding <- cds_metadata[cds_metadata$coding_status == "coding", ]
  if (!is.null(isoform_ids)) {
    coding <- coding[coding$isoform_id %in% isoform_ids, ]
  }

  results <- vector("list", nrow(coding))

  for (i in seq_len(nrow(coding))) {
    iso_id <- coding$isoform_id[i]
    str_idx <- match(iso_id, structures$isoform_id)
    if (is.na(str_idx)) next

    strand <- coding$strand[i]
    cds_start <- coding$cds_start[i]
    cds_stop  <- coding$cds_stop[i]
    exon_s <- structures$exon_starts[[str_idx]]
    exon_e <- structures$exon_ends[[str_idx]]

    # Strand-aware CDS boundaries in transcription direction
    cds5 <- if (strand == "+") cds_start else cds_stop
    cds3 <- if (strand == "+") cds_stop else cds_start

    utr5 <- 0L; utr3 <- 0L; cds_bp <- 0L
    for (j in seq_along(exon_s)) {
      if (strand == "+") {
        # 5'UTR: exonic bp before cds_start
        if (exon_s[j] < cds5) {
          utr5 <- utr5 + min(exon_e[j], cds5 - 1L) - exon_s[j] + 1L
        }
        # 3'UTR: exonic bp after cds_stop
        if (exon_e[j] > cds3) {
          utr3 <- utr3 + exon_e[j] - max(exon_s[j], cds3 + 1L) + 1L
        }
      } else {
        # Minus strand: 5'UTR is at high coordinates, 3'UTR at low
        if (exon_e[j] > cds5) {
          utr5 <- utr5 + exon_e[j] - max(exon_s[j], cds5 + 1L) + 1L
        }
        if (exon_s[j] < cds3) {
          utr3 <- utr3 + min(exon_e[j], cds3 - 1L) - exon_s[j] + 1L
        }
      }
    }

    tx_len <- sum(exon_e - exon_s + 1L)
    cds_bp <- tx_len - utr5 - utr3

    results[[i]] <- tibble::tibble(
      isoform_id = iso_id, utr5_bp = utr5, utr3_bp = utr3,
      cds_bp = cds_bp, tx_length = tx_len)
  }

  dplyr::bind_rows(results)
}


#' Extract All Detailed Events from Profiles
#'
#' Concatenates the \code{detailed_events} list-column from a profiles data
#' frame into a single data frame.
#'
#' @param profiles A data frame from \code{\link{buildProfiles}()} with a
#'   \code{detailed_events} list-column.
#' @param cols Optional character vector of columns to retain. If \code{NULL}
#'   (default), all columns are kept.
#' @return A data frame with all events concatenated.
#' @examples
#' # After building profiles:
#' # all_events <- extractAllEvents(profiles)
#' @export
extractAllEvents <- function(profiles, cols = NULL) {
  events_list <- lapply(seq_len(nrow(profiles)), function(i) {
    de <- profiles$detailed_events[[i]]
    if (is.null(de) || nrow(de) == 0L) return(NULL)
    if (!is.null(cols)) de <- de[, intersect(cols, names(de)), drop = FALSE]
    de
  })
  dplyr::bind_rows(events_list)
}
