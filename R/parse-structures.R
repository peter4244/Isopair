#' Parse Isoform Structures from a GTF File
#'
#' Extracts exon coordinates for isoforms from a GTF annotation file and
#' returns a nested tibble with one row per isoform. Optionally filters to
#' a specific set of isoform IDs.
#'
#' @param gtf_path Path to a GTF file (can be gzipped).
#' @param isoform_ids Optional character vector of isoform IDs to keep.
#'   If NULL (default), all isoforms in the GTF are returned.
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A tibble with columns: isoform_id, gene_id, chr, strand,
#'   n_exons, exon_starts (list), exon_ends (list), tx_start, tx_end,
#'   n_junctions.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' head(structures)
#' @export
#' @importFrom rtracklayer import
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr filter select group_by summarise mutate rename n first
#'   arrange n_distinct
#' @importFrom rlang .data
parseIsoformStructures <- function(gtf_path, isoform_ids = NULL,
                                   verbose = TRUE) {
  if (!file.exists(gtf_path)) {
    stop(sprintf("GTF file not found: %s", gtf_path))
  }

  if (verbose) message("Loading GTF: ", gtf_path)
  gtf <- rtracklayer::import(gtf_path)

  if (verbose) message("Filtering to exon features...")
  exons <- gtf |>
    tibble::as_tibble() |>
    dplyr::filter(.data$type == "exon") |>
    dplyr::select("transcript_id", "gene_id", "seqnames", "start", "end",
                  "strand")

  if (!is.null(isoform_ids)) {
    exons <- dplyr::filter(exons, .data$transcript_id %in% isoform_ids)
  }

  if (nrow(exons) == 0L) {
    warning("No exons found after filtering")
    return(tibble::tibble(
      isoform_id = character(0), gene_id = character(0),
      chr = character(0), strand = character(0),
      n_exons = integer(0),
      exon_starts = list(), exon_ends = list(),
      tx_start = integer(0), tx_end = integer(0),
      n_junctions = integer(0)
    ))
  }

  if (verbose) {
    message(sprintf("  %d exons for %d isoforms",
                    nrow(exons), dplyr::n_distinct(exons$transcript_id)))
  }

  # Build nested structure
  structures <- exons |>
    dplyr::group_by(.data$transcript_id) |>
    dplyr::arrange(.data$start) |>
    dplyr::summarise(
      gene_id = dplyr::first(.data$gene_id),
      chr = dplyr::first(as.character(.data$seqnames)),
      strand = dplyr::first(as.character(.data$strand)),
      n_exons = dplyr::n(),
      exon_starts = list(.data$start),
      exon_ends = list(.data$end),
      tx_start = min(.data$start),
      tx_end = max(.data$end),
      .groups = "drop"
    ) |>
    dplyr::rename(isoform_id = "transcript_id") |>
    dplyr::mutate(n_junctions = .data$n_exons - 1L)

  if (verbose) {
    message(sprintf("  Created structures for %d isoforms", nrow(structures)))
  }

  structures
}
