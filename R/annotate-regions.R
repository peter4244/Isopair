#' Extract CDS Annotations from a GTF
#'
#' Extracts coding sequence (CDS) coordinates for isoforms from a GTF
#' annotation file. Returns a tibble with CDS boundaries and ORF length
#' for each isoform.
#'
#' @param gtf A character path to a GTF file (can be gzipped) or a
#'   \code{GRanges} object already loaded via \code{rtracklayer::import()}.
#' @param isoform_ids Optional character vector of isoform IDs. If provided,
#'   any IDs not found in the GTF are returned with \code{coding_status =
#'   "unknown"} and NA coordinates.
#' @param verbose Logical; if TRUE (default), print progress messages.
#' @return A tibble with columns: \code{isoform_id}, \code{coding_status}
#'   ("coding", "unknown"), \code{cds_start} (min genomic CDS coordinate),
#'   \code{cds_stop} (max genomic CDS coordinate), \code{orf_length} (total
#'   CDS bp), \code{strand} ("+", "-", or NA for unknown isoforms).
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' cds <- extractCdsAnnotations(gtf_path, verbose = FALSE)
#' head(cds)
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr filter group_by summarise mutate bind_rows n
#' @importFrom rlang .data
extractCdsAnnotations <- function(gtf, isoform_ids = NULL, verbose = TRUE) {

  gtf_tbl <- .loadGtfData(gtf, verbose)

  # Filter to CDS-relevant features
  cds_features <- dplyr::filter(
    gtf_tbl,
    .data$type %in% c("CDS", "start_codon", "stop_codon")
  )

  if (nrow(cds_features) == 0L) {
    if (verbose) message("  No CDS features found in GTF")
    if (!is.null(isoform_ids) && length(isoform_ids) > 0L) {
      return(tibble::tibble(
        isoform_id = isoform_ids,
        coding_status = "unknown",
        cds_start = NA_integer_,
        cds_stop = NA_integer_,
        orf_length = NA_integer_,
        strand = NA_character_
      ))
    }
    return(tibble::tibble(
      isoform_id = character(0),
      coding_status = character(0),
      cds_start = integer(0),
      cds_stop = integer(0),
      orf_length = integer(0),
      strand = character(0)
    ))
  }

  if (verbose) {
    message(sprintf("  %d CDS/codon features for %d isoforms",
                    nrow(cds_features),
                    length(unique(cds_features$transcript_id))))
  }

  # Aggregate per isoform
  cds_only <- dplyr::filter(cds_features, .data$type == "CDS")

  result <- cds_only |>
    dplyr::group_by(.data$transcript_id) |>
    dplyr::summarise(
      coding_status = "coding",
      cds_start = min(.data$start),
      cds_stop = max(.data$end),
      orf_length = sum(.data$end - .data$start + 1L),
      strand = dplyr::first(as.character(.data$strand)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      isoform_id = .data$transcript_id
    )

  # Use min/max across all features (CDS + codons) to capture full CDS range
  all_ranges <- cds_features |>
    dplyr::group_by(.data$transcript_id) |>
    dplyr::summarise(
      full_start = min(.data$start),
      full_stop = max(.data$end),
      .groups = "drop"
    )

  idx <- match(result$transcript_id, all_ranges$transcript_id)
  result$cds_start <- all_ranges$full_start[idx]
  result$cds_stop <- all_ranges$full_stop[idx]

  result <- result[, c("isoform_id", "coding_status", "cds_start",
                        "cds_stop", "orf_length", "strand")]

  # Add missing isoforms as "unknown"
  if (!is.null(isoform_ids)) {
    missing_ids <- setdiff(isoform_ids, result$isoform_id)
    if (length(missing_ids) > 0L) {
      missing_df <- tibble::tibble(
        isoform_id = missing_ids,
        coding_status = "unknown",
        cds_start = NA_integer_,
        cds_stop = NA_integer_,
        orf_length = NA_integer_,
        strand = NA_character_
      )
      result <- dplyr::bind_rows(result, missing_df)
    }
    # Filter to only requested IDs
    result <- result[result$isoform_id %in% isoform_ids, ]
  }

  if (verbose) {
    n_coding <- sum(result$coding_status == "coding")
    message(sprintf("  Result: %d isoforms (%d coding, %d unknown)",
                    nrow(result), n_coding, nrow(result) - n_coding))
  }

  result
}


#' Extract Per-Exon CDS Intervals from a GTF
#'
#' Returns individual CDS exon intervals for each isoform, preserving the
#' per-row structure from the GTF rather than aggregating into min/max as
#' \code{\link{extractCdsAnnotations}()} does. This is needed for reading
#' frame phase computation at specific genomic positions.
#'
#' @param gtf A character path to a GTF file (can be gzipped) or a
#'   \code{GRanges} object already loaded via \code{rtracklayer::import()}.
#' @param isoform_ids Optional character vector of isoform IDs to include.
#' @param verbose Logical; if TRUE (default), print progress messages.
#' @return A tibble with columns: \code{isoform_id}, \code{cds_exon_start}
#'   (integer), \code{cds_exon_end} (integer), \code{strand} (character),
#'   \code{cds_exon_rank} (integer, 1-based in 5'-to-3' order: ascending
#'   genomic start for "+", descending for "-").
#' @note Must use the same GTF input as \code{extractCdsAnnotations()} to
#'   ensure consistency between coding status and CDS exon data.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' cds_exons <- extractCdsExons(gtf_path, verbose = FALSE)
#' head(cds_exons)
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr filter arrange group_by mutate ungroup n
#' @importFrom rlang .data
extractCdsExons <- function(gtf, isoform_ids = NULL, verbose = TRUE) {

  gtf_tbl <- .loadGtfData(gtf, verbose)

  cds_rows <- dplyr::filter(gtf_tbl, .data$type == "CDS")

  if (nrow(cds_rows) == 0L) {
    if (verbose) message("  No CDS features found in GTF")
    return(tibble::tibble(
      isoform_id = character(0),
      cds_exon_start = integer(0),
      cds_exon_end = integer(0),
      strand = character(0),
      cds_exon_rank = integer(0)
    ))
  }

  result <- tibble::tibble(
    isoform_id = cds_rows$transcript_id,
    cds_exon_start = as.integer(cds_rows$start),
    cds_exon_end = as.integer(cds_rows$end),
    strand = as.character(cds_rows$strand)
  )

  if (!is.null(isoform_ids)) {
    result <- result[result$isoform_id %in% isoform_ids, ]
  }

  if (nrow(result) == 0L) {
    result$cds_exon_rank <- integer(0)
    return(result)
  }

  # Rank exons 5'→3': ascending start for +, descending start for -
  result <- result |>
    dplyr::arrange(.data$isoform_id, .data$cds_exon_start) |>
    dplyr::group_by(.data$isoform_id) |>
    dplyr::mutate(
      cds_exon_rank = if (.data$strand[1L] == "-") {
        rev(seq_len(dplyr::n()))
      } else {
        seq_len(dplyr::n())
      }
    ) |>
    dplyr::ungroup()

  if (verbose) {
    n_iso <- length(unique(result$isoform_id))
    message(sprintf("  %d CDS exon intervals for %d isoforms",
                    nrow(result), n_iso))
  }

  result
}


#' Annotate Region Types for Isoform Exons
#'
#' Classifies each exon segment of an isoform into a genomic region type
#' based on CDS boundaries: 5'UTR, CDS, 3'UTR, contains_orf_start,
#' contains_orf_stop, contains_orf_start_stop, or non_coding.
#'
#' Classification uses the isoform's actual exon coordinates (not union exon
#' segments). The classification is strand-aware: for minus-strand genes,
#' the 5'UTR is genomically to the right of the CDS and the 3'UTR is
#' genomically to the left. The \code{strand} column in \code{cds_metadata}
#' is required for correct directional labeling.
#'
#' @param isoform_union_mapping A tibble with at least columns
#'   \code{isoform_id}, \code{isoform_exon_start}, \code{isoform_exon_end}.
#'   Typically from \code{buildUnionExons()$isoform_union_mapping}.
#' @param cds_metadata A tibble from \code{extractCdsAnnotations()} with
#'   columns \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{strand}.
#' @return The input tibble with an added \code{region_type} column.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' table(annotated$region_type)
#' @export
#' @importFrom dplyr left_join mutate case_when
#' @importFrom rlang .data
annotateRegionTypes <- function(isoform_union_mapping, cds_metadata) {

  if (!all(c("isoform_id", "isoform_exon_start", "isoform_exon_end") %in%
           names(isoform_union_mapping))) {
    stop("isoform_union_mapping must have columns: isoform_id, ",
         "isoform_exon_start, isoform_exon_end")
  }

  required_cds_cols <- c("isoform_id", "coding_status", "cds_start",
                         "cds_stop", "strand")
  if (!all(required_cds_cols %in% names(cds_metadata))) {
    stop("cds_metadata must have columns: isoform_id, coding_status, ",
         "cds_start, cds_stop, strand")
  }

  # Join CDS info
  annotated <- dplyr::left_join(
    isoform_union_mapping,
    cds_metadata[, c("isoform_id", "coding_status", "cds_start", "cds_stop", "strand")],
    by = "isoform_id"
  )

  # Classify region types (strand-aware: minus-strand swaps 5'/3' and start/stop)
  annotated <- dplyr::mutate(
    annotated,
    region_type = dplyr::case_when(
      is.na(.data$coding_status) | .data$coding_status != "coding"
        ~ "non_coding",
      .data$isoform_exon_end < .data$cds_start
        ~ dplyr::if_else(.data$strand == "-", "3'UTR", "5'UTR"),
      .data$isoform_exon_start > .data$cds_stop
        ~ dplyr::if_else(.data$strand == "-", "5'UTR", "3'UTR"),
      .data$isoform_exon_start <= .data$cds_start &
        .data$isoform_exon_end >= .data$cds_stop
        ~ "contains_orf_start_stop",
      .data$isoform_exon_start <= .data$cds_start &
        .data$isoform_exon_end >= .data$cds_start
        ~ dplyr::if_else(.data$strand == "-", "contains_orf_stop",
                         "contains_orf_start"),
      .data$isoform_exon_start <= .data$cds_stop &
        .data$isoform_exon_end >= .data$cds_stop
        ~ dplyr::if_else(.data$strand == "-", "contains_orf_start",
                         "contains_orf_stop"),
      .data$isoform_exon_start >= .data$cds_start &
        .data$isoform_exon_end <= .data$cds_stop
        ~ "CDS",
      TRUE ~ "unknown"
    )
  )

  annotated
}


#' Load GTF data into a tibble
#'
#' Internal helper that normalizes a GTF path or GRanges object into a
#' tibble with standardized column names.
#'
#' @param gtf Character path or GRanges object.
#' @param verbose Logical.
#' @return A tibble with columns: seqnames, start, end, strand, type,
#'   transcript_id, gene_id.
#' @keywords internal
#' @importFrom tibble as_tibble
.loadGtfData <- function(gtf, verbose = TRUE) {
  if (is.character(gtf)) {
    if (!file.exists(gtf)) {
      stop(sprintf("GTF file not found: %s", gtf))
    }
    if (verbose) message("Loading GTF: ", gtf)
    gtf_tbl <- .readGtf(gtf)
  } else if (inherits(gtf, "GRanges")) {
    if (verbose) message("Using provided GRanges object")
    gtf_tbl <- tibble::as_tibble(gtf)
  } else {
    stop("gtf must be a file path (character) or GRanges object")
  }

  # Ensure required columns exist
  required <- c("type", "transcript_id")
  missing <- setdiff(required, names(gtf_tbl))
  if (length(missing) > 0L) {
    stop(sprintf("GTF data missing required columns: %s",
                 paste(missing, collapse = ", ")))
  }

  gtf_tbl
}
