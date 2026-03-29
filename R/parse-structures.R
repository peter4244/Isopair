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
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr filter select group_by summarise mutate rename n first
#'   arrange n_distinct bind_rows
#' @importFrom rlang .data
parseIsoformStructures <- function(gtf_path, isoform_ids = NULL,
                                   verbose = TRUE) {
  if (!file.exists(gtf_path)) {
    stop(sprintf("GTF file not found: %s", gtf_path))
  }

  if (verbose) message("Loading GTF: ", gtf_path)
  gtf_tbl <- .readGtf(gtf_path)

  if (verbose) message("Filtering to exon features...")
  exons <- gtf_tbl |>
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


#' Deduplicate Merged Isoform Structures
#'
#' Merges two structures tibbles (e.g., reference and novel annotations),
#' detects isoforms with identical exon structures within each gene, and
#' keeps the preferred source's ID. Useful when combining GENCODE with
#' PacBio or other de novo assemblies.
#'
#' @param structures_a A tibble from [parseIsoformStructures()].
#' @param structures_b A tibble from [parseIsoformStructures()].
#' @param prefer Which source to keep when duplicates are found:
#'   \code{"a"} (default) or \code{"b"}.
#' @param strip_gene_version Logical; if TRUE, strip \code{.N} version
#'   suffix from gene_id before grouping. Useful when one GTF has
#'   versioned gene IDs and the other does not.
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A named list with:
#'   \describe{
#'     \item{structures}{Merged, deduplicated tibble with a \code{source}
#'       column (\code{"a"} or \code{"b"}).}
#'     \item{id_mapping}{Tibble with \code{kept_id}, \code{removed_id},
#'       \code{gene_id} for each duplicate pair.}
#'     \item{n_total}{Total unique isoforms retained.}
#'     \item{n_duplicates}{Count of duplicates removed.}
#'   }
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' s1 <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' s2 <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' dedup <- deduplicateStructures(s1, s2, verbose = FALSE)
#' dedup$n_duplicates  # all flagged as duplicates
#' @export
deduplicateStructures <- function(structures_a, structures_b,
                                  prefer = c("a", "b"),
                                  strip_gene_version = FALSE,
                                  verbose = TRUE) {
  prefer <- match.arg(prefer)

  # Tag sources
  structures_a$source <- "a"
  structures_b$source <- "b"
  merged <- dplyr::bind_rows(structures_a, structures_b)

  if (nrow(merged) == 0L) {
    return(list(
      structures = merged,
      id_mapping = tibble::tibble(
        kept_id = character(0), removed_id = character(0),
        gene_id = character(0)
      ),
      n_total = 0L, n_duplicates = 0L
    ))
  }

  # Compute fingerprint for each isoform
  merged$fingerprint <- vapply(seq_len(nrow(merged)), function(i) {
    paste(merged$chr[i], merged$strand[i],
          paste(merged$exon_starts[[i]], collapse = ","),
          paste(merged$exon_ends[[i]], collapse = ","),
          sep = "|")
  }, character(1))

  # Grouping key: gene_id (optionally strip version)
  if (strip_gene_version) {
    merged$group_key <- sub("\\.[0-9]+$", "", merged$gene_id)
  } else {
    merged$group_key <- merged$gene_id
  }

  # Check for same isoform_id in both sources with different structures
  id_counts <- table(merged$isoform_id)
  dup_ids <- names(id_counts[id_counts > 1L])
  if (length(dup_ids) > 0L) {
    for (did in dup_ids) {
      rows <- merged[merged$isoform_id == did, ]
      fps <- unique(rows$fingerprint)
      if (length(fps) > 1L) {
        warning(sprintf(
          "Isoform '%s' appears in both sources with different exon structures. Keeping '%s' source.",
          did, prefer
        ))
      }
    }
  }

  # Find duplicates within each gene group by fingerprint
  id_mapping <- list()
  to_remove <- character(0)
  keep_source <- prefer
  drop_source <- if (prefer == "a") "b" else "a"

  gene_groups <- split(seq_len(nrow(merged)), merged$group_key)

  for (gene in names(gene_groups)) {
    idxs <- gene_groups[[gene]]
    gene_rows <- merged[idxs, ]
    fp_groups <- split(seq_len(nrow(gene_rows)), gene_rows$fingerprint)

    for (fp_idxs in fp_groups) {
      fp_rows <- gene_rows[fp_idxs, ]
      if (nrow(fp_rows) < 2L) next

      sources_present <- unique(fp_rows$source)
      if (length(sources_present) < 2L) {
        # Same source duplicates — keep first, remove rest
        for (k in seq_len(nrow(fp_rows))[-1L]) {
          to_remove <- c(to_remove, paste(fp_rows$isoform_id[k],
                                           fp_rows$source[k], sep = "::"))
          id_mapping[[length(id_mapping) + 1L]] <- tibble::tibble(
            kept_id = fp_rows$isoform_id[1L],
            removed_id = fp_rows$isoform_id[k],
            gene_id = gene
          )
        }
        next
      }

      # Cross-source duplicates: keep preferred, remove others
      keep_rows <- fp_rows[fp_rows$source == keep_source, ]
      drop_rows <- fp_rows[fp_rows$source == drop_source, ]
      if (nrow(keep_rows) == 0L) next

      kept_id <- keep_rows$isoform_id[1L]
      for (k in seq_len(nrow(drop_rows))) {
        to_remove <- c(to_remove, paste(drop_rows$isoform_id[k],
                                         drop_rows$source[k], sep = "::"))
        id_mapping[[length(id_mapping) + 1L]] <- tibble::tibble(
          kept_id = kept_id,
          removed_id = drop_rows$isoform_id[k],
          gene_id = gene
        )
      }
    }
  }

  # Build removal index using isoform_id::source composite key
  merged$composite_key <- paste(merged$isoform_id, merged$source, sep = "::")
  result <- merged[!merged$composite_key %in% to_remove, ]

  # Clean up temporary columns
  result$fingerprint <- NULL
  result$group_key <- NULL
  result$composite_key <- NULL

  id_mapping_df <- if (length(id_mapping) > 0L) {
    dplyr::bind_rows(id_mapping)
  } else {
    tibble::tibble(kept_id = character(0), removed_id = character(0),
                   gene_id = character(0))
  }

  n_duplicates <- nrow(id_mapping_df)
  if (verbose) {
    message(sprintf(
      "Deduplication: %d total isoforms, %d duplicates removed, %d retained",
      nrow(merged), n_duplicates, nrow(result)
    ))
  }

  list(
    structures = tibble::as_tibble(result),
    id_mapping = id_mapping_df,
    n_total = nrow(result),
    n_duplicates = n_duplicates
  )
}
