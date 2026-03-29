#' Read a GTF File into a Tibble
#'
#' Reads a GTF annotation file and returns a tibble with standardized columns.
#' When rtracklayer is installed, uses its robust C-level parser. Otherwise
#' falls back to a lightweight base-R parser that requires no Bioconductor
#' dependencies. Supports plain-text and gzipped GTF files.
#'
#' @param path Path to a GTF file (can be gzipped).
#' @return A tibble with columns: seqnames, start, end, strand, type,
#'   transcript_id, gene_id.
#' @keywords internal
#' @importFrom tibble as_tibble tibble
.readGtf <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("GTF file not found: %s", path))
  }

  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gr <- rtracklayer::import(path)
    return(tibble::as_tibble(gr))
  }

  .readGtfBase(path)
}

#' Lightweight Base-R GTF Parser
#'
#' Parses a GTF file using only base R, extracting the columns needed by
#' Isopair. Used as a fallback when rtracklayer is not installed.
#'
#' @param path Path to a GTF file (can be gzipped).
#' @return A tibble with columns: seqnames, start, end, strand, type,
#'   transcript_id, gene_id.
#' @keywords internal
#' @importFrom tibble tibble
.readGtfBase <- function(path) {
  # Read lines, skipping comments
  con <- if (grepl("\\.gz$", path)) gzfile(path) else file(path)
  on.exit(close(con))
  lines <- readLines(con)
  lines <- lines[!grepl("^#", lines) & nchar(lines) > 0L]

  if (length(lines) == 0L) {
    return(tibble::tibble(
      seqnames = character(0), start = integer(0), end = integer(0),
      strand = character(0), type = character(0),
      transcript_id = character(0), gene_id = character(0)
    ))
  }

  # Parse the 9 tab-delimited GTF columns
  fields <- strsplit(lines, "\t", fixed = TRUE)
  n_fields <- vapply(fields, length, integer(1))
  if (any(n_fields < 9L)) {
    warning(sprintf("Skipping %d lines with fewer than 9 fields",
                    sum(n_fields < 9L)))
    fields <- fields[n_fields >= 9L]
  }

  seqnames <- vapply(fields, `[`, character(1), 1L)
  type     <- vapply(fields, `[`, character(1), 3L)
  start    <- as.integer(vapply(fields, `[`, character(1), 4L))
  end      <- as.integer(vapply(fields, `[`, character(1), 5L))
  strand   <- vapply(fields, `[`, character(1), 7L)
  attrs    <- vapply(fields, `[`, character(1), 9L)

  # Extract transcript_id and gene_id from the attributes column
  transcript_id <- sub(
    '.*transcript_id\\s+"([^"]+)".*', "\\1", attrs, perl = TRUE
  )
  transcript_id[!grepl("transcript_id", attrs)] <- NA_character_

  gene_id <- sub(
    '.*gene_id\\s+"([^"]+)".*', "\\1", attrs, perl = TRUE
  )
  gene_id[!grepl("gene_id", attrs)] <- NA_character_

  tibble::tibble(
    seqnames = seqnames,
    start = start,
    end = end,
    strand = strand,
    type = type,
    transcript_id = transcript_id,
    gene_id = gene_id
  )
}
