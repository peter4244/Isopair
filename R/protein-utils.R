#' Parse hmmscan Domain Table Output
#'
#' Reads an hmmscan \code{--domtblout} output file and returns a tidy tibble
#' in the \code{domain_annotations} format expected by
#' \code{\link{classifyDomainEffects}()}.
#'
#' @param domtblout_path Character; path to an hmmscan \code{--domtblout}
#'   output file.
#' @param evalue_threshold Numeric; domains with independent E-value above
#'   this threshold are excluded (default \code{1e-5}).
#' @return A tibble with columns \code{isoform_id}, \code{domain_id},
#'   \code{domain_name}, \code{domain_start} (integer, envelope start),
#'   \code{domain_end} (integer, envelope end), \code{domain_score} (numeric),
#'   \code{domain_evalue} (numeric, independent E-value),
#'   \code{description} (character). One row per domain hit passing the
#'   E-value filter.
#' @examples
#' # Create a small mock domtblout file
#' lines <- c(
#'   "#                                                                  --- full",
#'   "#  target   accession  tlen  query  accession  qlen  E-value  score  bias",
#'   paste("Tubulin  PF00091.32  200  ENST001  -  450  1e-50  180.5  0.1",
#'         "1  2  1e-30  1e-28  170.2  0.0  5  195  10  210  8  212  0.95",
#'         "Tubulin family"),
#'   paste("Kinase   PF00069.10  300  ENST001  -  450  1e-40  150.3  0.2",
#'         "1  1  1e-20  1e-18  145.0  0.1  1  280  50  340  48  342  0.90",
#'         "Protein kinase domain")
#' )
#' tmp <- tempfile(fileext = ".domtblout")
#' writeLines(lines, tmp)
#' result <- parseHmmscanOutput(tmp)
#' unlink(tmp)
#' @export
#' @importFrom tibble tibble
parseHmmscanOutput <- function(domtblout_path, evalue_threshold = 1e-5) {

  # --- Input validation ---
  if (!is.character(domtblout_path) || length(domtblout_path) != 1L) {
    stop("'domtblout_path' must be a single character string.", call. = FALSE)
  }
  if (!file.exists(domtblout_path)) {
    stop(sprintf("File not found: %s", domtblout_path), call. = FALSE)
  }
  if (!is.numeric(evalue_threshold) || length(evalue_threshold) != 1L ||
      is.na(evalue_threshold) || evalue_threshold <= 0) {
    stop("'evalue_threshold' must be a single positive number.", call. = FALSE)
  }

  # --- Read and parse ---
  raw_lines <- readLines(domtblout_path, warn = FALSE)

  # Remove comment lines
  data_lines <- raw_lines[!grepl("^#", raw_lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]

  # Empty output template
  empty_result <- tibble::tibble(
    isoform_id = character(0),
    domain_id = character(0),
    domain_name = character(0),
    domain_start = integer(0),
    domain_end = integer(0),
    domain_score = numeric(0),
    domain_evalue = numeric(0),
    description = character(0)
  )

  if (length(data_lines) == 0L) return(empty_result)

  # Parse each line: first 22 fields are space-delimited, rest is description
  results <- vector("list", length(data_lines))
  n_kept <- 0L

  for (i in seq_along(data_lines)) {
    # Split into at most 23 fields (22 fixed + description as remainder)
    fields <- strsplit(trimws(data_lines[i]), "\\s+")[[1L]]

    if (length(fields) < 22L) next

    domain_name <- fields[1L]
    domain_accession <- fields[2L]
    query_name <- fields[4L]
    i_evalue <- as.numeric(fields[13L])
    domain_score <- as.numeric(fields[14L])
    env_from <- as.integer(fields[20L])
    env_to <- as.integer(fields[21L])

    # Description: everything from field 23 onward
    desc_text <- if (length(fields) >= 23L) {
      paste(fields[23L:length(fields)], collapse = " ")
    } else {
      NA_character_
    }

    # Filter by E-value
    if (is.na(i_evalue) || i_evalue > evalue_threshold) next

    # Strip version from Pfam accession (e.g., "PF00091.32" -> "PF00091")
    domain_id <- sub("\\..*$", "", domain_accession)

    n_kept <- n_kept + 1L
    results[[n_kept]] <- tibble::tibble(
      isoform_id = query_name,
      domain_id = domain_id,
      domain_name = domain_name,
      domain_start = env_from,
      domain_end = env_to,
      domain_score = domain_score,
      domain_evalue = i_evalue,
      description = desc_text
    )
  }

  if (n_kept == 0L) return(empty_result)
  dplyr::bind_rows(results[seq_len(n_kept)])
}


#' Run hmmscan and Parse Results
#'
#' Wrapper that runs HMMER's \code{hmmscan} on protein sequences and returns
#' parsed domain annotations. Requires \code{hmmscan} to be installed and
#' available on the system PATH.
#'
#' @param sequences Either a named character vector (names = isoform IDs,
#'   values = amino acid sequences) or a character path to an existing
#'   FASTA file.
#' @param hmm_database Character; path to an \code{hmmpress}'d HMM database
#'   (e.g., \code{Pfam-A.hmm}).
#' @param n_cores Integer; number of CPU cores for hmmscan (default 4).
#' @param evalue_threshold Numeric; E-value threshold for both full-sequence
#'   and per-domain filtering (default \code{1e-5}).
#' @param temp_dir Character; directory for temporary files
#'   (default \code{tempdir()}).
#' @return A tibble in the same format as \code{\link{parseHmmscanOutput}()}.
#' @examples
#' \dontrun{
#' seqs <- c(iso1 = "MKTLLVVFLAG...", iso2 = "MKTLLVVFLAG...")
#' result <- runHmmscan(seqs, "/path/to/Pfam-A.hmm")
#' }
#' @export
#' @importFrom tibble tibble
runHmmscan <- function(sequences, hmm_database, n_cores = 4L,
                       evalue_threshold = 1e-5, temp_dir = tempdir()) {

  # --- Input validation ---
  if (!is.character(sequences) || length(sequences) == 0L) {
    stop("'sequences' must be a non-empty character vector or a path to a ",
         "FASTA file.", call. = FALSE)
  }
  if (!is.character(hmm_database) || length(hmm_database) != 1L) {
    stop("'hmm_database' must be a single character string.", call. = FALSE)
  }
  if (!file.exists(hmm_database)) {
    stop(sprintf("HMM database not found: %s", hmm_database), call. = FALSE)
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    stop("'n_cores' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(evalue_threshold) || length(evalue_threshold) != 1L ||
      is.na(evalue_threshold) || evalue_threshold <= 0) {
    stop("'evalue_threshold' must be a single positive number.", call. = FALSE)
  }

  # Check hmmscan availability
  hmmscan_path <- Sys.which("hmmscan")
  if (!nzchar(hmmscan_path)) {
    stop("'hmmscan' is not installed or not found on the system PATH. ",
         "Install HMMER (http://hmmer.org/) and ensure 'hmmscan' is ",
         "accessible.", call. = FALSE)
  }

  # Determine input: file path or named vector
  if (length(sequences) == 1L && is.null(names(sequences)) &&
      file.exists(sequences)) {
    input_fasta <- sequences
    cleanup_fasta <- FALSE
  } else if (length(sequences) == 1L && !is.null(names(sequences)) &&
             file.exists(sequences)) {
    # Ambiguous: named single-element vector where value is also a file path.
    # Treat as a sequence vector (names present = sequences).
    input_fasta <- .writeSequencesToFasta(sequences, temp_dir)
    cleanup_fasta <- TRUE
  } else {
    # Named character vector of sequences
    if (is.null(names(sequences))) {
      stop("'sequences' must be a named character vector (names = isoform ",
           "IDs) or a path to an existing FASTA file.", call. = FALSE)
    }
    input_fasta <- .writeSequencesToFasta(sequences, temp_dir)
    cleanup_fasta <- TRUE
  }

  # Prepare output path
  output_domtblout <- file.path(temp_dir,
                                 paste0("hmmscan_domtblout_",
                                        format(Sys.time(), "%Y%m%d%H%M%S"),
                                        "_", Sys.getpid(), ".txt"))

  # Build command
  cmd <- sprintf(
    "%s --cpu %d --domtblout %s --noali -E %g --domE %g %s %s",
    shQuote(hmmscan_path), n_cores, shQuote(output_domtblout),
    evalue_threshold, evalue_threshold,
    shQuote(hmm_database), shQuote(input_fasta)
  )

  # Run hmmscan
  exit_code <- tryCatch({
    system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)
  }, error = function(e) {
    stop(sprintf("hmmscan execution failed: %s", conditionMessage(e)),
         call. = FALSE)
  })

  # Parse results
  result <- tryCatch({
    if (file.exists(output_domtblout)) {
      parseHmmscanOutput(output_domtblout, evalue_threshold = evalue_threshold)
    } else {
      stop("hmmscan did not produce output file.", call. = FALSE)
    }
  }, finally = {
    # Cleanup
    if (cleanup_fasta && file.exists(input_fasta)) unlink(input_fasta)
    if (file.exists(output_domtblout)) unlink(output_domtblout)
  })

  result
}


#' In Silico Protein Digestion
#'
#' Performs in silico enzymatic digestion of protein sequences and returns
#' all peptides within specified length bounds.
#'
#' @param protein_sequences Named character vector; names are isoform IDs,
#'   values are amino acid sequences. Terminal stop codon markers
#'   (\code{"*"}) are removed before digestion.
#' @param enzyme Character; enzyme to use for digestion (default
#'   \code{"trypsin"}). Currently supported: \code{"trypsin"} (cleaves after
#'   K or R, except when followed by P).
#' @param min_length Integer; minimum peptide length to retain (default 7).
#' @param max_length Integer; maximum peptide length to retain (default 30).
#' @return A tibble with columns \code{isoform_id} (character),
#'   \code{peptide} (character), \code{start_aa} (integer, 1-based start
#'   position in the protein), \code{end_aa} (integer, 1-based end position).
#'   One row per peptide per isoform, filtered to peptides within the
#'   specified length bounds.
#' @examples
#' seqs <- c(iso1 = "MKTLLVKPFLAG*", iso2 = "MKRAATYR")
#' result <- digestProtein(seqs)
#' result
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
digestProtein <- function(protein_sequences, enzyme = "trypsin",
                          min_length = 7L, max_length = 30L) {

  # --- Input validation ---
  if (!is.character(protein_sequences) || length(protein_sequences) == 0L) {
    stop("'protein_sequences' must be a non-empty named character vector.",
         call. = FALSE)
  }
  if (is.null(names(protein_sequences))) {
    stop("'protein_sequences' must be named (names = isoform IDs).",
         call. = FALSE)
  }

  enzyme <- match.arg(enzyme, choices = "trypsin")

  min_length <- as.integer(min_length)
  max_length <- as.integer(max_length)
  if (is.na(min_length) || min_length < 1L) {
    stop("'min_length' must be a positive integer.", call. = FALSE)
  }
  if (is.na(max_length) || max_length < min_length) {
    stop("'max_length' must be an integer >= 'min_length'.", call. = FALSE)
  }

  # Empty output template
  empty_result <- tibble::tibble(
    isoform_id = character(0),
    peptide = character(0),
    start_aa = integer(0),
    end_aa = integer(0)
  )

  all_results <- vector("list", length(protein_sequences))

  for (idx in seq_along(protein_sequences)) {
    iso_id <- names(protein_sequences)[idx]
    seq <- protein_sequences[idx]

    # Remove terminal stop codon marker(s) and drop name
    seq <- unname(sub("\\*+$", "", seq))

    if (!nzchar(seq)) next

    # Trypsin digestion: cleave after K or R, unless followed by P
    # Find cleavage positions
    peptides <- .digestTrypsin(seq)

    if (length(peptides) == 0L) next

    # Compute start/end positions
    starts <- integer(length(peptides))
    ends <- integer(length(peptides))
    pos <- 1L
    for (j in seq_along(peptides)) {
      pep_len <- nchar(peptides[j])
      starts[j] <- pos
      ends[j] <- pos + pep_len - 1L
      pos <- pos + pep_len
    }

    # Filter by length
    pep_lengths <- nchar(peptides)
    keep <- pep_lengths >= min_length & pep_lengths <= max_length

    if (!any(keep)) next

    all_results[[idx]] <- tibble::tibble(
      isoform_id = iso_id,
      peptide = peptides[keep],
      start_aa = starts[keep],
      end_aa = ends[keep]
    )
  }

  result <- dplyr::bind_rows(all_results)
  if (nrow(result) == 0L) return(empty_result)
  result
}


#' Query Peptides Against a Local PeptideAtlas List
#'
#' Performs hash-based lookup of peptide sequences against a local
#' PeptideAtlas peptide list file. This provides fast membership testing
#' for large peptide sets.
#'
#' @param peptides Character vector of peptide sequences to query.
#' @param peptide_atlas_path Character; path to a text file with one peptide
#'   sequence per line (no headers).
#' @return A tibble with columns \code{peptide} (character) and
#'   \code{found} (logical).
#' @examples
#' # Create a small mock PeptideAtlas file
#' atlas <- tempfile(fileext = ".txt")
#' writeLines(c("MKTLLVK", "PFLAGR", "AATYR"), atlas)
#' result <- queryPeptideAtlas(c("MKTLLVK", "UNKNOWN", "AATYR"), atlas)
#' unlink(atlas)
#' result
#' @export
#' @importFrom tibble tibble
queryPeptideAtlas <- function(peptides, peptide_atlas_path) {

  # --- Input validation ---
  if (!is.character(peptides) || length(peptides) == 0L) {
    stop("'peptides' must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(peptide_atlas_path) || length(peptide_atlas_path) != 1L) {
    stop("'peptide_atlas_path' must be a single character string.",
         call. = FALSE)
  }
  if (!file.exists(peptide_atlas_path)) {
    stop(sprintf("File not found: %s", peptide_atlas_path), call. = FALSE)
  }

  # Load peptides into a hash-based lookup (environment)
  atlas_lines <- readLines(peptide_atlas_path, warn = FALSE)
  atlas_lines <- trimws(atlas_lines)
  atlas_lines <- atlas_lines[nzchar(atlas_lines)]

  atlas_env <- new.env(hash = TRUE, parent = emptyenv(),
                       size = length(atlas_lines))
  for (pep in atlas_lines) {
    atlas_env[[pep]] <- TRUE
  }

  # Query
  found <- vapply(peptides, function(p) exists(p, envir = atlas_env,
                                                inherits = FALSE),
                  logical(1), USE.NAMES = FALSE)

  tibble::tibble(
    peptide = peptides,
    found = found
  )
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Write Protein Sequences to a Temporary FASTA File
#'
#' @param sequences Named character vector of AA sequences.
#' @param temp_dir Directory for the temporary file.
#' @return Path to the written FASTA file.
#' @keywords internal
.writeSequencesToFasta <- function(sequences, temp_dir) {
  fasta_path <- file.path(temp_dir,
                           paste0("isopair_seqs_",
                                  format(Sys.time(), "%Y%m%d%H%M%S"),
                                  "_", Sys.getpid(), ".fasta"))
  con <- file(fasta_path, open = "w")
  on.exit(close(con))

  for (i in seq_along(sequences)) {
    writeLines(paste0(">", names(sequences)[i]), con)
    writeLines(sequences[i], con)
  }

  fasta_path
}


#' Trypsin In Silico Digestion
#'
#' Cleaves after K or R, except when followed by P.
#'
#' @param seq Character; single amino acid sequence.
#' @return Character vector of peptide fragments.
#' @keywords internal
.digestTrypsin <- function(seq) {
  chars <- strsplit(seq, "")[[1L]]
  n <- length(chars)

  if (n == 0L) return(character(0))

  # Find cleavage sites: positions where char is K or R and next char is not P

  cleavage_sites <- integer(0)
  for (i in seq_len(n - 1L)) {
    if (chars[i] %in% c("K", "R") && chars[i + 1L] != "P") {
      cleavage_sites <- c(cleavage_sites, i)
    }
  }

  # Also check the last character (K or R at end of sequence = cleavage)
  # No need; the last fragment includes everything after the last cleavage

  if (length(cleavage_sites) == 0L) {
    return(seq)
  }

  # Build peptides from cleavage sites
  peptides <- character(length(cleavage_sites) + 1L)
  prev_end <- 0L
  for (k in seq_along(cleavage_sites)) {
    peptides[k] <- substr(seq, prev_end + 1L, cleavage_sites[k])
    prev_end <- cleavage_sites[k]
  }
  # Last peptide
  peptides[length(cleavage_sites) + 1L] <- substr(seq, prev_end + 1L, n)

  # Remove empty fragments
  peptides[nzchar(peptides)]
}
