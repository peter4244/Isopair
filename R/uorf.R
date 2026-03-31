# ==============================================================================
# uorf.R — 5'UTR feature scanning, uORF detection, comparison, and attribution
# ==============================================================================

# --- Internal helpers ---------------------------------------------------------

#' Extract 5'UTR Genomic Coordinate Intervals
#'
#' Returns the genomic intervals comprising the 5'UTR for a single coding
#' isoform, based on its exon structure and CDS boundaries.
#'
#' @param strand Character; "+" or "-".
#' @param exon_starts Integer vector of exon start coordinates (1-based).
#' @param exon_ends Integer vector of exon end coordinates (1-based, closed).
#' @param cds_start Integer; min genomic CDS coordinate.
#' @param cds_stop Integer; max genomic CDS coordinate.
#' @return A data.frame with columns \code{start} and \code{end} (integer) of
#'   genomic intervals comprising the 5'UTR, in genomic order. Empty data.frame
#'   if no 5'UTR exists.
#' @keywords internal
.extract5UtrCoords <- function(strand, exon_starts, exon_ends,
                               cds_start, cds_stop) {
  starts <- integer(0)
  ends <- integer(0)

  if (strand == "+") {
    # + strand: translation starts at cds_start. 5'UTR = exonic bases < cds_start
    tstart <- cds_start
    for (i in seq_along(exon_starts)) {
      es <- exon_starts[i]
      ee <- exon_ends[i]
      if (ee < tstart) {
        # Entire exon is 5'UTR
        starts <- c(starts, es)
        ends <- c(ends, ee)
      } else if (es < tstart) {
        # Exon spans translation start — partial 5'UTR
        starts <- c(starts, es)
        ends <- c(ends, tstart - 1L)
      }
      # Once we've reached or passed tstart, no more 5'UTR
      if (es >= tstart) break
    }
  } else {
    # - strand: translation starts at cds_stop (highest genomic = 5' end)
    # Walk exons from highest to lowest (5' → 3' in transcript)
    tstart <- cds_stop
    for (i in rev(seq_along(exon_starts))) {
      es <- exon_starts[i]
      ee <- exon_ends[i]
      if (es > tstart) {
        # Entire exon is 5'UTR
        starts <- c(starts, es)
        ends <- c(ends, ee)
      } else if (ee > tstart) {
        # Exon spans translation start — partial 5'UTR
        starts <- c(starts, tstart + 1L)
        ends <- c(ends, ee)
      }
      if (ee <= tstart) break
    }
    # Sort intervals by genomic position (ascending)
    if (length(starts) > 1L) {
      ord <- order(starts)
      starts <- starts[ord]
      ends <- ends[ord]
    }
  }

  data.frame(start = starts, end = ends)
}


#' Map Transcript-Space Position to Genomic Coordinate
#'
#' Given a 1-based position in transcript space (where position 1 is the first
#' base of the first exon in 5'-to-3' order), returns the corresponding genomic
#' coordinate.
#'
#' @param transcript_pos Integer; 1-based position in transcript space.
#' @param strand Character; "+" or "-".
#' @param exon_starts Integer vector of exon start coordinates (1-based).
#' @param exon_ends Integer vector of exon end coordinates (1-based, closed).
#' @return Integer genomic position, or \code{NA_integer_} if the position
#'   exceeds the transcript length.
#' @keywords internal
.transcriptToGenomic <- function(transcript_pos, strand, exon_starts, exon_ends) {
  if (is.na(transcript_pos) || transcript_pos < 1L) return(NA_integer_)

  if (strand == "+") {
    # Walk exons in ascending order (5' → 3')
    ord <- order(exon_starts)
  } else {
    # Walk exons in descending order (5' → 3' for minus strand)
    ord <- order(exon_starts, decreasing = TRUE)
  }

  remaining <- as.integer(transcript_pos)
  for (idx in ord) {
    exon_len <- exon_ends[idx] - exon_starts[idx] + 1L
    if (remaining <= exon_len) {
      if (strand == "+") {
        return(exon_starts[idx] + remaining - 1L)
      } else {
        return(exon_ends[idx] - remaining + 1L)
      }
    }
    remaining <- remaining - exon_len
  }

  NA_integer_
}


#' Compute 5'UTR Length in Transcript Coordinates
#'
#' @param strand Character; "+" or "-".
#' @param exon_starts Integer vector of exon start coordinates.
#' @param exon_ends Integer vector of exon end coordinates.
#' @param cds_start Integer; min genomic CDS coordinate.
#' @param cds_stop Integer; max genomic CDS coordinate.
#' @return Integer 5'UTR length in nucleotides.
#' @keywords internal
.compute5UtrLength <- function(strand, exon_starts, exon_ends,
                               cds_start, cds_stop) {
  coords <- .extract5UtrCoords(strand, exon_starts, exon_ends,
                               cds_start, cds_stop)
  if (nrow(coords) == 0L) return(0L)
  sum(coords$end - coords$start + 1L)
}


#' Compute Total Transcript Length
#'
#' @param exon_starts Integer vector of exon start coordinates.
#' @param exon_ends Integer vector of exon end coordinates.
#' @return Integer total transcript length in nucleotides.
#' @keywords internal
.computeTranscriptLength <- function(exon_starts, exon_ends) {
  sum(exon_ends - exon_starts + 1L)
}


#' Score Kozak Context Around an ATG
#'
#' Simple Kozak scoring: checks -3 position (A/G = strong) and +4 position
#' (G = strong). Returns a score 0-2 and the context string.
#'
#' @param seq Character; full transcript sequence (uppercase).
#' @param atg_pos Integer; 1-based position of the 'A' of ATG in transcript.
#' @return A list with \code{score} (integer 0-2) and \code{context} (character
#'   string of positions -3 to +4, or NA if insufficient context).
#' @keywords internal
.scoreKozak <- function(seq, atg_pos) {
  seq_len <- nchar(seq)
  # Need position atg_pos-3 through atg_pos+5 (the G of ATG is at atg_pos+2,
  # +4 position is at atg_pos+3)
  if (atg_pos < 4L || (atg_pos + 3L) > seq_len) {
    return(list(score = NA_integer_, context = NA_character_))
  }

  context <- substr(seq, atg_pos - 3L, atg_pos + 3L)  # 7 nt: -3 to +4
  minus3 <- substr(seq, atg_pos - 3L, atg_pos - 3L)
  plus4 <- substr(seq, atg_pos + 3L, atg_pos + 3L)

  score <- 0L
  if (minus3 %in% c("A", "G")) score <- score + 1L
  if (plus4 == "G") score <- score + 1L

  list(score = score, context = context)
}


#' Default Kozak Position Weight Matrix
#'
#' Returns a log-odds PWM for positions -6 to -1 and +4 to +5 relative to the
#' ATG start codon. Weights are derived from nucleotide frequencies observed at
#' verified vertebrate translation initiation sites (Kozak 1987, Cavener & Ray
#' 1991) versus a uniform (25\%) background. Positive values indicate enrichment
#' (favorable for initiation); negative values indicate depletion.
#'
#' @return A numeric matrix with rows A, C, G, T and columns named by position
#'   (e.g. "-6", "-5", ..., "+4", "+5").
#' @references
#' Kozak M (1987) Nucleic Acids Res 15:8125-8148.
#' Cavener DR, Ray SC (1991) Nucleic Acids Res 19:3185-3192.
#' Noderer WL et al. (2014) Mol Syst Biol 10:748.
#' @keywords internal
.defaultKozakPWM <- function() {
  # Observed frequencies at vertebrate TIS sites (Cavener & Ray 1991,
  # 3012 vertebrate mRNAs, positions -6 to -1 and +4 to +5):
  #
  #         -6    -5    -4    -3    -2    -1    +4    +5
  #   A   0.23  0.25  0.25  0.37  0.27  0.19  0.23  0.17
  #   C   0.35  0.22  0.32  0.19  0.33  0.39  0.16  0.31
  #   G   0.23  0.33  0.25  0.34  0.17  0.23  0.46  0.34
  #   T   0.19  0.20  0.18  0.10  0.23  0.19  0.15  0.18
  #
  # Log-odds = log2(observed_freq / 0.25):
  freqs <- matrix(c(
    # -6     -5     -4     -3     -2     -1     +4     +5
    0.23,  0.25,  0.25,  0.37,  0.27,  0.19,  0.23,  0.17,  # A
    0.35,  0.22,  0.32,  0.19,  0.33,  0.39,  0.16,  0.31,  # C
    0.23,  0.33,  0.25,  0.34,  0.17,  0.23,  0.46,  0.34,  # G
    0.19,  0.20,  0.18,  0.10,  0.23,  0.19,  0.15,  0.18   # T
  ), nrow = 4, byrow = TRUE,
  dimnames = list(c("A", "C", "G", "T"),
                  c("-6", "-5", "-4", "-3", "-2", "-1", "+4", "+5")))

  log2(freqs / 0.25)
}


#' Score Kozak Context Using a Position Weight Matrix
#'
#' Computes a continuous Kozak initiation strength score for one or more ATG
#' sites using a log-odds position weight matrix (PWM). The default PWM is
#' derived from vertebrate translation initiation site frequencies (Cavener &
#' Ray 1991). Users may supply a custom PWM (e.g., from Noderer et al. 2014
#' FACS-seq measurements).
#'
#' The scored positions are -6 to -1 (upstream of ATG) and +4 to +5 (downstream
#' of ATG), where ATG occupies positions +1/+2/+3. The score is the sum of
#' log-odds weights across all scorable positions.
#'
#' @param sequences Character vector of transcript sequences (uppercase), or a
#'   Biostrings \code{DNAStringSet}. Recycled to match the length of
#'   \code{atg_positions} if length 1.
#' @param atg_positions Integer vector of 1-based positions of the 'A' of ATG
#'   in each sequence.
#' @param weights Numeric matrix (optional). Rows must be named A, C, G, T;
#'   columns named by position (e.g. "-3", "+4"). If \code{NULL} (default), uses
#'   the built-in vertebrate TIS PWM from \code{.defaultKozakPWM()}.
#' @return A data frame with columns:
#' \describe{
#'   \item{score}{Numeric; sum of log-odds weights across scored positions.
#'     Higher values indicate stronger initiation context.}
#'   \item{context}{Character; the extracted sequence context from position -6
#'     to +5 (12 nt), or \code{NA} if the window extends beyond the sequence.}
#'   \item{n_scored}{Integer; number of positions that contributed to the score
#'     (may be less than the full window at sequence edges).}
#' }
#' @examples
#' # Strong Kozak: A at -3, G at +4
#' scoreKozakPWM("GCCGCCACCATGGCG", atg_positions = 10L)
#'
#' # Score multiple ATGs in one sequence
#' seq <- "AAAATGCCCGCCACCATGGCGAAA"
#' scoreKozakPWM(seq, atg_positions = c(4L, 14L))
#'
#' # Use a custom weight matrix (must have row names A/C/G/T)
#' custom_pwm <- matrix(0, nrow = 4, ncol = 2,
#'                      dimnames = list(c("A","C","G","T"), c("-3", "+4")))
#' custom_pwm["A", "-3"] <- 1.0
#' custom_pwm["G", "-3"] <- 0.8
#' custom_pwm["G", "+4"] <- 1.0
#' scoreKozakPWM("GCCACCATGGCG", atg_positions = 7L, weights = custom_pwm)
#'
#' @references
#' Kozak M (1987) Nucleic Acids Res 15:8125-8148.
#' Cavener DR, Ray SC (1991) Nucleic Acids Res 19:3185-3192.
#' Noderer WL et al. (2014) Mol Syst Biol 10:748.
#' @export
scoreKozakPWM <- function(sequences, atg_positions, weights = NULL) {
  # Handle DNAStringSet input
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- as.character(sequences)
  }
  sequences <- toupper(sequences)

  # Recycle length-1 sequence to match atg_positions
  if (length(sequences) == 1L && length(atg_positions) > 1L) {
    sequences <- rep(sequences, length(atg_positions))
  }
  if (length(sequences) != length(atg_positions)) {
    stop("'sequences' and 'atg_positions' must have the same length ",
         "(or sequences length 1 for recycling).")
  }

  # Validate and set up weight matrix
  if (is.null(weights)) {
    weights <- .defaultKozakPWM()
  }
  if (!is.matrix(weights) || !all(c("A", "C", "G", "T") %in% rownames(weights))) {
    stop("'weights' must be a numeric matrix with row names A, C, G, T.")
  }

  # Parse scored positions from column names
  pos_names <- colnames(weights)
  if (is.null(pos_names)) {
    stop("'weights' must have column names indicating positions ",
         "(e.g. '-6', '-3', '+4').")
  }
  pos_offsets <- as.integer(gsub("\\+", "", pos_names))
  if (any(is.na(pos_offsets))) {
    stop("Column names of 'weights' must be parseable as integers ",
         "(e.g. '-6', '-3', '+4').")
  }
  # Offsets are relative to ATG position +1: position -3 means atg_pos - 3,
  # position +4 means atg_pos + 3 (since +1 = A, +2 = T, +3 = G, +4 = next)
  seq_offsets <- ifelse(pos_offsets < 0, pos_offsets, pos_offsets - 1L)

  n <- length(sequences)
  scores <- numeric(n)
  contexts <- character(n)
  n_scored <- integer(n)

  for (i in seq_len(n)) {
    seq_i <- sequences[i]
    atg_i <- atg_positions[i]
    seq_len_i <- nchar(seq_i)

    if (is.na(atg_i) || is.na(seq_i)) {
      scores[i] <- NA_real_
      contexts[i] <- NA_character_
      n_scored[i] <- NA_integer_
      next
    }

    # Full context window: -6 to +5 (12 nt)
    ctx_start <- atg_i - 6L
    ctx_end <- atg_i + 4L  # +5 position = atg_pos + 4
    if (ctx_start >= 1L && ctx_end <= seq_len_i) {
      contexts[i] <- substr(seq_i, ctx_start, ctx_end)
    } else {
      contexts[i] <- NA_character_
    }

    # Score each position in the weight matrix
    s <- 0
    ns <- 0L
    for (j in seq_along(pos_offsets)) {
      nt_pos <- atg_i + seq_offsets[j]
      if (nt_pos < 1L || nt_pos > seq_len_i) next
      nt <- substr(seq_i, nt_pos, nt_pos)
      if (nt %in% rownames(weights)) {
        s <- s + weights[nt, j]
        ns <- ns + 1L
      }
    }

    scores[i] <- if (ns > 0L) s else NA_real_
    n_scored[i] <- ns
  }

  data.frame(score = scores, context = contexts, n_scored = n_scored,
             stringsAsFactors = FALSE)
}


# --- Exported functions -------------------------------------------------------

#' Scan 5'UTR Sequence Features
#'
#' Extracts a comprehensive set of 5'UTR sequence features for coding isoforms,
#' including ATG/stop codon counts, frame-stratified counts, ORF detection,
#' density measures, and Kozak context scoring.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()} with
#'   columns \code{isoform_id}, \code{gene_id}, \code{strand},
#'   \code{exon_starts} (list), \code{exon_ends} (list).
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()} with
#'   columns \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{strand}.
#' @param sequences A named character vector (or Biostrings DNAStringSet) of
#'   full transcript sequences, keyed by isoform_id.
#' @param min_orf_nt Integer; minimum ORF length in nucleotides (default 9,
#'   i.e. 3 codons including start).
#' @param verbose Logical; if TRUE (default), print progress messages.
#' @return A tibble with one row per coding isoform and columns:
#' \describe{
#'   \item{isoform_id}{Isoform identifier}
#'   \item{gene_id}{Gene identifier}
#'   \item{strand}{Strand ("+"/"-")}
#'   \item{utr5_length}{5'UTR length in bp}
#'   \item{n_atg}{Total ATG count in 5'UTR}
#'   \item{n_stop_total}{Total stop codons (TAA+TAG+TGA) in 5'UTR}
#'   \item{n_taa, n_tag, n_tga}{Individual stop codon counts}
#'   \item{n_atg_frame0, n_atg_frame1, n_atg_frame2}{ATGs by frame relative to CDS}
#'   \item{n_stop_frame0, n_stop_frame1, n_stop_frame2}{Stops by frame relative to CDS}
#'   \item{n_orfs}{Qualifying ORFs (ATG + in-frame stop, >= min_orf_nt)}
#'   \item{n_orfs_inframe}{ORFs in-frame with main CDS (frame 0)}
#'   \item{n_orfs_outframe}{ORFs out-of-frame with main CDS}
#'   \item{n_orfs_overlapping}{ORFs whose stop falls at/past CDS start}
#'   \item{atg_density}{ATGs per 100 bp of 5'UTR}
#'   \item{stop_density}{Stops per 100 bp}
#'   \item{orf_density}{ORFs per 100 bp}
#'   \item{pct_utr5_in_orfs}{Fraction of 5'UTR occupied by ORFs}
#'   \item{n_strong_kozak_atg}{ATGs with strong Kozak (A/G at -3 AND G at +4)}
#'   \item{best_kozak_score}{Highest Kozak score among all uATGs (0-2)}
#'   \item{best_kozak_pwm_score}{Highest PWM-based Kozak score among all uATGs
#'     (continuous, from \code{\link{scoreKozakPWM}})}
#'   \item{first_atg_to_cds_nt}{Distance from most 5' ATG to CDS start (nt)}
#'   \item{max_overlap_nt}{For overlapping uORFs, max extension past CDS start}
#'   \item{longest_orf_nt}{Length of longest uORF}
#'   \item{total_orf_nt}{Total bp in all uORFs}
#'   \item{mean_orf_length}{Mean uORF length}
#'   \item{n_atg_without_orf}{ATGs not producing a qualifying ORF}
#'   \item{atg_validated}{TRUE if computed CDS start = ATG in FASTA}
#' }
#' @examples
#' \dontrun{
#' features <- scan5UtrFeatures(structures, cds, sequences)
#' }
#' @export
#' @importFrom tibble tibble
#' @importFrom stats setNames
scan5UtrFeatures <- function(structures, cds_metadata, sequences,
                             min_orf_nt = 9L, verbose = TRUE) {
  # Convert DNAStringSet to named character if needed
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- setNames(as.character(sequences), names(sequences))
  }

  # Filter to coding isoforms present in all three inputs
  coding <- cds_metadata[cds_metadata$coding_status == "coding", ]
  common_ids <- intersect(
    intersect(coding$isoform_id, structures$isoform_id),
    names(sequences)
  )

  if (verbose) message(sprintf("  Scanning 5'UTR features for %d isoforms",
                                length(common_ids)))

  # Build lookup tables
  struct_idx <- match(common_ids, structures$isoform_id)
  cds_idx <- match(common_ids, coding$isoform_id)

  # Pre-allocate result columns
  n <- length(common_ids)
  res <- list(
    isoform_id = common_ids,
    gene_id = structures$gene_id[struct_idx],
    strand = structures$strand[struct_idx],
    utr5_length = integer(n),
    n_atg = integer(n),
    n_stop_total = integer(n),
    n_taa = integer(n),
    n_tag = integer(n),
    n_tga = integer(n),
    n_atg_frame0 = integer(n),
    n_atg_frame1 = integer(n),
    n_atg_frame2 = integer(n),
    n_stop_frame0 = integer(n),
    n_stop_frame1 = integer(n),
    n_stop_frame2 = integer(n),
    n_orfs = integer(n),
    n_orfs_inframe = integer(n),
    n_orfs_outframe = integer(n),
    n_orfs_overlapping = integer(n),
    atg_density = numeric(n),
    stop_density = numeric(n),
    orf_density = numeric(n),
    pct_utr5_in_orfs = numeric(n),
    n_strong_kozak_atg = integer(n),
    best_kozak_score = integer(n),
    best_kozak_pwm_score = numeric(n),
    first_atg_to_cds_nt = integer(n),
    max_overlap_nt = integer(n),
    longest_orf_nt = integer(n),
    total_orf_nt = integer(n),
    mean_orf_length = numeric(n),
    n_atg_without_orf = integer(n),
    atg_validated = logical(n)
  )

  for (i in seq_len(n)) {
    si <- struct_idx[i]
    ci <- cds_idx[i]
    iso_id <- common_ids[i]

    str_i <- structures$strand[si]
    es <- structures$exon_starts[[si]]
    ee <- structures$exon_ends[[si]]
    cds_s <- coding$cds_start[ci]
    cds_e <- coding$cds_stop[ci]

    # 5'UTR length
    utr5_len <- .compute5UtrLength(str_i, es, ee, cds_s, cds_e)
    res$utr5_length[i] <- utr5_len

    # Get sequence
    full_seq <- toupper(sequences[iso_id])
    full_len <- nchar(full_seq)

    # ATG validation
    if (utr5_len > 0L && (utr5_len + 3L) <= full_len) {
      res$atg_validated[i] <- substr(full_seq, utr5_len + 1L,
                                     utr5_len + 3L) == "ATG"
    } else if (utr5_len == 0L && full_len >= 3L) {
      res$atg_validated[i] <- substr(full_seq, 1L, 3L) == "ATG"
    } else {
      res$atg_validated[i] <- NA
    }

    # Skip feature computation for zero-length 5'UTR
    if (utr5_len < 3L) {
      res$best_kozak_score[i] <- NA_integer_
      res$best_kozak_pwm_score[i] <- NA_real_
      res$first_atg_to_cds_nt[i] <- NA_integer_
      res$mean_orf_length[i] <- NA_real_
      next
    }

    utr5_seq <- substr(full_seq, 1L, utr5_len)

    # CDS ATG position in transcript space
    cds_atg_pos <- utr5_len + 1L

    # ----- Raw counts -----
    # Find all ATG positions in 5'UTR
    atg_matches <- gregexpr("ATG", utr5_seq, fixed = TRUE)[[1]]
    atg_positions <- if (atg_matches[1] == -1L) integer(0) else as.integer(atg_matches)
    res$n_atg[i] <- length(atg_positions)

    # Find all stop codons in 5'UTR
    taa_m <- gregexpr("TAA", utr5_seq, fixed = TRUE)[[1]]
    tag_m <- gregexpr("TAG", utr5_seq, fixed = TRUE)[[1]]
    tga_m <- gregexpr("TGA", utr5_seq, fixed = TRUE)[[1]]

    taa_pos <- if (taa_m[1] == -1L) integer(0) else as.integer(taa_m)
    tag_pos <- if (tag_m[1] == -1L) integer(0) else as.integer(tag_m)
    tga_pos <- if (tga_m[1] == -1L) integer(0) else as.integer(tga_m)

    res$n_taa[i] <- length(taa_pos)
    res$n_tag[i] <- length(tag_pos)
    res$n_tga[i] <- length(tga_pos)
    res$n_stop_total[i] <- length(taa_pos) + length(tag_pos) + length(tga_pos)

    # ----- Frame-stratified counts -----
    # Frame = (cds_atg_pos - codon_pos) %% 3
    # Frame 0 = in-frame with CDS
    if (length(atg_positions) > 0L) {
      atg_frames <- (cds_atg_pos - atg_positions) %% 3L
      res$n_atg_frame0[i] <- sum(atg_frames == 0L)
      res$n_atg_frame1[i] <- sum(atg_frames == 1L)
      res$n_atg_frame2[i] <- sum(atg_frames == 2L)
    }

    all_stop_pos <- c(taa_pos, tag_pos, tga_pos)
    if (length(all_stop_pos) > 0L) {
      stop_frames <- (cds_atg_pos - all_stop_pos) %% 3L
      res$n_stop_frame0[i] <- sum(stop_frames == 0L)
      res$n_stop_frame1[i] <- sum(stop_frames == 1L)
      res$n_stop_frame2[i] <- sum(stop_frames == 2L)
    }

    # ----- ORF detection -----
    n_orfs <- 0L
    n_inframe <- 0L
    n_outframe <- 0L
    n_overlapping <- 0L
    total_orf_bp <- 0L
    longest_orf <- 0L
    # Track uORF intervals within 5'UTR for coverage calculation (union, no double-counting)
    orf_intervals <- list()
    max_overlap <- 0L
    n_atg_with_orf <- 0L
    best_kozak <- NA_integer_
    best_kozak_pwm <- NA_real_
    n_strong <- 0L
    first_atg_dist <- NA_integer_

    for (atg_pos in atg_positions) {
      # Compute distance from first ATG to CDS
      if (is.na(first_atg_dist)) {
        first_atg_dist <- cds_atg_pos - atg_pos
      }

      # Walk in-frame from this ATG through full transcript
      found_stop <- FALSE
      stop_end <- NA_integer_
      codon_start <- atg_pos + 3L
      while (codon_start + 2L <= full_len) {
        codon <- substr(full_seq, codon_start, codon_start + 2L)
        if (codon %in% c("TAA", "TAG", "TGA")) {
          found_stop <- TRUE
          stop_end <- codon_start + 2L
          break
        }
        codon_start <- codon_start + 3L
      }

      if (!found_stop) next

      orf_len <- stop_end - atg_pos + 1L
      if (orf_len < min_orf_nt) next

      n_atg_with_orf <- n_atg_with_orf + 1L
      n_orfs <- n_orfs + 1L
      total_orf_bp <- total_orf_bp + orf_len
      # Record 5'UTR-clipped interval for union coverage calculation
      orf_end_clipped <- min(stop_end, utr5_len)
      orf_intervals[[length(orf_intervals) + 1L]] <- c(atg_pos, orf_end_clipped)
      if (orf_len > longest_orf) longest_orf <- orf_len

      # Frame relative to CDS
      frame <- (cds_atg_pos - atg_pos) %% 3L
      if (frame == 0L) {
        n_inframe <- n_inframe + 1L
      } else {
        n_outframe <- n_outframe + 1L
      }

      # Overlapping: stop extends past CDS start
      if (stop_end > utr5_len) {
        n_overlapping <- n_overlapping + 1L
        overlap_dist <- stop_end - utr5_len
        if (overlap_dist > max_overlap) max_overlap <- overlap_dist
      }

      # Kozak scoring (categorical 0-2)
      kz <- .scoreKozak(full_seq, atg_pos)
      if (!is.na(kz$score)) {
        if (is.na(best_kozak) || kz$score > best_kozak) {
          best_kozak <- kz$score
        }
        if (kz$score == 2L) n_strong <- n_strong + 1L
      }

      # PWM-based Kozak scoring (continuous)
      kz_pwm <- scoreKozakPWM(full_seq, atg_pos)
      if (!is.na(kz_pwm$score)) {
        if (is.na(best_kozak_pwm) || kz_pwm$score > best_kozak_pwm) {
          best_kozak_pwm <- kz_pwm$score
        }
      }
    }

    res$n_orfs[i] <- n_orfs
    res$n_orfs_inframe[i] <- n_inframe
    res$n_orfs_outframe[i] <- n_outframe
    res$n_orfs_overlapping[i] <- n_overlapping
    res$longest_orf_nt[i] <- longest_orf
    res$total_orf_nt[i] <- total_orf_bp
    res$mean_orf_length[i] <- if (n_orfs > 0L) total_orf_bp / n_orfs else NA_real_
    res$max_overlap_nt[i] <- max_overlap
    res$n_atg_without_orf[i] <- length(atg_positions) - n_atg_with_orf
    res$n_strong_kozak_atg[i] <- n_strong
    res$best_kozak_score[i] <- best_kozak
    res$best_kozak_pwm_score[i] <- best_kozak_pwm
    res$first_atg_to_cds_nt[i] <- first_atg_dist

    # ----- Density features -----
    res$atg_density[i] <- length(atg_positions) / utr5_len * 100
    res$stop_density[i] <- res$n_stop_total[i] / utr5_len * 100
    res$orf_density[i] <- n_orfs / utr5_len * 100

    # Compute coverage as union of intervals (no double-counting)
    if (length(orf_intervals) > 0L) {
      # Sort intervals by start position
      iv_mat <- do.call(rbind, orf_intervals)
      iv_mat <- iv_mat[order(iv_mat[, 1]), , drop = FALSE]
      # Merge overlapping intervals
      merged_bp <- 0L
      cur_start <- iv_mat[1, 1]
      cur_end <- iv_mat[1, 2]
      for (k in seq_len(nrow(iv_mat))[-1]) {
        if (iv_mat[k, 1] <= cur_end) {
          cur_end <- max(cur_end, iv_mat[k, 2])
        } else {
          merged_bp <- merged_bp + (cur_end - cur_start + 1L)
          cur_start <- iv_mat[k, 1]
          cur_end <- iv_mat[k, 2]
        }
      }
      merged_bp <- merged_bp + (cur_end - cur_start + 1L)
      res$pct_utr5_in_orfs[i] <- merged_bp / utr5_len
    } else {
      res$pct_utr5_in_orfs[i] <- 0
    }

    if (verbose && i %% 2000 == 0) {
      message(sprintf("    Processed %d / %d isoforms", i, n))
    }
  }

  tibble::as_tibble(res)
}


#' Convert Isoform Structures to GRangesList
#'
#' Converts an Isopair structures tibble to a \code{GRangesList} suitable for
#' use with ORFik and other Bioconductor tools. Each list element contains the
#' exon coordinates for one isoform.
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param isoform_ids Optional character vector of isoform IDs to include.
#' @return A named \code{GRangesList}, one element per isoform, with exon
#'   coordinates as \code{GRanges} entries.
#' @export
structuresToGRangesList <- function(structures, isoform_ids = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges is required for structuresToGRangesList()")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges is required for structuresToGRangesList()")
  }

  if (!is.null(isoform_ids)) {
    structures <- structures[structures$isoform_id %in% isoform_ids, ]
  }

  gr_list <- lapply(seq_len(nrow(structures)), function(i) {
    GenomicRanges::GRanges(
      seqnames = structures$chr[i],
      ranges = IRanges::IRanges(
        start = structures$exon_starts[[i]],
        end = structures$exon_ends[[i]]
      ),
      strand = structures$strand[i]
    )
  })

  names(gr_list) <- structures$isoform_id
  GenomicRanges::GRangesList(gr_list)
}


#' Detect Upstream Open Reading Frames in 5'UTR
#'
#' Finds all qualifying uORFs in the 5'UTR of coding isoforms. Each uORF is
#' defined by an ATG start codon in the 5'UTR and the first in-frame stop codon
#' (which may fall in or beyond the main CDS).
#'
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @param sequences A named character vector (or DNAStringSet) of full
#'   transcript sequences.
#' @param method Character; \code{"simple"} (default) for pure ATG-walk
#'   detection, or \code{"orfik"} to use ORFik's \code{findMapORFs()}.
#' @param min_length_codons Integer; minimum uORF length in codons including
#'   start codon (default 3).
#' @param kozak Logical; if TRUE (default), compute Kozak context scores.
#' @param verbose Logical; if TRUE (default), print progress.
#' @return A tibble with one row per uORF and columns:
#' \describe{
#'   \item{isoform_id, gene_id}{Identifiers}
#'   \item{uorf_id}{Sequential ID per isoform ("uORF_1", "uORF_2", ...)}
#'   \item{atg_transcript_pos, stop_transcript_pos}{Transcript coordinates}
#'   \item{atg_genomic_pos, stop_genomic_pos}{Genomic coordinates}
#'   \item{uorf_length_nt}{Total length in nucleotides (ATG through stop)}
#'   \item{is_overlapping}{Logical; stop at/past CDS start}
#'   \item{frame_relative_to_cds}{Phase (0/1/2) of uORF ATG relative to CDS}
#'   \item{kozak_score}{Kozak strength (0-2), NA if \code{kozak=FALSE}}
#'   \item{kozak_pwm_score}{Continuous PWM-based Kozak score (from
#'     \code{\link{scoreKozakPWM}}), NA if \code{kozak=FALSE}}
#'   \item{kozak_context}{7-nt context string around ATG}
#'   \item{uorf_sequence}{Full nucleotide sequence of the uORF}
#'   \item{stop_codon}{Which stop codon (TAA/TAG/TGA)}
#' }
#' @export
#' @importFrom tibble tibble
detectUorfs <- function(structures, cds_metadata, sequences,
                        method = c("simple", "orfik"),
                        min_length_codons = 3L, kozak = TRUE,
                        verbose = TRUE) {
  method <- match.arg(method)
  min_orf_nt <- min_length_codons * 3L

  if (method == "orfik") {
    if (!requireNamespace("ORFik", quietly = TRUE)) {
      stop("ORFik package is required for method='orfik'. ",
           "Install via BiocManager::install('ORFik') or use method='simple'.")
    }
    return(.detectUorfsOrfik(structures, cds_metadata, sequences,
                             min_length_codons, kozak, verbose))
  }

  # Convert DNAStringSet to named character if needed
  if (inherits(sequences, "DNAStringSet")) {
    sequences <- setNames(as.character(sequences), names(sequences))
  }

  coding <- cds_metadata[cds_metadata$coding_status == "coding", ]
  common_ids <- intersect(
    intersect(coding$isoform_id, structures$isoform_id),
    names(sequences)
  )

  if (verbose) message(sprintf("  Detecting uORFs (method='simple') for %d isoforms",
                                length(common_ids)))

  struct_idx <- match(common_ids, structures$isoform_id)
  cds_idx <- match(common_ids, coding$isoform_id)

  result_rows <- vector("list", length(common_ids))

  for (i in seq_along(common_ids)) {
    si <- struct_idx[i]
    ci <- cds_idx[i]
    iso_id <- common_ids[i]
    gene_id <- structures$gene_id[si]
    str_i <- structures$strand[si]
    es <- structures$exon_starts[[si]]
    ee <- structures$exon_ends[[si]]
    cds_s <- coding$cds_start[ci]
    cds_e <- coding$cds_stop[ci]

    utr5_len <- .compute5UtrLength(str_i, es, ee, cds_s, cds_e)
    if (utr5_len < 3L) next

    full_seq <- toupper(sequences[iso_id])
    full_len <- nchar(full_seq)
    utr5_seq <- substr(full_seq, 1L, utr5_len)
    cds_atg_pos <- utr5_len + 1L

    atg_matches <- gregexpr("ATG", utr5_seq, fixed = TRUE)[[1]]
    if (atg_matches[1] == -1L) next
    atg_positions <- as.integer(atg_matches)

    uorf_rows <- vector("list", length(atg_positions))
    uorf_count <- 0L

    for (atg_pos in atg_positions) {
      # Walk in-frame to first stop
      found_stop <- FALSE
      stop_codon_start <- NA_integer_
      stop_codon_type <- NA_character_
      codon_start <- atg_pos + 3L

      while (codon_start + 2L <= full_len) {
        codon <- substr(full_seq, codon_start, codon_start + 2L)
        if (codon %in% c("TAA", "TAG", "TGA")) {
          found_stop <- TRUE
          stop_codon_start <- codon_start
          stop_codon_type <- codon
          break
        }
        codon_start <- codon_start + 3L
      }

      if (!found_stop) next

      stop_end <- stop_codon_start + 2L
      orf_len <- stop_end - atg_pos + 1L
      if (orf_len < min_orf_nt) next

      uorf_count <- uorf_count + 1L
      frame <- (cds_atg_pos - atg_pos) %% 3L
      is_overlapping <- stop_end > utr5_len

      # Map to genomic coordinates
      atg_genomic <- .transcriptToGenomic(atg_pos, str_i, es, ee)
      stop_genomic <- .transcriptToGenomic(stop_end, str_i, es, ee)

      # Kozak
      kz_score <- NA_integer_
      kz_pwm_score <- NA_real_
      kz_context <- NA_character_
      if (kozak) {
        kz <- .scoreKozak(full_seq, atg_pos)
        kz_score <- kz$score
        kz_context <- kz$context
        kz_pwm <- scoreKozakPWM(full_seq, atg_pos)
        kz_pwm_score <- kz_pwm$score
      }

      uorf_seq <- substr(full_seq, atg_pos, stop_end)

      uorf_rows[[uorf_count]] <- data.frame(
        isoform_id = iso_id,
        gene_id = gene_id,
        uorf_id = paste0("uORF_", uorf_count),
        atg_transcript_pos = atg_pos,
        stop_transcript_pos = stop_end,
        atg_genomic_pos = atg_genomic,
        stop_genomic_pos = stop_genomic,
        uorf_length_nt = orf_len,
        is_overlapping = is_overlapping,
        frame_relative_to_cds = frame,
        kozak_score = kz_score,
        kozak_pwm_score = kz_pwm_score,
        kozak_context = kz_context,
        uorf_sequence = uorf_seq,
        stop_codon = stop_codon_type,
        stringsAsFactors = FALSE
      )
    }

    if (uorf_count > 0L) {
      result_rows[[i]] <- do.call(rbind, uorf_rows[seq_len(uorf_count)])
    }

    if (verbose && i %% 2000 == 0) {
      message(sprintf("    Processed %d / %d isoforms", i, length(common_ids)))
    }
  }

  result <- do.call(rbind, result_rows)
  if (is.null(result)) {
    return(tibble::tibble(
      isoform_id = character(0), gene_id = character(0),
      uorf_id = character(0),
      atg_transcript_pos = integer(0), stop_transcript_pos = integer(0),
      atg_genomic_pos = integer(0), stop_genomic_pos = integer(0),
      uorf_length_nt = integer(0), is_overlapping = logical(0),
      frame_relative_to_cds = integer(0),
      kozak_score = integer(0), kozak_pwm_score = numeric(0),
      kozak_context = character(0),
      uorf_sequence = character(0), stop_codon = character(0)
    ))
  }

  tibble::as_tibble(result)
}


#' ORFik-based uORF Detection (Internal)
#'
#' @keywords internal
.detectUorfsOrfik <- function(structures, cds_metadata, sequences,
                              min_length_codons, kozak, verbose) {
  # Placeholder for ORFik integration — Phase 2

  stop("ORFik method not yet implemented. Use method='simple'.")
}


#' Compare uORFs Between Reference and Comparator Isoforms
#'
#' For each isoform pair in profiles, performs set operations on uORFs: shared
#' (same ATG genomic position), gained (ATG in comparator only), lost (ATG in
#' reference only). For shared uORFs, checks whether the stop codon position
#' differs (indicating a frameshift between ATG and stop).
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with columns
#'   \code{gene_id}, \code{reference_isoform_id}, \code{comparator_isoform_id}.
#' @param uorf_table A tibble from \code{\link{detectUorfs}()}.
#' @param match_tolerance_bp Integer; tolerance for ATG genomic position
#'   matching (default 0 = exact match).
#' @return A list with two elements:
#' \describe{
#'   \item{pair_summary}{Tibble (one row per pair) with counts of shared,
#'     gained, lost uORFs}
#'   \item{uorf_detail}{Tibble (one row per uORF) with status and per-isoform
#'     properties}
#' }
#' @export
#' @importFrom tibble tibble
compareUorfs <- function(profiles, uorf_table, match_tolerance_bp = 0L) {
  # Get unique pairs
  pairs <- unique(profiles[, c("gene_id", "reference_isoform_id",
                                "comparator_isoform_id")])

  pair_summaries <- vector("list", nrow(pairs))
  uorf_details <- vector("list", nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    gene <- pairs$gene_id[i]
    ref_id <- pairs$reference_isoform_id[i]
    comp_id <- pairs$comparator_isoform_id[i]

    ref_uorfs <- uorf_table[uorf_table$isoform_id == ref_id, ]
    comp_uorfs <- uorf_table[uorf_table$isoform_id == comp_id, ]

    n_ref <- nrow(ref_uorfs)
    n_comp <- nrow(comp_uorfs)

    # Match by ATG genomic position
    if (n_ref == 0L && n_comp == 0L) {
      pair_summaries[[i]] <- data.frame(
        gene_id = gene,
        reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        n_ref_uorfs = 0L, n_comp_uorfs = 0L,
        n_shared = 0L, n_gained = 0L, n_lost = 0L,
        n_shared_same_stop = 0L, n_shared_diff_stop = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    ref_matched <- rep(FALSE, n_ref)
    comp_matched <- rep(FALSE, n_comp)
    shared_rows <- list()

    if (n_ref > 0L && n_comp > 0L) {
      for (ri in seq_len(n_ref)) {
        for (ci in seq_len(n_comp)) {
          if (comp_matched[ci]) next
          if (abs(ref_uorfs$atg_genomic_pos[ri] -
                  comp_uorfs$atg_genomic_pos[ci]) <= match_tolerance_bp) {
            ref_matched[ri] <- TRUE
            comp_matched[ci] <- TRUE

            stop_changed <- !identical(ref_uorfs$stop_genomic_pos[ri],
                                       comp_uorfs$stop_genomic_pos[ci])

            shared_rows <- c(shared_rows, list(data.frame(
              gene_id = gene,
              reference_isoform_id = ref_id,
              comparator_isoform_id = comp_id,
              atg_genomic_pos = ref_uorfs$atg_genomic_pos[ri],
              status = "shared",
              ref_stop_genomic_pos = ref_uorfs$stop_genomic_pos[ri],
              comp_stop_genomic_pos = comp_uorfs$stop_genomic_pos[ci],
              ref_uorf_length_nt = ref_uorfs$uorf_length_nt[ri],
              comp_uorf_length_nt = comp_uorfs$uorf_length_nt[ci],
              ref_is_overlapping = ref_uorfs$is_overlapping[ri],
              comp_is_overlapping = comp_uorfs$is_overlapping[ci],
              ref_frame_to_cds = ref_uorfs$frame_relative_to_cds[ri],
              comp_frame_to_cds = comp_uorfs$frame_relative_to_cds[ci],
              stop_changed = stop_changed,
              stringsAsFactors = FALSE
            )))
            break
          }
        }
      }
    }

    # Lost uORFs (in ref but not comp)
    lost_rows <- list()
    if (n_ref > 0L) {
      for (ri in which(!ref_matched)) {
        lost_rows <- c(lost_rows, list(data.frame(
          gene_id = gene,
          reference_isoform_id = ref_id,
          comparator_isoform_id = comp_id,
          atg_genomic_pos = ref_uorfs$atg_genomic_pos[ri],
          status = "lost",
          ref_stop_genomic_pos = ref_uorfs$stop_genomic_pos[ri],
          comp_stop_genomic_pos = NA_integer_,
          ref_uorf_length_nt = ref_uorfs$uorf_length_nt[ri],
          comp_uorf_length_nt = NA_integer_,
          ref_is_overlapping = ref_uorfs$is_overlapping[ri],
          comp_is_overlapping = NA,
          ref_frame_to_cds = ref_uorfs$frame_relative_to_cds[ri],
          comp_frame_to_cds = NA_integer_,
          stop_changed = NA,
          stringsAsFactors = FALSE
        )))
      }
    }

    # Gained uORFs (in comp but not ref)
    gained_rows <- list()
    if (n_comp > 0L) {
      for (ci in which(!comp_matched)) {
        gained_rows <- c(gained_rows, list(data.frame(
          gene_id = gene,
          reference_isoform_id = ref_id,
          comparator_isoform_id = comp_id,
          atg_genomic_pos = comp_uorfs$atg_genomic_pos[ci],
          status = "gained",
          ref_stop_genomic_pos = NA_integer_,
          comp_stop_genomic_pos = comp_uorfs$stop_genomic_pos[ci],
          ref_uorf_length_nt = NA_integer_,
          comp_uorf_length_nt = comp_uorfs$uorf_length_nt[ci],
          ref_is_overlapping = NA,
          comp_is_overlapping = comp_uorfs$is_overlapping[ci],
          ref_frame_to_cds = NA_integer_,
          comp_frame_to_cds = comp_uorfs$frame_relative_to_cds[ci],
          stop_changed = NA,
          stringsAsFactors = FALSE
        )))
      }
    }

    all_detail <- do.call(rbind, c(shared_rows, lost_rows, gained_rows))
    uorf_details[[i]] <- all_detail

    n_shared <- sum(ref_matched)
    n_shared_same <- if (length(shared_rows) > 0L) {
      sum(!vapply(shared_rows, function(x) x$stop_changed, logical(1)))
    } else 0L
    n_shared_diff <- n_shared - n_shared_same

    pair_summaries[[i]] <- data.frame(
      gene_id = gene,
      reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      n_ref_uorfs = n_ref,
      n_comp_uorfs = n_comp,
      n_shared = n_shared,
      n_gained = sum(!comp_matched),
      n_lost = sum(!ref_matched),
      n_shared_same_stop = n_shared_same,
      n_shared_diff_stop = n_shared_diff,
      stringsAsFactors = FALSE
    )
  }

  pair_summary <- do.call(rbind, pair_summaries)
  uorf_detail <- do.call(rbind, uorf_details)

  if (is.null(pair_summary)) {
    pair_summary <- tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      n_ref_uorfs = integer(0), n_comp_uorfs = integer(0),
      n_shared = integer(0), n_gained = integer(0), n_lost = integer(0),
      n_shared_same_stop = integer(0), n_shared_diff_stop = integer(0)
    )
  }
  if (is.null(uorf_detail)) {
    uorf_detail <- tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      atg_genomic_pos = integer(0), status = character(0),
      ref_stop_genomic_pos = integer(0), comp_stop_genomic_pos = integer(0),
      ref_uorf_length_nt = integer(0), comp_uorf_length_nt = integer(0),
      ref_is_overlapping = logical(0), comp_is_overlapping = logical(0),
      ref_frame_to_cds = integer(0), comp_frame_to_cds = integer(0),
      stop_changed = logical(0)
    )
  }

  list(
    pair_summary = tibble::as_tibble(pair_summary),
    uorf_detail = tibble::as_tibble(uorf_detail)
  )
}


#' Attribute uORF Changes to Splice Events
#'
#' For each gained, lost, or stop-changed uORF, identifies the causal splice
#' event using a frameshift-first attribution hierarchy paralleling PTC
#' attribution.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with a
#'   \code{detailed_events} list-column.
#' @param uorf_comparison A list from \code{\link{compareUorfs}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @return A tibble with one row per attributed uORF-event pair:
#' \describe{
#'   \item{gene_id, reference_isoform_id, comparator_isoform_id}{Identifiers}
#'   \item{atg_genomic_pos}{ATG genomic position}
#'   \item{uorf_status}{"gained", "lost", "stop_changed", "frame_shifted"}
#'   \item{attribution_type}{"frameshift", "containment",
#'     "frame_relationship_shift", "inframe_splice", "unresolved"}
#'   \item{event_type, event_direction, event_five_prime, event_three_prime,
#'     event_bp_diff}{Event details from profiles}
#' }
#' @export
#' @importFrom tibble tibble as_tibble
attributeUorfsToEvents <- function(profiles, uorf_comparison, cds_metadata) {
  detail <- uorf_comparison$uorf_detail

  # Filter to actionable uORFs
  actionable <- detail[
    detail$status %in% c("gained", "lost") |
    (!is.na(detail$stop_changed) & detail$stop_changed) |
    (!is.na(detail$ref_frame_to_cds) & !is.na(detail$comp_frame_to_cds) &
     detail$ref_frame_to_cds != detail$comp_frame_to_cds),
  ]

  if (nrow(actionable) == 0L) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      atg_genomic_pos = integer(0), uorf_status = character(0),
      attribution_type = character(0),
      event_type = character(0), event_direction = character(0),
      event_five_prime = integer(0), event_three_prime = integer(0),
      event_bp_diff = integer(0)
    ))
  }

  result_rows <- vector("list", nrow(actionable))

  for (i in seq_len(nrow(actionable))) {
    row <- actionable[i, ]
    gene <- row$gene_id
    ref_id <- row$reference_isoform_id
    comp_id <- row$comparator_isoform_id
    atg_pos <- row$atg_genomic_pos

    # Classify uORF status
    uorf_status <- if (row$status == "gained") {
      "gained"
    } else if (row$status == "lost") {
      "lost"
    } else if (!is.na(row$stop_changed) && row$stop_changed) {
      "stop_changed"
    } else {
      "frame_shifted"
    }

    # Get events for this pair
    pair_row <- profiles[
      profiles$gene_id == gene &
      profiles$reference_isoform_id == ref_id &
      profiles$comparator_isoform_id == comp_id,
    ]

    if (nrow(pair_row) == 0L || is.null(pair_row$detailed_events[[1]])) {
      result_rows[[i]] <- data.frame(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        atg_genomic_pos = atg_pos, uorf_status = uorf_status,
        attribution_type = "unresolved",
        event_type = NA_character_, event_direction = NA_character_,
        event_five_prime = NA_integer_, event_three_prime = NA_integer_,
        event_bp_diff = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    events <- pair_row$detailed_events[[1]]
    if (nrow(events) == 0L) {
      result_rows[[i]] <- data.frame(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        atg_genomic_pos = atg_pos, uorf_status = uorf_status,
        attribution_type = "unresolved",
        event_type = NA_character_, event_direction = NA_character_,
        event_five_prime = NA_integer_, event_three_prime = NA_integer_,
        event_bp_diff = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Stop position region (use whichever side has data)
    stop_pos <- if (!is.na(row$ref_stop_genomic_pos)) {
      row$ref_stop_genomic_pos
    } else if (!is.na(row$comp_stop_genomic_pos)) {
      row$comp_stop_genomic_pos
    } else {
      NA_integer_
    }

    # Get CDS ATG genomic position for frame-relationship attribution
    cds_row <- cds_metadata[cds_metadata$isoform_id == ref_id, ]
    cds_atg_genomic <- if (nrow(cds_row) > 0L) {
      str_i <- cds_row$strand[1]
      if (str_i == "+") cds_row$cds_start[1] else cds_row$cds_stop[1]
    } else {
      NA_integer_
    }

    attributed <- FALSE
    attr_type <- "unresolved"
    attr_event_idx <- NA_integer_

    # --- Priority 1: Frameshift ---
    # For stop_changed or gained/lost with shared ATG: event between ATG and
    # stop with bp_diff %% 3 != 0
    if (!attributed && !is.na(stop_pos)) {
      region_min <- min(atg_pos, stop_pos, na.rm = TRUE)
      region_max <- max(atg_pos, stop_pos, na.rm = TRUE)

      for (ei in seq_len(nrow(events))) {
        ev_5p <- events$five_prime[ei]
        ev_3p <- events$three_prime[ei]
        ev_min <- min(ev_5p, ev_3p)
        ev_max <- max(ev_5p, ev_3p)
        bp_diff <- events$bp_diff[ei]

        if (ev_min >= region_min && ev_max <= region_max &&
            !is.na(bp_diff) && bp_diff %% 3L != 0L) {
          attr_type <- "frameshift"
          attr_event_idx <- ei
          attributed <- TRUE
          break
        }
      }
    }

    # --- Priority 2: Containment ---
    # For gained uORFs: event boundaries contain the ATG
    if (!attributed && uorf_status == "gained") {
      for (ei in seq_len(nrow(events))) {
        ev_5p <- events$five_prime[ei]
        ev_3p <- events$three_prime[ei]
        ev_min <- min(ev_5p, ev_3p)
        ev_max <- max(ev_5p, ev_3p)

        if (atg_pos >= ev_min && atg_pos <= ev_max) {
          attr_type <- "containment"
          attr_event_idx <- ei
          attributed <- TRUE
          break
        }
      }
    }

    # --- Priority 2b: Containment for lost uORFs ---
    if (!attributed && uorf_status == "lost") {
      for (ei in seq_len(nrow(events))) {
        ev_5p <- events$five_prime[ei]
        ev_3p <- events$three_prime[ei]
        ev_min <- min(ev_5p, ev_3p)
        ev_max <- max(ev_5p, ev_3p)

        if (atg_pos >= ev_min && atg_pos <= ev_max) {
          attr_type <- "containment"
          attr_event_idx <- ei
          attributed <- TRUE
          break
        }
      }
    }

    # --- Priority 3: Frame-relationship shift ---
    # Shared uORFs where frame relative to CDS changed
    if (!attributed && uorf_status == "frame_shifted" &&
        !is.na(cds_atg_genomic)) {
      region_min <- min(atg_pos, cds_atg_genomic)
      region_max <- max(atg_pos, cds_atg_genomic)

      for (ei in seq_len(nrow(events))) {
        ev_5p <- events$five_prime[ei]
        ev_3p <- events$three_prime[ei]
        ev_min <- min(ev_5p, ev_3p)
        ev_max <- max(ev_5p, ev_3p)

        if (ev_min >= region_min && ev_max <= region_max) {
          attr_type <- "frame_relationship_shift"
          attr_event_idx <- ei
          attributed <- TRUE
          break
        }
      }
    }

    # --- Priority 4: In-frame splice ---
    if (!attributed && !is.na(stop_pos)) {
      region_min <- min(atg_pos, stop_pos, na.rm = TRUE)
      region_max <- max(atg_pos, stop_pos, na.rm = TRUE)

      for (ei in seq_len(nrow(events))) {
        ev_5p <- events$five_prime[ei]
        ev_3p <- events$three_prime[ei]
        ev_min <- min(ev_5p, ev_3p)
        ev_max <- max(ev_5p, ev_3p)
        bp_diff <- events$bp_diff[ei]

        if (ev_min >= region_min && ev_max <= region_max &&
            !is.na(bp_diff) && bp_diff %% 3L == 0L) {
          attr_type <- "inframe_splice"
          attr_event_idx <- ei
          attributed <- TRUE
          break
        }
      }
    }

    # Build result row
    if (!is.na(attr_event_idx)) {
      ev <- events[attr_event_idx, ]
      result_rows[[i]] <- data.frame(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        atg_genomic_pos = atg_pos, uorf_status = uorf_status,
        attribution_type = attr_type,
        event_type = ev$event_type,
        event_direction = ev$direction,
        event_five_prime = ev$five_prime,
        event_three_prime = ev$three_prime,
        event_bp_diff = ev$bp_diff,
        stringsAsFactors = FALSE
      )
    } else {
      result_rows[[i]] <- data.frame(
        gene_id = gene, reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        atg_genomic_pos = atg_pos, uorf_status = uorf_status,
        attribution_type = "unresolved",
        event_type = NA_character_, event_direction = NA_character_,
        event_five_prime = NA_integer_, event_three_prime = NA_integer_,
        event_bp_diff = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, result_rows)
  if (is.null(result)) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      atg_genomic_pos = integer(0), uorf_status = character(0),
      attribution_type = character(0),
      event_type = character(0), event_direction = character(0),
      event_five_prime = integer(0), event_three_prime = integer(0),
      event_bp_diff = integer(0)
    ))
  }

  tibble::as_tibble(result)
}
