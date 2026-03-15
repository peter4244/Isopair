#' Map Splice Events to Protein-Level Affected Regions
#'
#' For each isoform pair, maps each splice event's genomic coordinates to
#' protein-level affected regions using the reference isoform's CDS and exon
#' structure. Events entirely in UTR regions are excluded. Events that
#' partially overlap the CDS are clipped to CDS boundaries.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with columns
#'   \code{reference_isoform_id}, \code{comparator_isoform_id},
#'   \code{gene_id}, \code{detailed_events} (list column of event data frames).
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()} with
#'   columns \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{strand}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()} with
#'   columns \code{isoform_id}, \code{exon_starts} (list), \code{exon_ends}
#'   (list), \code{strand}.
#' @param buffer_aa Integer; buffer in amino acids around event boundaries for
#'   divergence confirmation (default 5).
#' @return A tibble with one row per CDS-overlapping event per pair. Columns:
#'   \code{gene_id}, \code{reference_isoform_id}, \code{comparator_isoform_id},
#'   \code{event_type}, \code{direction}, \code{genomic_start},
#'   \code{genomic_end}, \code{cds_nt_start}, \code{cds_nt_end},
#'   \code{protein_start_aa}, \code{protein_end_aa},
#'   \code{protein_start_aa_buffered}, \code{protein_end_aa_buffered},
#'   \code{event_affects_frame}. Returns a zero-row tibble if no events
#'   overlap any CDS.
#' @examples
#' # Build a simple + strand example
#' profiles <- tibble::tibble(
#'   gene_id = "gene1",
#'   reference_isoform_id = "tx1",
#'   comparator_isoform_id = "tx2",
#'   detailed_events = list(data.frame(
#'     event_type = "SE", direction = "LOSS",
#'     five_prime = 350L, three_prime = 450L,
#'     bp_diff = 101L, stringsAsFactors = FALSE
#'   ))
#' )
#' cds_metadata <- tibble::tibble(
#'   isoform_id = "tx1", coding_status = "coding",
#'   cds_start = 200L, cds_stop = 800L, orf_length = 400L,
#'   strand = "+"
#' )
#' structures <- tibble::tibble(
#'   isoform_id = "tx1", gene_id = "gene1", chr = "chr1",
#'   strand = "+", n_exons = 3L,
#'   exon_starts = list(c(100L, 300L, 600L)),
#'   exon_ends = list(c(250L, 500L, 900L)),
#'   tx_start = 100L, tx_end = 900L, n_junctions = 2L
#' )
#' result <- mapSpliceToProtein(profiles, cds_metadata, structures)
#' result
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
mapSpliceToProtein <- function(profiles, cds_metadata, structures,
                               buffer_aa = 5L) {

  # --- Input validation ---
  if (!is.data.frame(profiles)) {
    stop("'profiles' must be a data frame or tibble.", call. = FALSE)
  }
  required_prof <- c("gene_id", "reference_isoform_id",
                     "comparator_isoform_id", "detailed_events")
  missing_prof <- setdiff(required_prof, names(profiles))
  if (length(missing_prof) > 0L) {
    stop(sprintf("'profiles' is missing required columns: %s",
                 paste(missing_prof, collapse = ", ")), call. = FALSE)
  }

  if (!is.data.frame(cds_metadata)) {
    stop("'cds_metadata' must be a data frame or tibble.", call. = FALSE)
  }
  required_cds <- c("isoform_id", "coding_status", "cds_start",
                    "cds_stop", "strand")
  missing_cds <- setdiff(required_cds, names(cds_metadata))
  if (length(missing_cds) > 0L) {
    stop(sprintf("'cds_metadata' is missing required columns: %s",
                 paste(missing_cds, collapse = ", ")), call. = FALSE)
  }

  if (!is.data.frame(structures)) {
    stop("'structures' must be a data frame or tibble.", call. = FALSE)
  }
  required_str <- c("isoform_id", "exon_starts", "exon_ends", "strand")
  missing_str <- setdiff(required_str, names(structures))
  if (length(missing_str) > 0L) {
    stop(sprintf("'structures' is missing required columns: %s",
                 paste(missing_str, collapse = ", ")), call. = FALSE)
  }

  buffer_aa <- as.integer(buffer_aa)
  if (is.na(buffer_aa) || buffer_aa < 0L) {
    stop("'buffer_aa' must be a non-negative integer.", call. = FALSE)
  }

  # --- Build lookups ---
  coding_cds <- cds_metadata[cds_metadata$coding_status == "coding", ]
  cds_lookup <- split(coding_cds, coding_cds$isoform_id)
  struct_lookup <- split(structures, structures$isoform_id)

  # --- Empty output template ---
  empty_result <- tibble::tibble(
    gene_id = character(0),
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    event_type = character(0),
    direction = character(0),
    genomic_start = integer(0),
    genomic_end = integer(0),
    cds_nt_start = integer(0),
    cds_nt_end = integer(0),
    protein_start_aa = integer(0),
    protein_end_aa = integer(0),
    protein_start_aa_buffered = integer(0),
    protein_end_aa_buffered = integer(0),
    event_affects_frame = logical(0)
  )

  if (nrow(profiles) == 0L) return(empty_result)

  all_results <- list()

  for (i in seq_len(nrow(profiles))) {
    gene <- profiles$gene_id[i]
    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]
    events <- profiles$detailed_events[[i]]

    # Skip if no events
    if (is.null(events) || nrow(events) == 0L) next

    # Skip if reference is not coding
    ref_cds <- cds_lookup[[ref_id]]
    if (is.null(ref_cds) || nrow(ref_cds) == 0L) next

    # Get reference structure
    ref_struct <- struct_lookup[[ref_id]]
    if (is.null(ref_struct) || nrow(ref_struct) == 0L) next

    cds_start <- ref_cds$cds_start[1L]
    cds_stop <- ref_cds$cds_stop[1L]
    strand <- ref_cds$strand[1L]
    exon_starts <- ref_struct$exon_starts[[1L]]
    exon_ends <- ref_struct$exon_ends[[1L]]

    # Compute total CDS length in nt (for clamping AA positions)
    cds_nt_total <- .computeCdsLength(exon_starts, exon_ends,
                                       cds_start, cds_stop)
    if (cds_nt_total <= 0L) next
    max_aa <- as.integer(ceiling(cds_nt_total / 3))

    for (j in seq_len(nrow(events))) {
      evt <- events[j, ]

      # Genomic boundaries (always start < end regardless of strand)
      genomic_start <- as.integer(pmin(evt$five_prime, evt$three_prime))
      genomic_end <- as.integer(pmax(evt$five_prime, evt$three_prime))

      # Check CDS overlap
      overlap_start <- max(genomic_start, cds_start)
      overlap_end <- min(genomic_end, cds_stop)
      if (overlap_start > overlap_end) next  # entirely in UTR

      # Convert clipped boundaries to CDS-relative nt positions
      nt_start <- .genomicToCdsPosition(overlap_start, exon_starts,
                                         exon_ends, cds_start, cds_stop,
                                         strand)
      nt_end <- .genomicToCdsPosition(overlap_end, exon_starts,
                                       exon_ends, cds_start, cds_stop,
                                       strand)

      # Skip if either position falls outside exons (intronic)
      if (is.na(nt_start) || is.na(nt_end)) next

      # Ensure nt_start <= nt_end in CDS space
      cds_nt_start <- min(nt_start, nt_end)
      cds_nt_end <- max(nt_start, nt_end)

      # Convert to AA positions
      protein_start_aa <- as.integer(ceiling(cds_nt_start / 3))
      protein_end_aa <- as.integer(ceiling(cds_nt_end / 3))

      # Apply buffer (clamp to valid range)
      protein_start_buffered <- max(1L, protein_start_aa - buffer_aa)
      protein_end_buffered <- min(max_aa, protein_end_aa + buffer_aa)

      # Does this event affect reading frame?
      bp_diff <- if (!is.null(evt$bp_diff) && !is.na(evt$bp_diff)) {
        evt$bp_diff
      } else {
        abs(genomic_end - genomic_start) + 1L
      }
      event_affects_frame <- (bp_diff %% 3L) != 0L

      all_results[[length(all_results) + 1L]] <- tibble::tibble(
        gene_id = gene,
        reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        event_type = evt$event_type,
        direction = evt$direction,
        genomic_start = genomic_start,
        genomic_end = genomic_end,
        cds_nt_start = as.integer(cds_nt_start),
        cds_nt_end = as.integer(cds_nt_end),
        protein_start_aa = protein_start_aa,
        protein_end_aa = protein_end_aa,
        protein_start_aa_buffered = protein_start_buffered,
        protein_end_aa_buffered = protein_end_buffered,
        event_affects_frame = event_affects_frame
      )
    }
  }

  if (length(all_results) == 0L) return(empty_result)
  dplyr::bind_rows(all_results)
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Convert a Genomic Position to CDS-Relative Nucleotide Position
#'
#' Walks through exons in transcript order (5' to 3') counting only exonic
#' bases between the translation start and the target genomic position.
#' Returns a 1-based CDS-relative nucleotide position, or NA if the
#' position falls outside exonic regions of the CDS.
#'
#' @param genomic_pos Integer; the genomic coordinate to convert.
#' @param exon_starts Integer vector of exon start coordinates (1-based).
#' @param exon_ends Integer vector of exon end coordinates (1-based, closed).
#' @param cds_start Integer; min genomic CDS coordinate.
#' @param cds_stop Integer; max genomic CDS coordinate.
#' @param strand Character; "+" or "-".
#' @return Integer CDS-relative nucleotide position (1-based), or NA if the
#'   position is not within a CDS exon.
#' @keywords internal
.genomicToCdsPosition <- function(genomic_pos, exon_starts, exon_ends,
                                   cds_start, cds_stop, strand) {

  # Position must be within CDS genomic range
  if (genomic_pos < cds_start || genomic_pos > cds_stop) return(NA_integer_)

  if (strand == "+") {
    # Walk exons in ascending order (5' to 3')
    ord <- order(exon_starts)
    cumulative <- 0L

    for (idx in ord) {
      es <- exon_starts[idx]
      ee <- exon_ends[idx]

      # Clip exon to CDS boundaries
      exon_cds_start <- max(es, cds_start)
      exon_cds_end <- min(ee, cds_stop)

      # Skip exons with no CDS overlap
      if (exon_cds_start > exon_cds_end) next

      if (genomic_pos >= exon_cds_start && genomic_pos <= exon_cds_end) {
        # Target is in this exon
        cumulative <- cumulative + as.integer(genomic_pos - exon_cds_start) + 1L
        return(as.integer(cumulative))
      }

      # Entire CDS portion of this exon is upstream
      cumulative <- cumulative + as.integer(exon_cds_end - exon_cds_start) + 1L
    }

  } else {
    # Minus strand: walk exons in descending order (5' to 3')
    ord <- order(exon_starts, decreasing = TRUE)
    cumulative <- 0L

    for (idx in ord) {
      es <- exon_starts[idx]
      ee <- exon_ends[idx]

      # Clip exon to CDS boundaries
      exon_cds_start <- max(es, cds_start)
      exon_cds_end <- min(ee, cds_stop)

      # Skip exons with no CDS overlap
      if (exon_cds_start > exon_cds_end) next

      if (genomic_pos >= exon_cds_start && genomic_pos <= exon_cds_end) {
        # Target is in this exon — on minus strand, 5' is at high coordinates
        cumulative <- cumulative + as.integer(exon_cds_end - genomic_pos) + 1L
        return(as.integer(cumulative))
      }

      # Entire CDS portion of this exon is upstream (in transcript order)
      cumulative <- cumulative + as.integer(exon_cds_end - exon_cds_start) + 1L
    }
  }

  # Position not found in any CDS exon (intronic)
  NA_integer_
}


#' Compute Total CDS Length in Nucleotides
#'
#' Sums up exonic bases that fall within the CDS boundaries.
#'
#' @param exon_starts Integer vector of exon starts.
#' @param exon_ends Integer vector of exon ends.
#' @param cds_start Integer; min genomic CDS coordinate.
#' @param cds_stop Integer; max genomic CDS coordinate.
#' @return Integer total CDS length in nucleotides.
#' @keywords internal
.computeCdsLength <- function(exon_starts, exon_ends, cds_start, cds_stop) {
  total <- 0L
  for (i in seq_along(exon_starts)) {
    overlap_start <- max(exon_starts[i], cds_start)
    overlap_end <- min(exon_ends[i], cds_stop)
    if (overlap_start <= overlap_end) {
      total <- total + as.integer(overlap_end - overlap_start) + 1L
    }
  }
  total
}
