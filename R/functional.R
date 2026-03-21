#' Map Events to Genomic Regions
#'
#' Maps each splicing event to overlapping genomic regions (5'UTR, CDS,
#' 3'UTR, etc.) using annotated union exon mappings. Events can map to
#' multiple region types when they span region boundaries.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param annotated_ue A tibble from \code{\link{annotateRegionTypes}()}
#'   with columns: \code{isoform_id}, \code{union_exon_id},
#'   \code{isoform_exon_start}, \code{isoform_exon_end},
#'   \code{region_type}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with \code{isoform_id}, \code{cds_start}, \code{cds_stop}.
#' @return A tibble with columns: \code{event_row_id}, \code{gene_id},
#'   \code{reference_isoform_id}, \code{comparator_isoform_id},
#'   \code{event_type}, \code{direction}, \code{genomic_start},
#'   \code{genomic_end}, \code{region_type}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' er <- mapEventsToRegions(example_profiles, annotated, cds)
#' head(er)
#' @export
#' @importFrom dplyr bind_rows distinct
#' @importFrom tibble tibble
mapEventsToRegions <- function(profiles, annotated_ue, cds_metadata) {

  # Filter annotated_ue to reference isoforms in profiles
  ref_ids <- unique(profiles$reference_isoform_id)
  ref_ue <- annotated_ue[annotated_ue$isoform_id %in% ref_ids, ]

  all_mapped <- list()
  event_counter <- 0L

  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next

    ref_id <- profiles$reference_isoform_id[i]
    gene_id <- profiles$gene_id[i]
    comp_id <- profiles$comparator_isoform_id[i]

    # Get annotated union exons for this reference isoform
    ref_annotations <- ref_ue[ref_ue$isoform_id == ref_id, ]
    if (nrow(ref_annotations) == 0L) next

    for (j in seq_len(nrow(events))) {
      event_counter <- event_counter + 1L
      evt <- events[j, ]

      genomic_start <- pmin(evt$five_prime, evt$three_prime)
      genomic_end <- pmax(evt$five_prime, evt$three_prime)

      # Find overlapping union exon annotations
      overlaps <- ref_annotations[
        ref_annotations$isoform_exon_start <= genomic_end &
          ref_annotations$isoform_exon_end >= genomic_start, ]

      if (nrow(overlaps) > 0L) {
        regions <- unique(overlaps$region_type)
        for (r in regions) {
          all_mapped[[length(all_mapped) + 1L]] <- tibble::tibble(
            event_row_id = event_counter,
            gene_id = gene_id,
            reference_isoform_id = ref_id,
            comparator_isoform_id = comp_id,
            event_type = evt$event_type,
            direction = evt$direction,
            genomic_start = as.integer(genomic_start),
            genomic_end = as.integer(genomic_end),
            region_type = r
          )
        }
      } else {
        # No overlap with annotated exons → intronic
        all_mapped[[length(all_mapped) + 1L]] <- tibble::tibble(
          event_row_id = event_counter,
          gene_id = gene_id,
          reference_isoform_id = ref_id,
          comparator_isoform_id = comp_id,
          event_type = evt$event_type,
          direction = evt$direction,
          genomic_start = as.integer(genomic_start),
          genomic_end = as.integer(genomic_end),
          region_type = "intronic"
        )
      }
    }
  }

  if (length(all_mapped) == 0L) {
    return(tibble::tibble(
      event_row_id = integer(0), gene_id = character(0),
      reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      event_type = character(0), direction = character(0),
      genomic_start = integer(0), genomic_end = integer(0),
      region_type = character(0)
    ))
  }

  dplyr::bind_rows(all_mapped)
}


#' Quantify Regional Impact of Each Splicing Event
#'
#' For each event in each profile, computes exact basepair overlap with
#' 5'UTR, CDS, 3'UTR, and intronic regions from both the reference and
#' comparator isoform perspectives. Boundary exons (contains_orf_start,
#' contains_orf_stop, contains_orf_start_stop) are decomposed into their
#' UTR and CDS components using CDS coordinates.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()} with
#'   \code{detailed_events}.
#' @param annotated_ue A tibble from \code{\link{annotateRegionTypes}()}
#'   with columns: \code{isoform_id}, \code{isoform_exon_start},
#'   \code{isoform_exon_end}, \code{region_type}, \code{cds_start},
#'   \code{cds_stop}, \code{strand}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with \code{isoform_id}, \code{cds_start}, \code{cds_stop},
#'   \code{strand}.
#' @return A tibble with one row per event, columns: \code{event_id},
#'   \code{gene_id}, \code{reference_isoform_id},
#'   \code{comparator_isoform_id}, \code{event_type}, \code{direction},
#'   \code{genomic_start}, \code{genomic_end}, \code{event_span_bp},
#'   \code{ref_utr5_bp}, \code{ref_cds_bp}, \code{ref_utr3_bp},
#'   \code{ref_intronic_bp}, \code{comp_utr5_bp}, \code{comp_cds_bp},
#'   \code{comp_utr3_bp}, \code{comp_intronic_bp}, \code{ref_category},
#'   \code{comp_category}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' ri <- quantifyRegionalImpact(example_profiles, annotated, cds)
#' head(ri)
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
quantifyRegionalImpact <- function(profiles, annotated_ue, cds_metadata) {

  # Pre-compute lookup: split annotated_ue by isoform_id for O(1) access.
  # Deduplicate by physical coordinates first — multiple union exon IDs can
  # map to the same isoform exon coordinates (region_type is consistent).
  dedup_key <- paste(annotated_ue$isoform_id,
                     annotated_ue$isoform_exon_start,
                     annotated_ue$isoform_exon_end)
  aue_dedup <- annotated_ue[!duplicated(dedup_key), , drop = FALSE]
  aue_by_iso <- split(aue_dedup, aue_dedup$isoform_id)

  # CDS lookup for boundary decomposition
  cds_lookup <- cds_metadata[cds_metadata$coding_status == "coding",
                              c("isoform_id", "cds_start", "cds_stop", "strand")]

  all_rows <- list()
  event_counter <- 0L

  for (i in seq_len(nrow(profiles))) {
    events <- profiles$detailed_events[[i]]
    if (is.null(events) || nrow(events) == 0L) next

    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]
    gene_id <- profiles$gene_id[i]

    ref_segs <- aue_by_iso[[ref_id]]
    comp_segs <- aue_by_iso[[comp_id]]
    ref_cds <- cds_lookup[cds_lookup$isoform_id == ref_id, ]
    comp_cds <- cds_lookup[cds_lookup$isoform_id == comp_id, ]

    for (j in seq_len(nrow(events))) {
      event_counter <- event_counter + 1L
      evt <- events[j, ]

      genomic_start <- as.integer(pmin(evt$five_prime, evt$three_prime))
      genomic_end <- as.integer(pmax(evt$five_prime, evt$three_prime))
      event_span <- genomic_end - genomic_start + 1L

      ref_bp <- .decompose_event_bp(genomic_start, genomic_end, ref_segs,
                                     ref_cds, event_span)
      comp_bp <- .decompose_event_bp(genomic_start, genomic_end, comp_segs,
                                      comp_cds, event_span)

      all_rows[[event_counter]] <- tibble::tibble(
        event_id = event_counter,
        gene_id = gene_id,
        reference_isoform_id = ref_id,
        comparator_isoform_id = comp_id,
        event_type = evt$event_type,
        direction = evt$direction,
        genomic_start = genomic_start,
        genomic_end = genomic_end,
        event_span_bp = event_span,
        ref_utr5_bp = ref_bp$utr5, ref_cds_bp = ref_bp$cds,
        ref_utr3_bp = ref_bp$utr3, ref_intronic_bp = ref_bp$intronic,
        comp_utr5_bp = comp_bp$utr5, comp_cds_bp = comp_bp$cds,
        comp_utr3_bp = comp_bp$utr3, comp_intronic_bp = comp_bp$intronic,
        ref_category = ref_bp$category,
        comp_category = comp_bp$category
      )
    }
  }

  if (length(all_rows) == 0L) {
    return(tibble::tibble(
      event_id = integer(0), gene_id = character(0),
      reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      event_type = character(0), direction = character(0),
      genomic_start = integer(0), genomic_end = integer(0),
      event_span_bp = integer(0),
      ref_utr5_bp = integer(0), ref_cds_bp = integer(0),
      ref_utr3_bp = integer(0), ref_intronic_bp = integer(0),
      comp_utr5_bp = integer(0), comp_cds_bp = integer(0),
      comp_utr3_bp = integer(0), comp_intronic_bp = integer(0),
      ref_category = character(0), comp_category = character(0)
    ))
  }

  dplyr::bind_rows(all_rows)
}


#' Quantify Structural Divergence Between Isoform Pairs
#'
#' For each pair in \code{profiles}, computes per-isoform exonic content
#' (total and by 5'UTR, CDS, 3'UTR), inter-isoform differences, and
#' shared/unique exonic basepairs via physical interval overlap of exon sets.
#'
#' @param profiles A tibble from \code{\link{buildProfiles}()}.
#' @param structures A tibble from \code{\link{parseIsoformStructures}()}
#'   with list columns \code{exon_starts} and \code{exon_ends}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}
#'   with \code{isoform_id}, \code{coding_status}, \code{cds_start},
#'   \code{cds_stop}, \code{strand}.
#' @return A tibble with one row per pair. See Details for column descriptions.
#'
#' @details
#' Per-isoform columns: \code{exonic_bp_ref}, \code{exonic_bp_comp} (total
#' exonic bp); \code{utr5_bp_ref}, \code{cds_bp_ref}, \code{utr3_bp_ref},
#' \code{noncoding_bp_ref} (and \code{_comp} variants).
#'
#' Difference columns: \code{exonic_bp_diff}, \code{utr5_bp_diff},
#' \code{cds_bp_diff}, \code{utr3_bp_diff} (all ref minus comp).
#'
#' Overlap columns: \code{shared_exonic_bp}, \code{ref_unique_bp},
#' \code{comp_unique_bp}, \code{union_exonic_bp}, \code{pct_shared}
#' (100 * shared / union).
#'
#' Invariants (asserted):
#' \itemize{
#'   \item \code{utr5 + cds + utr3 + noncoding == exonic_bp} per isoform
#'   \item \code{shared + ref_unique + comp_unique == union}
#' }
#'
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' data(example_profiles)
#' div <- quantifyPairDivergence(example_profiles, structures, cds)
#' head(div)
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
quantifyPairDivergence <- function(profiles, structures, cds_metadata) {

  if (nrow(profiles) == 0L) {
    return(tibble::tibble(
      gene_id = character(0),
      reference_isoform_id = character(0),
      comparator_isoform_id = character(0),
      exonic_bp_ref = integer(0), exonic_bp_comp = integer(0),
      utr5_bp_ref = integer(0), utr5_bp_comp = integer(0),
      cds_bp_ref = integer(0), cds_bp_comp = integer(0),
      utr3_bp_ref = integer(0), utr3_bp_comp = integer(0),
      noncoding_bp_ref = integer(0), noncoding_bp_comp = integer(0),
      exonic_bp_diff = integer(0),
      utr5_bp_diff = integer(0), cds_bp_diff = integer(0),
      utr3_bp_diff = integer(0),
      shared_exonic_bp = integer(0), ref_unique_bp = integer(0),
      comp_unique_bp = integer(0), union_exonic_bp = integer(0),
      pct_shared = numeric(0)
    ))
  }

  # Pre-build lookups for O(1) access
  struct_lookup <- split(structures, structures$isoform_id)
  cds_coding <- cds_metadata[cds_metadata$coding_status == "coding",
                              c("isoform_id", "cds_start", "cds_stop", "strand")]
  cds_lookup <- split(cds_coding, cds_coding$isoform_id)

  all_rows <- vector("list", nrow(profiles))

  for (i in seq_len(nrow(profiles))) {
    ref_id <- profiles$reference_isoform_id[i]
    comp_id <- profiles$comparator_isoform_id[i]
    gene_id <- profiles$gene_id[i]

    # Extract exon coordinates
    ref_struct <- struct_lookup[[ref_id]]
    comp_struct <- struct_lookup[[comp_id]]
    if (is.null(ref_struct) || is.null(comp_struct)) next

    ref_starts <- ref_struct$exon_starts[[1L]]
    ref_ends <- ref_struct$exon_ends[[1L]]
    comp_starts <- comp_struct$exon_starts[[1L]]
    comp_ends <- comp_struct$exon_ends[[1L]]

    # Total exonic bp
    exonic_ref <- sum(ref_ends - ref_starts + 1L)
    exonic_comp <- sum(comp_ends - comp_starts + 1L)

    # Regional bp
    ref_cds_row <- cds_lookup[[ref_id]]
    comp_cds_row <- cds_lookup[[comp_id]]

    ref_cds_start <- if (!is.null(ref_cds_row)) ref_cds_row$cds_start[1L] else NA_integer_
    ref_cds_stop <- if (!is.null(ref_cds_row)) ref_cds_row$cds_stop[1L] else NA_integer_
    ref_strand <- if (!is.null(ref_cds_row)) ref_cds_row$strand[1L] else NA_character_
    comp_cds_start <- if (!is.null(comp_cds_row)) comp_cds_row$cds_start[1L] else NA_integer_
    comp_cds_stop <- if (!is.null(comp_cds_row)) comp_cds_row$cds_stop[1L] else NA_integer_
    comp_strand <- if (!is.null(comp_cds_row)) comp_cds_row$strand[1L] else NA_character_

    ref_reg <- .computeRegionalBp(ref_starts, ref_ends,
                                   ref_cds_start, ref_cds_stop, ref_strand)
    comp_reg <- .computeRegionalBp(comp_starts, comp_ends,
                                    comp_cds_start, comp_cds_stop, comp_strand)

    # Shared/unique exonic bp
    shared <- .computeExonOverlapBp(ref_starts, ref_ends, comp_starts, comp_ends)
    ref_unique <- as.integer(exonic_ref - shared)
    comp_unique <- as.integer(exonic_comp - shared)
    union_bp <- shared + ref_unique + comp_unique
    pct_shared <- if (union_bp > 0L) 100 * shared / union_bp else NA_real_

    all_rows[[i]] <- tibble::tibble(
      gene_id = gene_id,
      reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,
      exonic_bp_ref = as.integer(exonic_ref),
      exonic_bp_comp = as.integer(exonic_comp),
      utr5_bp_ref = ref_reg$utr5, utr5_bp_comp = comp_reg$utr5,
      cds_bp_ref = ref_reg$cds, cds_bp_comp = comp_reg$cds,
      utr3_bp_ref = ref_reg$utr3, utr3_bp_comp = comp_reg$utr3,
      noncoding_bp_ref = ref_reg$noncoding,
      noncoding_bp_comp = comp_reg$noncoding,
      exonic_bp_diff = as.integer(exonic_ref - exonic_comp),
      utr5_bp_diff = as.integer(ref_reg$utr5 - comp_reg$utr5),
      cds_bp_diff = as.integer(ref_reg$cds - comp_reg$cds),
      utr3_bp_diff = as.integer(ref_reg$utr3 - comp_reg$utr3),
      shared_exonic_bp = shared,
      ref_unique_bp = ref_unique,
      comp_unique_bp = comp_unique,
      union_exonic_bp = as.integer(union_bp),
      pct_shared = pct_shared
    )
  }

  result <- dplyr::bind_rows(all_rows)

  # Assert invariants
  if (nrow(result) > 0L) {
    ref_sum <- result$utr5_bp_ref + result$cds_bp_ref +
      result$utr3_bp_ref + result$noncoding_bp_ref
    stopifnot(all(ref_sum == result$exonic_bp_ref))
    comp_sum <- result$utr5_bp_comp + result$cds_bp_comp +
      result$utr3_bp_comp + result$noncoding_bp_comp
    stopifnot(all(comp_sum == result$exonic_bp_comp))
    overlap_sum <- result$shared_exonic_bp + result$ref_unique_bp +
      result$comp_unique_bp
    stopifnot(all(overlap_sum == result$union_exonic_bp))
  }

  result
}


#' Compute regional bp decomposition for one isoform
#'
#' Given exon coordinates and CDS boundaries, decomposes total exonic bp
#' into 5'UTR, CDS, 3'UTR, and non-coding components. Strand-aware.
#'
#' @param exon_starts Integer vector of exon start coordinates.
#' @param exon_ends Integer vector of exon end coordinates.
#' @param cds_start Integer; CDS start (NA for non-coding).
#' @param cds_stop Integer; CDS stop (NA for non-coding).
#' @param strand Character; "+" or "-" (NA for non-coding).
#' @return Named list: \code{utr5}, \code{cds}, \code{utr3}, \code{noncoding}
#'   (all integer).
#' @keywords internal
.computeRegionalBp <- function(exon_starts, exon_ends,
                                cds_start, cds_stop, strand) {
  total_exonic <- sum(exon_ends - exon_starts + 1L)

  # Non-coding isoform
  if (is.na(cds_start) || is.na(cds_stop)) {
    return(list(utr5 = 0L, cds = 0L, utr3 = 0L,
                noncoding = as.integer(total_exonic)))
  }

  utr5 <- 0L
  cds_bp <- 0L
  utr3 <- 0L

  for (k in seq_along(exon_starts)) {
    es <- exon_starts[k]
    ee <- exon_ends[k]

    # Portion before CDS
    if (es < cds_start) {
      pre_end <- min(ee, cds_start - 1L)
      pre_bp <- as.integer(pre_end - es + 1L)
      if (!is.na(strand) && strand == "-") {
        utr3 <- utr3 + pre_bp
      } else {
        utr5 <- utr5 + pre_bp
      }
    }

    # Portion within CDS
    ov_start <- max(es, cds_start)
    ov_end <- min(ee, cds_stop)
    if (ov_start <= ov_end) {
      cds_bp <- cds_bp + as.integer(ov_end - ov_start + 1L)
    }

    # Portion after CDS
    if (ee > cds_stop) {
      post_start <- max(es, cds_stop + 1L)
      post_bp <- as.integer(ee - post_start + 1L)
      if (!is.na(strand) && strand == "-") {
        utr5 <- utr5 + post_bp
      } else {
        utr3 <- utr3 + post_bp
      }
    }
  }

  list(utr5 = utr5, cds = cds_bp, utr3 = utr3, noncoding = 0L)
}


#' Compute exonic basepair overlap between two isoforms
#'
#' Sorted-interval-intersection algorithm. Returns the total number of
#' basepairs where both isoforms have exonic sequence.
#'
#' @param ref_starts Integer vector of reference exon starts.
#' @param ref_ends Integer vector of reference exon ends.
#' @param comp_starts Integer vector of comparator exon starts.
#' @param comp_ends Integer vector of comparator exon ends.
#' @return Integer; total shared exonic bp.
#' @keywords internal
.computeExonOverlapBp <- function(ref_starts, ref_ends,
                                   comp_starts, comp_ends) {
  # Sort both exon sets by start
  ref_ord <- order(ref_starts)
  comp_ord <- order(comp_starts)
  rs <- ref_starts[ref_ord]; re <- ref_ends[ref_ord]
  cs <- comp_starts[comp_ord]; ce <- comp_ends[comp_ord]

  n_ref <- length(rs)
  n_comp <- length(cs)
  shared <- 0L
  ri <- 1L
  ci <- 1L

  while (ri <= n_ref && ci <= n_comp) {
    ov_start <- max(rs[ri], cs[ci])
    ov_end <- min(re[ri], ce[ci])
    if (ov_start <= ov_end) {
      shared <- shared + as.integer(ov_end - ov_start + 1L)
    }
    # Advance the interval that ends first
    if (re[ri] <= ce[ci]) {
      ri <- ri + 1L
    } else {
      ci <- ci + 1L
    }
  }

  shared
}


#' Decompose event bp into regional components for one isoform
#'
#' @param event_start Integer; genomic start of event.
#' @param event_end Integer; genomic end of event.
#' @param iso_segs Data frame of annotated union exon segments for this
#'   isoform (from \code{aue_by_iso}). May be NULL if isoform not in mapping.
#' @param iso_cds One-row data frame with \code{cds_start}, \code{cds_stop},
#'   \code{strand}, or zero-row if non-coding.
#' @param event_span Integer; total event span bp.
#' @return Named list: \code{utr5}, \code{cds}, \code{utr3}, \code{intronic},
#'   \code{category}.
#' @keywords internal
.decompose_event_bp <- function(event_start, event_end, iso_segs, iso_cds,
                                 event_span) {
  utr5 <- 0L; cds_bp <- 0L; utr3 <- 0L; exonic_total <- 0L

  has_cds <- nrow(iso_cds) > 0L
  i_cds_start <- if (has_cds) iso_cds$cds_start[1L] else NA_integer_
  i_cds_stop <- if (has_cds) iso_cds$cds_stop[1L] else NA_integer_
  i_strand <- if (has_cds) iso_cds$strand[1L] else NA_character_

  if (!is.null(iso_segs) && nrow(iso_segs) > 0L) {
    # Find overlapping exon segments
    overlapping <- iso_segs[iso_segs$isoform_exon_start <= event_end &
                            iso_segs$isoform_exon_end >= event_start, ,
                            drop = FALSE]

    for (k in seq_len(nrow(overlapping))) {
      seg <- overlapping[k, ]
      ov_start <- max(seg$isoform_exon_start, event_start)
      ov_end <- min(seg$isoform_exon_end, event_end)
      ov_bp <- as.integer(ov_end - ov_start + 1L)
      if (ov_bp <= 0L) next

      rt <- seg$region_type

      # non_coding and unknown exonic regions are NOT counted as exonic
      # for decomposition purposes — their bp flow to intronic
      if (rt == "non_coding" || rt == "unknown") {
        next
      }

      exonic_total <- exonic_total + ov_bp

      if (rt == "CDS") {
        cds_bp <- cds_bp + ov_bp
      } else if (rt == "5'UTR") {
        utr5 <- utr5 + ov_bp
      } else if (rt == "3'UTR") {
        utr3 <- utr3 + ov_bp
      } else if (has_cds && grepl("contains_orf", rt)) {
        # Decompose boundary exon using CDS interval [cds_start, cds_stop]
        # This is strand-independent because cds_start <= cds_stop always
        # and the region_type labels tell us which side is 5'UTR vs 3'UTR
        cds_ov_start <- max(ov_start, i_cds_start)
        cds_ov_end <- min(ov_end, i_cds_stop)
        cds_portion <- max(0L, as.integer(cds_ov_end - cds_ov_start + 1L))
        cds_bp <- cds_bp + cds_portion

        # Portion before CDS interval
        if (ov_start < i_cds_start) {
          pre_bp <- as.integer(min(ov_end, i_cds_start - 1L) - ov_start + 1L)
          if (pre_bp > 0L) {
            # On + strand: before CDS = 5'UTR; on - strand: before CDS = 3'UTR
            if (!is.na(i_strand) && i_strand == "-") {
              utr3 <- utr3 + pre_bp
            } else {
              utr5 <- utr5 + pre_bp
            }
          }
        }
        # Portion after CDS interval
        if (ov_end > i_cds_stop) {
          post_bp <- as.integer(ov_end - max(ov_start, i_cds_stop + 1L) + 1L)
          if (post_bp > 0L) {
            if (!is.na(i_strand) && i_strand == "-") {
              utr5 <- utr5 + post_bp
            } else {
              utr3 <- utr3 + post_bp
            }
          }
        }
      }
    }
  }

  intronic <- as.integer(event_span - exonic_total)
  if (intronic < 0L) intronic <- 0L

  # Build category string from non-zero regions
  parts <- character(0)
  if (utr5 > 0L) parts <- c(parts, "5'UTR")
  if (cds_bp > 0L) parts <- c(parts, "CDS")
  if (utr3 > 0L) parts <- c(parts, "3'UTR")
  if (intronic > 0L) parts <- c(parts, "intronic")
  category <- if (length(parts) == 0L) "non_coding" else paste(parts, collapse = " + ")

  list(utr5 = utr5, cds = cds_bp, utr3 = utr3, intronic = intronic,
       category = category)
}


#' Bootstrap Regional Comparison Between Two Profile Sets
#'
#' For each event-type x region-type cell, tests whether the proportion
#' differs between two sets of mapped event regions using permutation
#' bootstrap. Provides cell-level enrichment detail beyond the chi-square
#' in [comparePairSets()].
#'
#' @param event_regions_a A tibble from [mapEventsToRegions()] for set A.
#' @param event_regions_b A tibble from [mapEventsToRegions()] for set B.
#' @param n_boot Integer; number of bootstrap iterations (default 1000).
#' @param min_events Integer; minimum events per event_type in **both** sets
#'   to run bootstrap. Below this threshold, the cell is skipped with
#'   \code{p_value = NA} (default 5).
#' @param seed Optional integer for reproducibility. Sets the global RNG
#'   state via \code{set.seed()}.
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A tibble with columns: \code{event_type}, \code{region_type},
#'   \code{n_events_a}, \code{n_events_b}, \code{prop_a}, \code{prop_b},
#'   \code{diff}, \code{ci_lower}, \code{ci_upper}, \code{p_value},
#'   \code{adj_p_value}, \code{significant}, \code{skipped}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' er <- mapEventsToRegions(example_profiles, annotated, cds)
#' # Self-comparison (expect no significance)
#' boot <- bootstrapRegionalComparison(er, er, n_boot = 100,
#'   seed = 42, verbose = FALSE)
#' head(boot)
#' @export
#' @importFrom stats quantile p.adjust
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
bootstrapRegionalComparison <- function(event_regions_a, event_regions_b,
                                        n_boot = 1000L, min_events = 5L,
                                        seed = NULL, verbose = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  # Empty input handling
  if (nrow(event_regions_a) == 0L && nrow(event_regions_b) == 0L) {
    return(tibble::tibble(
      event_type = character(0), region_type = character(0),
      n_events_a = integer(0), n_events_b = integer(0),
      prop_a = numeric(0), prop_b = numeric(0), diff = numeric(0),
      ci_lower = numeric(0), ci_upper = numeric(0),
      p_value = numeric(0), adj_p_value = numeric(0),
      significant = logical(0), skipped = logical(0)
    ))
  }

  # Deduplicate events to one row per event_row_id × region_type
  dedup_a <- unique(event_regions_a[, c("event_row_id", "event_type",
                                         "region_type")])
  dedup_b <- unique(event_regions_b[, c("event_row_id", "event_type",
                                         "region_type")])

  # Count events per event_type in each set (distinct event_row_id)
  count_events <- function(dedup) {
    tapply(dedup$event_row_id, dedup$event_type,
           function(x) length(unique(x)))
  }
  n_a_by_et <- count_events(dedup_a)
  n_b_by_et <- count_events(dedup_b)

  all_event_types <- union(names(n_a_by_et), names(n_b_by_et))
  all_results <- list()

  for (et in all_event_types) {
    n_a <- if (et %in% names(n_a_by_et)) n_a_by_et[[et]] else 0L
    n_b <- if (et %in% names(n_b_by_et)) n_b_by_et[[et]] else 0L

    # Get region-level data for this event type
    et_a <- dedup_a[dedup_a$event_type == et, ]
    et_b <- dedup_b[dedup_b$event_type == et, ]
    all_regions <- union(et_a$region_type, et_b$region_type)

    if (n_a < min_events || n_b < min_events) {
      # Skip: too few events
      for (rt in all_regions) {
        all_results[[length(all_results) + 1L]] <- tibble::tibble(
          event_type = et, region_type = rt,
          n_events_a = as.integer(n_a), n_events_b = as.integer(n_b),
          prop_a = NA_real_, prop_b = NA_real_, diff = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_,
          p_value = NA_real_, adj_p_value = NA_real_,
          significant = NA, skipped = TRUE
        )
      }
      if (length(all_regions) == 0L) {
        all_results[[length(all_results) + 1L]] <- tibble::tibble(
          event_type = et, region_type = NA_character_,
          n_events_a = as.integer(n_a), n_events_b = as.integer(n_b),
          prop_a = NA_real_, prop_b = NA_real_, diff = NA_real_,
          ci_lower = NA_real_, ci_upper = NA_real_,
          p_value = NA_real_, adj_p_value = NA_real_,
          significant = NA, skipped = TRUE
        )
      }
      next
    }

    # Count per-region proportions
    count_region <- function(et_data, n_total) {
      vapply(all_regions, function(rt) {
        sum(et_data$region_type == rt) / n_total
      }, numeric(1))
    }

    obs_prop_a <- count_region(et_a, n_a)
    obs_prop_b <- count_region(et_b, n_b)
    obs_diff <- obs_prop_a - obs_prop_b

    # Pool events, permute labels.
    # Prefix event_row_id with source to avoid collisions when both sets
    # start counting from 1 (separate mapEventsToRegions calls).
    pool_a <- et_a[, c("event_row_id", "region_type")]
    pool_b <- et_b[, c("event_row_id", "region_type")]
    pool_a$pool_id <- paste0("a_", pool_a$event_row_id)
    pool_b$pool_id <- paste0("b_", pool_b$event_row_id)
    pooled <- dplyr::bind_rows(pool_a, pool_b)
    # Unique event IDs for permutation (using disambiguated pool_id)
    unique_ids_a <- unique(pool_a$pool_id)
    unique_ids_b <- unique(pool_b$pool_id)
    all_unique_ids <- c(unique_ids_a, unique_ids_b)
    n_total <- length(all_unique_ids)

    boot_diffs <- matrix(NA_real_, nrow = n_boot, ncol = length(all_regions))

    for (b in seq_len(n_boot)) {
      perm <- sample(n_total)
      perm_a_ids <- all_unique_ids[perm[seq_len(n_a)]]
      perm_b_ids <- all_unique_ids[perm[(n_a + 1L):n_total]]

      perm_a <- pooled[pooled$pool_id %in% perm_a_ids, ]
      perm_b <- pooled[pooled$pool_id %in% perm_b_ids, ]

      perm_prop_a <- vapply(all_regions, function(rt) {
        sum(perm_a$region_type == rt) / n_a
      }, numeric(1))
      perm_prop_b <- vapply(all_regions, function(rt) {
        sum(perm_b$region_type == rt) / n_b
      }, numeric(1))

      boot_diffs[b, ] <- perm_prop_a - perm_prop_b
    }

    for (r_idx in seq_along(all_regions)) {
      rt <- all_regions[r_idx]
      boot_col <- boot_diffs[, r_idx]
      ci <- stats::quantile(boot_col, probs = c(0.025, 0.975), na.rm = TRUE)
      p_val <- (sum(abs(boot_col) >= abs(obs_diff[r_idx])) + 1L) /
        (n_boot + 1L)

      all_results[[length(all_results) + 1L]] <- tibble::tibble(
        event_type = et, region_type = rt,
        n_events_a = as.integer(n_a), n_events_b = as.integer(n_b),
        prop_a = obs_prop_a[r_idx], prop_b = obs_prop_b[r_idx],
        diff = obs_diff[r_idx],
        ci_lower = as.numeric(ci[1L]),
        ci_upper = as.numeric(ci[2L]),
        p_value = p_val, adj_p_value = NA_real_,
        significant = NA, skipped = FALSE
      )
    }

    if (verbose) message(sprintf("  Bootstrap: %s (%d regions)", et,
                                  length(all_regions)))
  }

  result <- dplyr::bind_rows(all_results)

  # FDR correction across non-skipped cells
  testable <- !result$skipped
  if (any(testable)) {
    result$adj_p_value[testable] <- stats::p.adjust(
      result$p_value[testable], method = "fdr"
    )
    result$significant[testable] <- result$adj_p_value[testable] < 0.05
  }

  result
}


#' Test Regional Enrichment of Events
#'
#' Tests whether events are enriched or depleted in specific genomic regions
#' relative to the expected proportion based on region bp coverage. Uses
#' binomial tests with FDR correction.
#'
#' @param event_regions A tibble from \code{\link{mapEventsToRegions}()}.
#' @param annotated_ue A tibble from \code{\link{annotateRegionTypes}()}.
#' @param profiles A tibble from \code{\link{buildProfiles}()} (used to
#'   extract reference isoform IDs for null distribution).
#' @return A tibble with columns: \code{event_type}, \code{region_type},
#'   \code{observed_count}, \code{total_events}, \code{observed_prop},
#'   \code{expected_prop}, \code{enrichment_ratio}, \code{p_value},
#'   \code{adj_p_value}, \code{significant}, \code{direction_enrich}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' er <- mapEventsToRegions(example_profiles, annotated, cds)
#' enrich <- testRegionalEnrichment(er, annotated, example_profiles)
#' enrich
#' @export
#' @importFrom stats binom.test p.adjust
#' @importFrom dplyr distinct
#' @importFrom tibble tibble
testRegionalEnrichment <- function(event_regions, annotated_ue, profiles) {

  if (nrow(event_regions) == 0L) {
    return(tibble::tibble(
      event_type = character(0), region_type = character(0),
      observed_count = integer(0), total_events = integer(0),
      observed_prop = numeric(0), expected_prop = numeric(0),
      enrichment_ratio = numeric(0), p_value = numeric(0),
      adj_p_value = numeric(0), significant = logical(0),
      direction_enrich = character(0)
    ))
  }

  # Compute expected proportions: bp per region in reference isoforms
  ref_ids <- unique(profiles$reference_isoform_id)
  ref_ue <- annotated_ue[annotated_ue$isoform_id %in% ref_ids, ]

  # Deduplicate by isoform + union_exon to avoid double-counting
  ref_ue_dedup <- dplyr::distinct(
    ref_ue, .data$isoform_id, .data$union_exon_id, .keep_all = TRUE
  )

  ref_ue_dedup$bp <- ref_ue_dedup$isoform_exon_end -
    ref_ue_dedup$isoform_exon_start + 1L

  # Expected proportions by region type
  region_bp <- tapply(ref_ue_dedup$bp, ref_ue_dedup$region_type, sum)
  total_bp <- sum(region_bp)
  expected_props <- region_bp / total_bp

  # Observed counts per event_type x region_type
  event_types <- unique(event_regions$event_type)
  all_results <- list()

  for (et in event_types) {
    et_data <- event_regions[event_regions$event_type == et, ]
    total_ev <- length(unique(et_data$event_row_id))

    # Count distinct events per region type
    regions <- unique(et_data$region_type)
    for (rt in regions) {
      n_obs <- length(unique(
        et_data$event_row_id[et_data$region_type == rt]
      ))
      obs_prop <- n_obs / total_ev
      exp_prop <- if (rt %in% names(expected_props)) {
        as.numeric(expected_props[rt])
      } else {
        0
      }

      bt <- tryCatch(
        stats::binom.test(n_obs, total_ev, p = max(exp_prop, 1e-6)),
        error = function(e) NULL
      )

      all_results[[length(all_results) + 1L]] <- tibble::tibble(
        event_type = et,
        region_type = rt,
        observed_count = n_obs,
        total_events = total_ev,
        observed_prop = obs_prop,
        expected_prop = exp_prop,
        enrichment_ratio = if (exp_prop > 0) obs_prop / exp_prop else NA_real_,
        p_value = if (!is.null(bt)) bt$p.value else NA_real_,
        adj_p_value = NA_real_,
        significant = NA,
        direction_enrich = NA_character_
      )
    }
  }

  result <- dplyr::bind_rows(all_results)
  result$adj_p_value <- stats::p.adjust(result$p_value, method = "fdr")
  result$significant <- result$adj_p_value < 0.05
  result$direction_enrich <- ifelse(
    !result$significant, "NS",
    ifelse(result$enrichment_ratio > 1, "Enriched", "Depleted")
  )

  result
}


#' Test ORF Boundary Susceptibility
#'
#' Tests whether A5SS/A3SS events are under- or over-represented at ORF
#' boundary exons (contains_orf_start, contains_orf_stop) compared to
#' regular CDS exons, using Fisher's exact test.
#'
#' @param event_regions A tibble from \code{\link{mapEventsToRegions}()}.
#' @param annotated_ue A tibble from \code{\link{annotateRegionTypes}()}.
#' @return A tibble with columns: \code{boundary_type}, \code{n_exons},
#'   \code{n_with_events}, \code{rate}, \code{odds_ratio}, \code{p_value}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' er <- mapEventsToRegions(example_profiles, annotated, cds)
#' orf <- testOrfBoundarySusceptibility(er, annotated)
#' orf
#' @export
#' @importFrom stats fisher.test
#' @importFrom tibble tibble
testOrfBoundarySusceptibility <- function(event_regions, annotated_ue) {

  # Filter to A5SS/A3SS events
  a5a3 <- event_regions[event_regions$event_type %in% c("A5SS", "A3SS"), ]

  # Get unique reference isoform IDs from event_regions
  ref_ids <- unique(event_regions$reference_isoform_id)
  ref_ue <- annotated_ue[annotated_ue$isoform_id %in% ref_ids, ]
  ref_ue <- dplyr::distinct(
    ref_ue, .data$isoform_id, .data$union_exon_id, .keep_all = TRUE
  )

  # Count exons by type
  orf_boundary_types <- c("contains_orf_start", "contains_orf_stop",
                           "contains_orf_start_stop")
  n_orf_boundary <- sum(ref_ue$region_type %in% orf_boundary_types)
  n_regular_cds <- sum(ref_ue$region_type == "CDS")

  if (n_orf_boundary == 0L || n_regular_cds == 0L) {
    return(tibble::tibble(
      boundary_type = c("ORF_boundary", "Regular_CDS"),
      n_exons = c(n_orf_boundary, n_regular_cds),
      n_with_events = c(0L, 0L),
      rate = c(NA_real_, NA_real_),
      odds_ratio = NA_real_,
      p_value = NA_real_
    ))
  }

  # Count unique exons (isoform × union_exon_id) with at least one A5SS/A3SS
  # overlapping event. Both dimensions of the contingency table count exons,
  # making the Fisher's test formally valid.
  hit_exon_keys <- character(0)
  for (k in seq_len(nrow(a5a3))) {
    overlaps <- ref_ue[
      ref_ue$isoform_id == a5a3$reference_isoform_id[k] &
      ref_ue$isoform_exon_start <= a5a3$genomic_end[k] &
      ref_ue$isoform_exon_end >= a5a3$genomic_start[k], ]
    if (nrow(overlaps) > 0L) {
      hit_exon_keys <- c(hit_exon_keys,
        paste(overlaps$isoform_id, overlaps$union_exon_id, sep = "|"))
    }
  }
  hit_exon_keys <- unique(hit_exon_keys)

  # Build exon keys for each class
  all_exon_keys <- paste(ref_ue$isoform_id, ref_ue$union_exon_id, sep = "|")
  orf_exon_keys <- unique(all_exon_keys[ref_ue$region_type %in% orf_boundary_types])
  cds_exon_keys <- unique(all_exon_keys[ref_ue$region_type == "CDS"])

  orf_hit <- sum(orf_exon_keys %in% hit_exon_keys)
  cds_hit <- sum(cds_exon_keys %in% hit_exon_keys)

  # Fisher's exact test: exons with/without A5SS/A3SS events
  contingency <- matrix(
    c(orf_hit, n_orf_boundary - orf_hit,
      cds_hit, n_regular_cds - cds_hit),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("ORF_boundary", "Regular_CDS"),
                    c("Has_event", "No_event"))
  )

  ft <- tryCatch(stats::fisher.test(contingency), error = function(e) NULL)

  tibble::tibble(
    boundary_type = c("ORF_boundary", "Regular_CDS"),
    n_exons = c(n_orf_boundary, n_regular_cds),
    n_with_events = c(orf_hit, cds_hit),
    rate = c(orf_hit / n_orf_boundary, cds_hit / n_regular_cds),
    odds_ratio = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
    p_value = if (!is.null(ft)) ft$p.value else NA_real_
  )
}


#' Summarize ORF Impact of Terminal Events
#'
#' For Alt_TSS and Alt_TES events, determines what fraction overlap the ORF
#' start or stop codon in the reference isoform.
#'
#' @param event_regions A tibble from \code{\link{mapEventsToRegions}()}.
#' @param cds_metadata A tibble from \code{\link{extractCdsAnnotations}()}.
#' @return A tibble with columns: \code{event_type}, \code{n_events},
#'   \code{n_affecting_orf_start}, \code{pct_affecting_orf_start},
#'   \code{n_affecting_orf_stop}, \code{pct_affecting_orf_stop}.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
#'   verbose = FALSE)
#' annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
#' data(example_profiles)
#' er <- mapEventsToRegions(example_profiles, annotated, cds)
#' orf_impact <- summarizeOrfImpact(er, cds)
#' orf_impact
#' @export
#' @importFrom tibble tibble
summarizeOrfImpact <- function(event_regions, cds_metadata) {

  terminal_types <- c("Alt_TSS", "Alt_TES")
  results <- list()

  for (et in terminal_types) {
    et_events <- event_regions[event_regions$event_type == et, ]
    total <- length(unique(et_events$event_row_id))

    if (total == 0L) {
      results[[length(results) + 1L]] <- tibble::tibble(
        event_type = et,
        n_events = 0L,
        n_affecting_orf_start = 0L,
        pct_affecting_orf_start = NA_real_,
        n_affecting_orf_stop = 0L,
        pct_affecting_orf_stop = NA_real_
      )
      next
    }

    # Events overlapping ORF start region
    orf_start_regions <- c("contains_orf_start", "contains_orf_start_stop")
    n_orf_start <- length(unique(
      et_events$event_row_id[et_events$region_type %in% orf_start_regions]
    ))

    # Events overlapping ORF stop region
    orf_stop_regions <- c("contains_orf_stop", "contains_orf_start_stop")
    n_orf_stop <- length(unique(
      et_events$event_row_id[et_events$region_type %in% orf_stop_regions]
    ))

    results[[length(results) + 1L]] <- tibble::tibble(
      event_type = et,
      n_events = total,
      n_affecting_orf_start = n_orf_start,
      pct_affecting_orf_start = 100 * n_orf_start / total,
      n_affecting_orf_stop = n_orf_stop,
      pct_affecting_orf_stop = 100 * n_orf_stop / total
    )
  }

  dplyr::bind_rows(results)
}
