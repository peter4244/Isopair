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

  # Count A5SS/A3SS events hitting each exon class
  # An event hits an ORF boundary exon if its region_type includes boundary types
  orf_hit <- length(unique(
    a5a3$event_row_id[a5a3$region_type %in% orf_boundary_types]
  ))
  cds_hit <- length(unique(
    a5a3$event_row_id[a5a3$region_type == "CDS"]
  ))

  # Fisher's exact test
  contingency <- matrix(
    c(orf_hit, n_orf_boundary - orf_hit,
      cds_hit, n_regular_cds - cds_hit),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("ORF_boundary", "Regular_CDS"),
                    c("Has_event", "No_event"))
  )

  # Ensure no negative values
  contingency <- pmax(contingency, 0L)

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
