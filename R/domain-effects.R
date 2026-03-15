#' Classify Domain Effects of Splice Events
#'
#' For each splice event mapped to protein coordinates by
#' \code{\link{mapSpliceToProtein}()}, checks whether any protein domains
#' overlap the affected region and classifies the effect as \code{disrupted},
#' \code{lost}, \code{gained_intact}, or \code{gained_partial}.
#'
#' Domain lookup is direction-aware:
#' \itemize{
#'   \item For \strong{LOSS} events, domains are looked up on the
#'     \strong{reference} isoform (the comparator lost this region).
#'   \item For \strong{GAIN} events, domains are looked up on the
#'     \strong{comparator} isoform (the comparator gained this region).
#' }
#'
#' Effect classification:
#' \itemize{
#'   \item \code{lost}: domain entirely within a LOSS event region.
#'   \item \code{disrupted}: domain partially overlaps a LOSS event
#'     (spans event boundary).
#'   \item \code{gained_intact}: domain entirely within a GAIN event region.
#'   \item \code{gained_partial}: domain partially overlaps a GAIN event.
#' }
#'
#' @param protein_events A tibble from \code{\link{mapSpliceToProtein}()} with
#'   columns including \code{reference_isoform_id},
#'   \code{comparator_isoform_id}, \code{direction},
#'   \code{protein_start_aa}, \code{protein_end_aa}.
#' @param domain_annotations A data.frame or tibble with at minimum:
#'   \code{isoform_id}, \code{domain_id}, \code{domain_name},
#'   \code{domain_start} (AA, 1-based), \code{domain_end} (AA, 1-based).
#'   Additional columns are preserved but not required.
#' @return A tibble with one row per domain-event overlap. Contains all
#'   columns from \code{protein_events} for the matched event, plus
#'   \code{domain_id}, \code{domain_name}, \code{domain_start},
#'   \code{domain_end}, \code{isoform_source} (\code{"reference"} or
#'   \code{"comparator"}), \code{effect}, \code{overlap_start},
#'   \code{overlap_end}. Returns a zero-row tibble with correct columns if
#'   no overlaps are found.
#' @examples
#' protein_events <- tibble::tibble(
#'   gene_id = "gene1",
#'   reference_isoform_id = "tx1",
#'   comparator_isoform_id = "tx2",
#'   event_type = "SE",
#'   direction = "LOSS",
#'   genomic_start = 350L,
#'   genomic_end = 450L,
#'   cds_nt_start = 100L,
#'   cds_nt_end = 200L,
#'   protein_start_aa = 34L,
#'   protein_end_aa = 67L,
#'   protein_start_aa_buffered = 29L,
#'   protein_end_aa_buffered = 72L,
#'   event_affects_frame = FALSE
#' )
#' domains <- tibble::tibble(
#'   isoform_id = "tx1",
#'   domain_id = "PF00091",
#'   domain_name = "Tubulin",
#'   domain_start = 40L,
#'   domain_end = 60L
#' )
#' result <- classifyDomainEffects(protein_events, domains)
#' result
#' @export
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
classifyDomainEffects <- function(protein_events, domain_annotations) {

  # --- Input validation ---
  if (!is.data.frame(protein_events)) {
    stop("'protein_events' must be a data frame or tibble.", call. = FALSE)
  }
  if (!is.data.frame(domain_annotations)) {
    stop("'domain_annotations' must be a data frame or tibble.", call. = FALSE)
  }

  required_pe <- c("reference_isoform_id", "comparator_isoform_id",
                   "direction", "protein_start_aa", "protein_end_aa")
  missing_pe <- setdiff(required_pe, names(protein_events))
  if (length(missing_pe) > 0L) {
    stop(sprintf("'protein_events' is missing required columns: %s",
                 paste(missing_pe, collapse = ", ")), call. = FALSE)
  }

  required_da <- c("isoform_id", "domain_id", "domain_name",
                   "domain_start", "domain_end")
  missing_da <- setdiff(required_da, names(domain_annotations))
  if (length(missing_da) > 0L) {
    stop(sprintf("'domain_annotations' is missing required columns: %s",
                 paste(missing_da, collapse = ", ")), call. = FALSE)
  }

  # --- Empty output template ---
  # All protein_events columns + domain columns
  pe_cols <- names(protein_events)
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
    event_affects_frame = logical(0),
    domain_id = character(0),
    domain_name = character(0),
    domain_start = integer(0),
    domain_end = integer(0),
    isoform_source = character(0),
    effect = character(0),
    overlap_start = integer(0),
    overlap_end = integer(0)
  )

  if (nrow(protein_events) == 0L || nrow(domain_annotations) == 0L) {
    return(empty_result)
  }

  # Remove domain rows with NA positions
  domain_annotations <- domain_annotations[
    !is.na(domain_annotations$domain_start) &
    !is.na(domain_annotations$domain_end), , drop = FALSE]

  if (nrow(domain_annotations) == 0L) return(empty_result)

  # --- Build domain lookup by isoform ---
  domain_lookup <- split(domain_annotations, domain_annotations$isoform_id)

  all_results <- list()

  for (i in seq_len(nrow(protein_events))) {
    evt <- protein_events[i, ]

    # Skip events with NA protein coordinates
    if (is.na(evt$protein_start_aa) || is.na(evt$protein_end_aa)) next

    evt_start <- evt$protein_start_aa
    evt_end <- evt$protein_end_aa

    # Determine which isoform's domains to check
    if (evt$direction == "LOSS") {
      lookup_id <- evt$reference_isoform_id
      isoform_source <- "reference"
    } else if (evt$direction == "GAIN") {
      lookup_id <- evt$comparator_isoform_id
      isoform_source <- "comparator"
    } else {
      next
    }

    domains <- domain_lookup[[lookup_id]]
    if (is.null(domains) || nrow(domains) == 0L) next

    for (j in seq_len(nrow(domains))) {
      dom <- domains[j, ]
      dom_start <- dom$domain_start
      dom_end <- dom$domain_end

      # Check overlap
      ov_start <- max(evt_start, dom_start)
      ov_end <- min(evt_end, dom_end)

      if (ov_start > ov_end) next  # no overlap

      # Classify effect
      domain_within_event <- dom_start >= evt_start && dom_end <= evt_end

      if (evt$direction == "LOSS") {
        if (domain_within_event) {
          effect <- "lost"
        } else {
          effect <- "disrupted"
        }
      } else {
        # GAIN
        if (domain_within_event) {
          effect <- "gained_intact"
        } else {
          effect <- "gained_partial"
        }
      }

      # Build result row: all protein_events columns + domain columns
      row <- evt
      row$domain_id <- dom$domain_id
      row$domain_name <- dom$domain_name
      row$domain_start <- as.integer(dom_start)
      row$domain_end <- as.integer(dom_end)
      row$isoform_source <- isoform_source
      row$effect <- effect
      row$overlap_start <- as.integer(ov_start)
      row$overlap_end <- as.integer(ov_end)

      all_results[[length(all_results) + 1L]] <- row
    }
  }

  if (length(all_results) == 0L) return(empty_result)
  dplyr::bind_rows(all_results)
}


#' Test Domain Enrichment in Disrupted/Lost Events
#'
#' Tests whether specific protein domains are disrupted or lost more often
#' than expected given their prevalence among reference isoforms. Uses
#' Fisher's exact test per domain with multiple-testing correction.
#'
#' @param domain_effects A tibble from \code{\link{classifyDomainEffects}()},
#'   or at minimum a data.frame with columns \code{domain_id},
#'   \code{domain_name}, \code{reference_isoform_id},
#'   \code{comparator_isoform_id}, \code{effect}.
#' @param domain_annotations A data.frame or tibble with at minimum
#'   \code{isoform_id}, \code{domain_id}, \code{domain_name}. Used to
#'   compute reference prevalence (how many reference isoforms carry each
#'   domain).
#' @param reference_isoform_ids Character vector of all reference isoform IDs
#'   in the analysis (denominator for prevalence calculation).
#' @param n_pairs Integer; total number of isoform pairs tested (denominator
#'   for disruption rate).
#' @param p_adjust_method Character; method for \code{\link[stats]{p.adjust}}
#'   (default \code{"BH"}).
#' @param effects_to_test Character vector of effect types to count as
#'   "affected" (default \code{c("disrupted", "lost")}).
#' @return A tibble with one row per tested domain, sorted by \code{fisher_p}
#'   ascending. Columns: \code{domain_id}, \code{domain_name},
#'   \code{n_ref_with_domain}, \code{n_affected}, \code{n_pairs_total},
#'   \code{fisher_OR}, \code{fisher_p}, \code{padj}, \code{significant}.
#' @examples
#' domain_effects <- tibble::tibble(
#'   domain_id = c("PF00091", "PF00091"),
#'   domain_name = c("Tubulin", "Tubulin"),
#'   reference_isoform_id = c("tx1", "tx2"),
#'   comparator_isoform_id = c("tx1c", "tx2c"),
#'   effect = c("lost", "disrupted")
#' )
#' domain_annotations <- tibble::tibble(
#'   isoform_id = c("tx1", "tx2", "tx3"),
#'   domain_id = c("PF00091", "PF00091", "PF00091"),
#'   domain_name = c("Tubulin", "Tubulin", "Tubulin"),
#'   domain_start = c(10L, 10L, 10L),
#'   domain_end = c(50L, 50L, 50L)
#' )
#' ref_ids <- c("tx1", "tx2", "tx3", "tx4", "tx5")
#' result <- testDomainEnrichment(domain_effects, domain_annotations,
#'   ref_ids, n_pairs = 5L)
#' result
#' @export
#' @importFrom stats fisher.test p.adjust
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
testDomainEnrichment <- function(domain_effects, domain_annotations,
                                  reference_isoform_ids, n_pairs,
                                  p_adjust_method = "BH",
                                  effects_to_test = c("disrupted", "lost")) {

  # --- Input validation ---
  if (!is.data.frame(domain_effects)) {
    stop("'domain_effects' must be a data frame or tibble.", call. = FALSE)
  }
  if (!is.data.frame(domain_annotations)) {
    stop("'domain_annotations' must be a data frame or tibble.", call. = FALSE)
  }
  if (!is.character(reference_isoform_ids) || length(reference_isoform_ids) == 0L) {
    stop("'reference_isoform_ids' must be a non-empty character vector.",
         call. = FALSE)
  }

  n_pairs <- as.integer(n_pairs)
  if (is.na(n_pairs) || n_pairs < 1L) {
    stop("'n_pairs' must be a positive integer.", call. = FALSE)
  }

  required_de <- c("domain_id", "reference_isoform_id",
                   "comparator_isoform_id", "effect")
  missing_de <- setdiff(required_de, names(domain_effects))
  if (length(missing_de) > 0L) {
    stop(sprintf("'domain_effects' is missing required columns: %s",
                 paste(missing_de, collapse = ", ")), call. = FALSE)
  }

  required_da <- c("isoform_id", "domain_id", "domain_name")
  missing_da <- setdiff(required_da, names(domain_annotations))
  if (length(missing_da) > 0L) {
    stop(sprintf("'domain_annotations' is missing required columns: %s",
                 paste(missing_da, collapse = ", ")), call. = FALSE)
  }

  # --- Empty output template ---
  empty_result <- tibble::tibble(
    domain_id = character(0),
    domain_name = character(0),
    n_ref_with_domain = integer(0),
    n_affected = integer(0),
    n_pairs_total = integer(0),
    fisher_OR = numeric(0),
    fisher_p = numeric(0),
    padj = numeric(0),
    significant = logical(0)
  )

  if (nrow(domain_effects) == 0L) return(empty_result)

  # --- 1. Reference prevalence: how many reference isoforms have each domain ---
  ref_domains <- domain_annotations[
    domain_annotations$isoform_id %in% reference_isoform_ids, , drop = FALSE]

  if (nrow(ref_domains) == 0L) return(empty_result)

  # Unique domain_id -> domain_name mapping (take first name per domain)
  domain_names <- ref_domains[!duplicated(ref_domains$domain_id),
                               c("domain_id", "domain_name"), drop = FALSE]

  # Count unique reference isoforms per domain
  ref_prevalence <- tapply(
    ref_domains$isoform_id, ref_domains$domain_id,
    function(x) length(unique(x))
  )

  # --- 2. Count affected pairs per domain ---
  affected <- domain_effects[domain_effects$effect %in% effects_to_test,
                              , drop = FALSE]

  if (nrow(affected) == 0L) {
    # All domains tested but none affected
    all_results <- list()
    for (k in seq_along(ref_prevalence)) {
      did <- names(ref_prevalence)[k]
      dname <- domain_names$domain_name[domain_names$domain_id == did][1L]
      n_ref <- as.integer(ref_prevalence[k])
      n_ref_without <- as.integer(length(reference_isoform_ids) - n_ref)

      contingency <- matrix(
        c(0L, n_pairs,
          n_ref, as.integer(length(reference_isoform_ids)) - n_ref),
        nrow = 2, byrow = TRUE
      )

      ft <- tryCatch(
        stats::fisher.test(contingency),
        error = function(e) NULL
      )

      all_results[[length(all_results) + 1L]] <- tibble::tibble(
        domain_id = did,
        domain_name = if (!is.na(dname)) dname else NA_character_,
        n_ref_with_domain = n_ref,
        n_affected = 0L,
        n_pairs_total = n_pairs,
        fisher_OR = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
        fisher_p = if (!is.null(ft)) ft$p.value else NA_real_,
        padj = NA_real_,
        significant = NA
      )
    }
    result <- dplyr::bind_rows(all_results)
    if (nrow(result) > 0L) {
      result$padj <- stats::p.adjust(result$fisher_p, method = p_adjust_method)
      result$significant <- result$padj < 0.05
      result <- result[order(result$fisher_p), , drop = FALSE]
    }
    return(result)
  }

  # Create unique pair keys for counting
  affected$pair_key <- paste(affected$reference_isoform_id,
                              affected$comparator_isoform_id, sep = "||")

  affected_counts <- tapply(
    affected$pair_key, affected$domain_id,
    function(x) length(unique(x))
  )

  # --- 3. Fisher's exact test per domain ---
  n_total_refs <- length(reference_isoform_ids)
  all_domain_ids <- unique(c(names(ref_prevalence), names(affected_counts)))

  all_results <- list()

  for (did in all_domain_ids) {
    n_ref <- if (did %in% names(ref_prevalence)) {
      as.integer(ref_prevalence[did])
    } else {
      0L
    }
    n_aff <- if (did %in% names(affected_counts)) {
      as.integer(affected_counts[did])
    } else {
      0L
    }
    dname <- domain_names$domain_name[domain_names$domain_id == did]
    dname <- if (length(dname) > 0L) dname[1L] else NA_character_

    # 2x2 table:
    #                  Has domain in ref   No domain in ref
    # Affected              a                    b
    # Not affected          c                    d
    a <- n_aff
    b <- n_pairs - n_aff  # affected pairs without this domain (approx)
    # But more precisely: we want to test whether pairs whose reference
    # has this domain are more likely to be affected.
    # Affected among domain-bearing refs vs affected among non-domain refs.
    # However, we don't have per-pair domain status in the simple case.
    # The spec says: 2x2 of (affected/not) x (has domain in ref/doesn't)
    # So rows = affected/not, cols = has domain/doesn't

    n_without <- as.integer(n_total_refs - n_ref)

    # Pairs where ref has domain and domain was affected
    # vs pairs where ref has domain but domain was not affected
    # vs pairs where ref doesn't have domain (cannot be affected for this domain)
    # This is a classic enrichment test:
    #                  Has domain    No domain
    # Affected            a           n_aff - a  (but affected for other reasons)
    #
    # Actually, reading the spec more carefully:
    # rows: affected/not affected pairs
    # cols: ref has domain / ref doesn't have domain
    # a = pairs where ref has domain AND domain was affected
    # b = pairs where ref does NOT have domain AND pair was somehow affected
    #     (this doesn't apply to this specific domain)
    # The natural 2x2 is:
    #                  Has domain in ref    No domain in ref   Total
    # Affected              n_aff                 0            n_aff
    # Not affected      n_ref - n_aff      n_pairs - n_ref    n_pairs - n_aff
    # Total                n_ref            n_pairs - n_ref    n_pairs

    contingency <- matrix(
      c(n_aff,
        n_ref - n_aff,
        0L,
        n_pairs - n_ref),
      nrow = 2, byrow = FALSE,
      dimnames = list(c("affected", "not_affected"),
                      c("has_domain", "no_domain"))
    )

    ft <- tryCatch(
      stats::fisher.test(contingency),
      error = function(e) NULL
    )

    all_results[[length(all_results) + 1L]] <- tibble::tibble(
      domain_id = did,
      domain_name = dname,
      n_ref_with_domain = n_ref,
      n_affected = n_aff,
      n_pairs_total = n_pairs,
      fisher_OR = if (!is.null(ft)) as.numeric(ft$estimate) else NA_real_,
      fisher_p = if (!is.null(ft)) ft$p.value else NA_real_,
      padj = NA_real_,
      significant = NA
    )
  }

  result <- dplyr::bind_rows(all_results)

  if (nrow(result) > 0L) {
    testable <- !is.na(result$fisher_p)
    if (any(testable)) {
      result$padj[testable] <- stats::p.adjust(
        result$fisher_p[testable], method = p_adjust_method
      )
      result$significant[testable] <- result$padj[testable] < 0.05
    }
    result <- result[order(result$fisher_p, na.last = TRUE), , drop = FALSE]
  }

  result
}
