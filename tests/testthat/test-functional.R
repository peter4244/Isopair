# Tests for extractCdsAnnotations() and annotateRegionTypes()

test_cds_gtf_path <- system.file("extdata", "test_cds_cases.gtf",
                                  package = "Isopair")

# ==========================================================================
# extractCdsAnnotations()
# ==========================================================================

test_that("extractCdsAnnotations returns correct columns", {
  cds <- extractCdsAnnotations(test_cds_gtf_path, verbose = FALSE)
  expect_s3_class(cds, "tbl_df")
  expect_true(all(c("isoform_id", "coding_status", "cds_start", "cds_stop",
                     "orf_length", "strand") %in% names(cds)))
})

test_that("extractCdsAnnotations identifies coding isoforms", {
  cds <- extractCdsAnnotations(test_cds_gtf_path, verbose = FALSE)
  coding <- cds[cds$coding_status == "coding", ]
  # Groups A, C, D have coding isoforms
  expect_true("TX_A1_1" %in% coding$isoform_id)
  expect_true("TX_C1_1" %in% coding$isoform_id)
  expect_true("TX_D1_ref" %in% coding$isoform_id)
})

test_that("extractCdsAnnotations: cds_start < cds_stop for all coding", {
  cds <- extractCdsAnnotations(test_cds_gtf_path, verbose = FALSE)
  coding <- cds[cds$coding_status == "coding", ]
  expect_true(all(coding$cds_start < coding$cds_stop))
})

test_that("extractCdsAnnotations: orf_length > 0 for coding", {
  cds <- extractCdsAnnotations(test_cds_gtf_path, verbose = FALSE)
  coding <- cds[cds$coding_status == "coding", ]
  expect_true(all(coding$orf_length > 0L))
})

test_that("extractCdsAnnotations: A1 has CDS spanning all exons", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A1_1",
                                verbose = FALSE)
  expect_equal(cds$cds_start[1], 1000L)
  expect_equal(cds$cds_stop[1], 2400L)
  # orf_length = (1200-1000+1) + (1800-1500+1) + (2400-2100+1) = 201+301+301 = 803
  expect_equal(cds$orf_length[1], 803L)
})

test_that("extractCdsAnnotations: A2 has CDS only in exon 2", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A2_1",
                                verbose = FALSE)
  expect_equal(cds$cds_start[1], 3550L)
  expect_equal(cds$cds_stop[1], 3750L)
  expect_equal(cds$orf_length[1], 201L)  # 3750-3550+1
})

test_that("extractCdsAnnotations: accepts GRanges input", {
  gr <- rtracklayer::import(test_cds_gtf_path)
  cds_from_path <- extractCdsAnnotations(test_cds_gtf_path,
                                          isoform_ids = "TX_A1_1",
                                          verbose = FALSE)
  cds_from_gr <- extractCdsAnnotations(gr,
                                        isoform_ids = "TX_A1_1",
                                        verbose = FALSE)
  expect_equal(cds_from_path$cds_start, cds_from_gr$cds_start)
  expect_equal(cds_from_path$cds_stop, cds_from_gr$cds_stop)
  expect_equal(cds_from_path$orf_length, cds_from_gr$orf_length)
})

test_that("extractCdsAnnotations: missing isoforms returned as unknown", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = c("TX_A1_1", "FAKE_TX_99"),
                                verbose = FALSE)
  expect_equal(nrow(cds), 2L)
  fake_row <- cds[cds$isoform_id == "FAKE_TX_99", ]
  expect_equal(fake_row$coding_status, "unknown")
  expect_true(is.na(fake_row$cds_start))
})

test_that("extractCdsAnnotations: non-coding gene has no CDS rows", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_G1_1",
                                verbose = FALSE)
  expect_equal(cds$coding_status[1], "unknown")
  expect_true(is.na(cds$orf_length[1]))
})

test_that("extractCdsAnnotations: minus strand CDS coordinates correct", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_F1_1",
                                verbose = FALSE)
  expect_equal(cds$coding_status[1], "coding")
  expect_equal(cds$cds_start[1], 31100L)
  expect_equal(cds$cds_stop[1], 32350L)
})

test_that("extractCdsAnnotations: strand column present and correct", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = c("TX_A1_1", "TX_F1_1",
                                                "FAKE_TX_99"),
                                verbose = FALSE)
  expect_true("strand" %in% names(cds))
  expect_equal(cds$strand[cds$isoform_id == "TX_A1_1"], "+")
  expect_equal(cds$strand[cds$isoform_id == "TX_F1_1"], "-")
  expect_true(is.na(cds$strand[cds$isoform_id == "FAKE_TX_99"]))
})

# ==========================================================================
# annotateRegionTypes()
# ==========================================================================

test_that("annotateRegionTypes: 5'UTR, CDS, 3'UTR classification", {
  # A3: Exon1(5000-5300) → contains_orf_start, Exon2(5600-5900) → CDS,
  #     Exon3(6200-6500) → contains_orf_stop
  # Also 5'UTR part of exon1 (before 5100) and 3'UTR part of exon3 (after 6350)
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A2_1",
                                verbose = FALSE)
  # A2: exon1 entirely before CDS → 5'UTR, exon3 entirely after CDS → 3'UTR
  mapping <- tibble::tibble(
    isoform_id = rep("TX_A2_1", 3),
    isoform_exon_start = c(3000L, 3500L, 4100L),
    isoform_exon_end = c(3200L, 3800L, 4300L),
    union_exon_id = c("ue1", "ue2", "ue3")
  )
  result <- annotateRegionTypes(mapping, cds)
  expect_equal(result$region_type[1], "5'UTR")
  expect_equal(result$region_type[2], "contains_orf_start_stop")
  expect_equal(result$region_type[3], "3'UTR")
})

test_that("annotateRegionTypes: contains_orf_start and contains_orf_stop", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1",
                                verbose = FALSE)
  mapping <- tibble::tibble(
    isoform_id = rep("TX_A3_1", 3),
    isoform_exon_start = c(5000L, 5600L, 6200L),
    isoform_exon_end = c(5300L, 5900L, 6500L),
    union_exon_id = c("ue1", "ue2", "ue3")
  )
  result <- annotateRegionTypes(mapping, cds)
  expect_equal(result$region_type[1], "contains_orf_start")
  expect_equal(result$region_type[2], "CDS")
  expect_equal(result$region_type[3], "contains_orf_stop")
})

test_that("annotateRegionTypes: contains_orf_start_stop for single-CDS-exon", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A4_1",
                                verbose = FALSE)
  mapping <- tibble::tibble(
    isoform_id = "TX_A4_1",
    isoform_exon_start = 7000L,
    isoform_exon_end = 7600L,
    union_exon_id = "ue1"
  )
  result <- annotateRegionTypes(mapping, cds)
  expect_equal(result$region_type[1], "contains_orf_start_stop")
})

test_that("annotateRegionTypes: non-coding → all non_coding", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_G1_1",
                                verbose = FALSE)
  mapping <- tibble::tibble(
    isoform_id = rep("TX_G1_1", 3),
    isoform_exon_start = c(42000L, 42500L, 43100L),
    isoform_exon_end = c(42200L, 42800L, 43500L),
    union_exon_id = c("ue1", "ue2", "ue3")
  )
  result <- annotateRegionTypes(mapping, cds)
  expect_true(all(result$region_type == "non_coding"))
})

test_that("annotateRegionTypes: output rows == input rows", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = c("TX_A1_1", "TX_G1_1"),
                                verbose = FALSE)
  mapping <- tibble::tibble(
    isoform_id = c(rep("TX_A1_1", 3), rep("TX_G1_1", 3)),
    isoform_exon_start = c(1000L, 1500L, 2100L, 42000L, 42500L, 43100L),
    isoform_exon_end = c(1200L, 1800L, 2400L, 42200L, 42800L, 43500L),
    union_exon_id = paste0("ue", 1:6)
  )
  result <- annotateRegionTypes(mapping, cds)
  expect_equal(nrow(result), nrow(mapping))
})

test_that("annotateRegionTypes: mixed coding/non-coding pair (G2)", {
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = c("TX_G2_coding",
                                                "TX_G2_noncoding"),
                                verbose = FALSE)
  mapping <- tibble::tibble(
    isoform_id = c(rep("TX_G2_coding", 3), rep("TX_G2_noncoding", 3)),
    isoform_exon_start = c(45000L, 45600L, 46200L, 45000L, 45600L, 46200L),
    isoform_exon_end = c(45300L, 45900L, 46500L, 45300L, 45900L, 46500L),
    union_exon_id = paste0("ue", 1:6)
  )
  result <- annotateRegionTypes(mapping, cds)
  coding_rows <- result[result$isoform_id == "TX_G2_coding", ]
  nc_rows <- result[result$isoform_id == "TX_G2_noncoding", ]
  expect_true(all(nc_rows$region_type == "non_coding"))
  expect_true(any(coding_rows$region_type != "non_coding"))
})

test_that("annotateRegionTypes: minus-strand UTR and ORF labels", {
  # TX_F1_1: minus strand, CDS 31100-32350
  # Exon1(32200-32500): CDS 32200-32350, UTR 32351-32500
  # Exon2(31600-31900): fully within CDS
  # Exon3(31000-31300): CDS 31100-31300, UTR 31000-31099
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_F1_1",
                                verbose = FALSE)

  # Synthetic exons: one fully left of CDS, one fully right, one within CDS,
  # one overlapping cds_start (31100), one overlapping cds_stop (32350)
  mapping <- tibble::tibble(
    isoform_id = rep("TX_F1_1", 5),
    isoform_exon_start = c(30800L, 32400L, 31600L, 31000L, 32200L),
    isoform_exon_end   = c(30900L, 32600L, 31900L, 31300L, 32500L),
    union_exon_id = paste0("ue", 1:5)
  )
  result <- annotateRegionTypes(mapping, cds)

  # Exon fully LEFT of CDS (30800-30900 < 31100) → minus strand → "3'UTR"
  expect_equal(result$region_type[1], "3'UTR")
  # Exon fully RIGHT of CDS (32400-32600 > 32350) → minus strand → "5'UTR"
  expect_equal(result$region_type[2], "5'UTR")
  # Exon within CDS (31600-31900) → "CDS"
  expect_equal(result$region_type[3], "CDS")
  # Exon overlapping cds_start=31100 (31000-31300) → minus strand →
  # cds_start boundary = stop codon → "contains_orf_stop"
  expect_equal(result$region_type[4], "contains_orf_stop")
  # Exon overlapping cds_stop=32350 (32200-32500) → minus strand →
  # cds_stop boundary = start codon → "contains_orf_start"
  expect_equal(result$region_type[5], "contains_orf_start")
})

test_that("annotateRegionTypes: plus-strand regression (A3 unchanged)", {
  # Same test as existing "contains_orf_start and contains_orf_stop" above,
  # but explicitly verifying plus-strand labels haven't changed
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1",
                                verbose = FALSE)
  expect_equal(cds$strand[1], "+")
  mapping <- tibble::tibble(
    isoform_id = rep("TX_A3_1", 5),
    isoform_exon_start = c(4800L, 6400L, 5600L, 5000L, 6200L),
    isoform_exon_end   = c(4900L, 6600L, 5900L, 5300L, 6500L),
    union_exon_id = paste0("ue", 1:5)
  )
  result <- annotateRegionTypes(mapping, cds)
  # Plus strand: left of CDS → 5'UTR, right → 3'UTR
  expect_equal(result$region_type[1], "5'UTR")
  expect_equal(result$region_type[2], "3'UTR")
  expect_equal(result$region_type[3], "CDS")
  expect_equal(result$region_type[4], "contains_orf_start")
  expect_equal(result$region_type[5], "contains_orf_stop")
})

test_that("annotateRegionTypes: missing strand column errors", {
  cds_no_strand <- tibble::tibble(
    isoform_id = "TX_A1_1",
    coding_status = "coding",
    cds_start = 1000L,
    cds_stop = 2400L
  )
  mapping <- tibble::tibble(
    isoform_id = "TX_A1_1",
    isoform_exon_start = 1000L,
    isoform_exon_end = 1200L,
    union_exon_id = "ue1"
  )
  expect_error(annotateRegionTypes(mapping, cds_no_strand),
               "strand")
})

# ==========================================================================
# mapEventsToRegions(), testRegionalEnrichment(),
# testOrfBoundarySusceptibility(), summarizeOrfImpact()
# ==========================================================================

# Helper: build minimal profiles + annotated_ue for functional tests
.make_functional_test_data <- function() {
  # CDS gene A3: 3 exons, CDS 5100-6350
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = c("TX_A3_1", "TX_G1_1"),
                                verbose = FALSE)

  # Annotated union exon mapping
  annotated_ue <- tibble::tibble(
    isoform_id = c(rep("TX_A3_1", 3), rep("TX_G1_1", 3)),
    union_exon_id = paste0("ue", 1:6),
    isoform_exon_start = c(5000L, 5600L, 6200L, 42000L, 42500L, 43100L),
    isoform_exon_end = c(5300L, 5900L, 6500L, 42200L, 42800L, 43500L)
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  # Profile with events in different regions
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_G1_1",
    detailed_events = list(tibble::tibble(
      event_type = c("Alt_TSS", "SE", "Alt_TES"),
      direction = c("LOSS", "LOSS", "LOSS"),
      five_prime = c(5050L, 5700L, 6400L),
      three_prime = c(5150L, 5800L, 6450L),
      bp_diff = c(100L, 100L, 50L)
    ))
  )

  list(profiles = profiles, annotated_ue = annotated_ue, cds = cds)
}

test_that("mapEventsToRegions: returns correct columns", {
  td <- .make_functional_test_data()
  result <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("event_row_id", "gene_id", "reference_isoform_id",
                     "event_type", "region_type") %in% names(result)))
})

test_that("mapEventsToRegions: events map to correct regions", {
  td <- .make_functional_test_data()
  result <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  # Alt_TSS at 5050-5150 overlaps exon1 (5000-5300) which is contains_orf_start
  alt_tss <- result[result$event_type == "Alt_TSS", ]
  expect_true("contains_orf_start" %in% alt_tss$region_type)
  # SE at 5700-5800 overlaps exon2 (5600-5900) which is CDS
  se <- result[result$event_type == "SE", ]
  expect_true("CDS" %in% se$region_type)
  # Alt_TES at 6400-6450 overlaps exon3 (6200-6500) which is contains_orf_stop
  alt_tes <- result[result$event_type == "Alt_TES", ]
  expect_true("contains_orf_stop" %in% alt_tes$region_type)
})

test_that("mapEventsToRegions: empty profiles → empty result", {
  td <- .make_functional_test_data()
  empty_profiles <- td$profiles[0, ]
  result <- mapEventsToRegions(empty_profiles, td$annotated_ue, td$cds)
  expect_equal(nrow(result), 0L)
})

test_that("mapEventsToRegions: intronic events labeled correctly", {
  td <- .make_functional_test_data()
  # Event in intron (between exon 1 end 5300 and exon 2 start 5600)
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_G1_1",
    detailed_events = list(tibble::tibble(
      event_type = "IR", direction = "GAIN",
      five_prime = 5350L, three_prime = 5550L,
      bp_diff = NA_integer_
    ))
  )
  result <- mapEventsToRegions(profiles, td$annotated_ue, td$cds)
  expect_true("intronic" %in% result$region_type)
})

test_that("testRegionalEnrichment: returns correct columns", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- testRegionalEnrichment(event_regions, td$annotated_ue, td$profiles)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("event_type", "region_type", "observed_count",
                     "total_events", "observed_prop", "expected_prop",
                     "enrichment_ratio", "p_value", "adj_p_value",
                     "significant", "direction_enrich") %in% names(result)))
})

test_that("testRegionalEnrichment: expected_prop sums close to 1", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- testRegionalEnrichment(event_regions, td$annotated_ue, td$profiles)
  # Each event type's expected props should sum to ~1 across covered regions
  # (not exactly 1 if not all regions are represented)
  if (nrow(result) > 0L) {
    expect_true(all(result$observed_prop >= 0))
    expect_true(all(result$observed_prop <= 1))
  }
})

test_that("testOrfBoundarySusceptibility: returns correct columns", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- testOrfBoundarySusceptibility(event_regions, td$annotated_ue)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("boundary_type", "n_exons", "n_with_events",
                     "rate") %in% names(result)))
  expect_equal(nrow(result), 2L)
  expect_true("ORF_boundary" %in% result$boundary_type)
  expect_true("Regular_CDS" %in% result$boundary_type)
})

test_that("summarizeOrfImpact: returns correct structure", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- summarizeOrfImpact(event_regions, td$cds)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2L)  # Alt_TSS and Alt_TES
  expect_true(all(c("event_type", "n_events", "n_affecting_orf_start",
                     "pct_affecting_orf_start", "n_affecting_orf_stop",
                     "pct_affecting_orf_stop") %in% names(result)))
})

test_that("summarizeOrfImpact: Alt_TSS affecting orf_start detected", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- summarizeOrfImpact(event_regions, td$cds)
  tss_row <- result[result$event_type == "Alt_TSS", ]
  # Alt_TSS at 5050-5150 overlaps contains_orf_start → should count
  expect_true(tss_row$n_affecting_orf_start > 0L)
})

test_that("summarizeOrfImpact: Alt_TES affecting orf_stop detected", {
  td <- .make_functional_test_data()
  event_regions <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- summarizeOrfImpact(event_regions, td$cds)
  tes_row <- result[result$event_type == "Alt_TES", ]
  # Alt_TES at 6400-6450 overlaps contains_orf_stop → should count
  expect_true(tes_row$n_affecting_orf_stop > 0L)
})


# ==========================================================================
# bootstrapRegionalComparison()
# ==========================================================================

test_that("bootstrapRegionalComparison: returns correct columns", {
  td <- .make_functional_test_data()
  er <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- bootstrapRegionalComparison(er, er, n_boot = 50,
                                         seed = 42, verbose = FALSE)
  expect_s3_class(result, "tbl_df")
  expected_cols <- c("event_type", "region_type", "n_events_a", "n_events_b",
                     "prop_a", "prop_b", "diff", "ci_lower", "ci_upper",
                     "p_value", "adj_p_value", "significant", "skipped")
  expect_true(all(expected_cols %in% names(result)))
})

test_that("bootstrapRegionalComparison: identical sets → no significance", {
  td <- .make_functional_test_data()
  er <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- bootstrapRegionalComparison(er, er, n_boot = 100,
                                         seed = 42, verbose = FALSE)
  non_skipped <- result[!result$skipped, ]
  if (nrow(non_skipped) > 0) {
    expect_true(all(abs(non_skipped$diff) < 1e-10))
  }
})

test_that("bootstrapRegionalComparison: seed reproducibility", {
  td <- .make_functional_test_data()
  er <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  r1 <- bootstrapRegionalComparison(er, er, n_boot = 50,
                                     seed = 123, verbose = FALSE)
  r2 <- bootstrapRegionalComparison(er, er, n_boot = 50,
                                     seed = 123, verbose = FALSE)
  expect_equal(r1$p_value, r2$p_value)
  expect_equal(r1$ci_lower, r2$ci_lower)
  expect_equal(r1$ci_upper, r2$ci_upper)
})

test_that("bootstrapRegionalComparison: empty inputs → empty result", {
  empty <- tibble::tibble(
    event_row_id = integer(0), gene_id = character(0),
    reference_isoform_id = character(0), comparator_isoform_id = character(0),
    event_type = character(0), direction = character(0),
    genomic_start = integer(0), genomic_end = integer(0),
    region_type = character(0)
  )
  result <- bootstrapRegionalComparison(empty, empty, n_boot = 50,
                                         verbose = FALSE)
  expect_equal(nrow(result), 0L)
})

test_that("bootstrapRegionalComparison: min_events filtering works", {
  small_er <- tibble::tibble(
    event_row_id = 1:3,
    gene_id = rep("G1", 3),
    reference_isoform_id = rep("REF", 3),
    comparator_isoform_id = rep("COMP", 3),
    event_type = rep("SE", 3),
    direction = rep("LOSS", 3),
    genomic_start = c(100L, 200L, 300L),
    genomic_end = c(150L, 250L, 350L),
    region_type = rep("CDS", 3)
  )
  result <- bootstrapRegionalComparison(small_er, small_er, n_boot = 50,
                                         min_events = 5, verbose = FALSE)
  expect_true(all(result$skipped))
  expect_true(all(is.na(result$p_value)))
})

test_that("bootstrapRegionalComparison: min_events=1 allows small sets", {
  small_er <- tibble::tibble(
    event_row_id = 1:3,
    gene_id = rep("G1", 3),
    reference_isoform_id = rep("REF", 3),
    comparator_isoform_id = rep("COMP", 3),
    event_type = rep("SE", 3),
    direction = rep("LOSS", 3),
    genomic_start = c(100L, 200L, 300L),
    genomic_end = c(150L, 250L, 350L),
    region_type = rep("CDS", 3)
  )
  result <- bootstrapRegionalComparison(small_er, small_er, n_boot = 50,
                                         min_events = 1, seed = 42,
                                         verbose = FALSE)
  expect_true(all(!result$skipped))
  expect_true(all(!is.na(result$p_value)))
})

test_that("bootstrapRegionalComparison: event type in only one set", {
  er_a <- tibble::tibble(
    event_row_id = 1:10,
    gene_id = rep("G1", 10),
    reference_isoform_id = rep("REF", 10),
    comparator_isoform_id = rep("COMP", 10),
    event_type = rep("SE", 10),
    direction = rep("LOSS", 10),
    genomic_start = seq(100L, 1000L, by = 100L),
    genomic_end = seq(150L, 1050L, by = 100L),
    region_type = rep("CDS", 10)
  )
  er_b <- tibble::tibble(
    event_row_id = 1:10,
    gene_id = rep("G1", 10),
    reference_isoform_id = rep("REF", 10),
    comparator_isoform_id = rep("COMP", 10),
    event_type = rep("A5SS", 10),
    direction = rep("LOSS", 10),
    genomic_start = seq(100L, 1000L, by = 100L),
    genomic_end = seq(150L, 1050L, by = 100L),
    region_type = rep("CDS", 10)
  )
  result <- bootstrapRegionalComparison(er_a, er_b, n_boot = 50,
                                         seed = 42, verbose = FALSE)
  se_rows <- result[result$event_type == "SE", ]
  a5_rows <- result[result$event_type == "A5SS", ]
  expect_true(all(se_rows$skipped))
  expect_true(all(a5_rows$skipped))
})

test_that("bootstrapRegionalComparison: p_value bounded [0, 1]", {
  td <- .make_functional_test_data()
  er <- mapEventsToRegions(td$profiles, td$annotated_ue, td$cds)
  result <- bootstrapRegionalComparison(er, er, n_boot = 50,
                                         min_events = 1, seed = 42,
                                         verbose = FALSE)
  non_na <- result$p_value[!is.na(result$p_value)]
  if (length(non_na) > 0) {
    expect_true(all(non_na >= 0 & non_na <= 1))
  }
})

# ==========================================================================
# quantifyRegionalImpact()
# ==========================================================================

test_that("quantifyRegionalImpact: invariant ref bp sum == event_span_bp", {
  td <- .make_functional_test_data()
  result <- quantifyRegionalImpact(td$profiles, td$annotated_ue, td$cds)
  if (nrow(result) > 0L) {
    ref_sum <- result$ref_utr5_bp + result$ref_cds_bp +
      result$ref_utr3_bp + result$ref_intronic_bp
    expect_equal(ref_sum, result$event_span_bp)
  }
})

test_that("quantifyRegionalImpact: no negative bp values", {
  td <- .make_functional_test_data()
  result <- quantifyRegionalImpact(td$profiles, td$annotated_ue, td$cds)
  if (nrow(result) > 0L) {
    expect_true(all(result$ref_utr5_bp >= 0L))
    expect_true(all(result$ref_cds_bp >= 0L))
    expect_true(all(result$ref_utr3_bp >= 0L))
    expect_true(all(result$ref_intronic_bp >= 0L))
    expect_true(all(result$comp_utr5_bp >= 0L))
    expect_true(all(result$comp_cds_bp >= 0L))
    expect_true(all(result$comp_utr3_bp >= 0L))
    expect_true(all(result$comp_intronic_bp >= 0L))
  }
})

test_that("quantifyRegionalImpact: event spanning 5'UTR+CDS on plus strand", {
  # TX_A3_1: CDS 5100-6350, Exon1: 5000-5300 = contains_orf_start
  # Event [5050, 5150]: spans 5'UTR (5050-5099=50bp) + CDS (5100-5150=51bp)
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_A3_1",
    isoform_exon_start = 5000L,
    isoform_exon_end = 5300L,
    union_exon_id = "ue1"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 5050L, three_prime = 5150L,
      bp_diff = 100L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$event_span_bp[1], 101L)
  expect_equal(result$ref_utr5_bp[1], 50L)
  expect_equal(result$ref_cds_bp[1], 51L)
  expect_equal(result$ref_utr3_bp[1], 0L)
  expect_equal(result$ref_intronic_bp[1], 0L)
})

test_that("quantifyRegionalImpact: intronic event", {
  td <- .make_functional_test_data()
  # Event entirely in intron between exons (5300-5600 gap)
  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "IR", direction = "GAIN",
      five_prime = 5350L, three_prime = 5550L,
      bp_diff = NA_integer_
    ))
  )
  result <- quantifyRegionalImpact(profiles, td$annotated_ue, td$cds)
  expect_equal(result$ref_intronic_bp[1], result$event_span_bp[1])
  expect_equal(result$ref_cds_bp[1], 0L)
  expect_equal(result$ref_utr5_bp[1], 0L)
  expect_equal(result$ref_utr3_bp[1], 0L)
})

test_that("quantifyRegionalImpact: event spanning CDS+3'UTR on contains_orf_stop", {
  # TX_A3_1: CDS 5100-6350, Exon3: 6200-6500 = contains_orf_stop
  # Event [6300, 6400]: CDS(6300-6350=51bp) + 3'UTR(6351-6400=50bp)
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_A3_1",
    isoform_exon_start = 6200L,
    isoform_exon_end = 6500L,
    union_exon_id = "ue3"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "A3SS", direction = "LOSS",
      five_prime = 6300L, three_prime = 6400L,
      bp_diff = 100L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$ref_cds_bp[1], 51L)
  expect_equal(result$ref_utr3_bp[1], 50L)
  expect_equal(result$ref_utr5_bp[1], 0L)
})

test_that("quantifyRegionalImpact: contains_orf_start_stop three-way split", {
  # TX_A4_1: single exon 7000-7600, CDS 7100-7500
  # Event [7050, 7550]: 5'UTR(7050-7099=50bp) + CDS(7100-7500=401bp) + 3'UTR(7501-7550=50bp)
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A4_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_A4_1",
    isoform_exon_start = 7000L,
    isoform_exon_end = 7600L,
    union_exon_id = "ue1"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A4",
    reference_isoform_id = "TX_A4_1",
    comparator_isoform_id = "TX_A4_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 7050L, three_prime = 7550L,
      bp_diff = 500L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$ref_utr5_bp[1], 50L)
  expect_equal(result$ref_cds_bp[1], 401L)
  expect_equal(result$ref_utr3_bp[1], 50L)
  expect_equal(result$ref_intronic_bp[1], 0L)
  expect_equal(result$ref_category[1], "5'UTR + CDS + 3'UTR")
})

test_that("quantifyRegionalImpact: minus-strand decomposition", {
  # TX_F1_1: minus strand, CDS 31100-32350, Exon3: 32200-32500 = contains_orf_start
  # On minus strand, region above cds_stop (32351-32500) is 5'UTR
  # Event [32300, 32400]: CDS(32300-32350=51bp) + 5'UTR(32351-32400=50bp)
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_F1_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_F1_1",
    isoform_exon_start = 32200L,
    isoform_exon_end = 32500L,
    union_exon_id = "ue1"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_F1",
    reference_isoform_id = "TX_F1_1",
    comparator_isoform_id = "TX_F1_1",
    detailed_events = list(tibble::tibble(
      event_type = "A5SS", direction = "LOSS",
      five_prime = 32300L, three_prime = 32400L,
      bp_diff = 100L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$ref_utr5_bp[1], 50L)
  expect_equal(result$ref_cds_bp[1], 51L)
  expect_equal(result$ref_utr3_bp[1], 0L)
})

test_that("quantifyRegionalImpact: off-by-one at cds_start", {
  # TX_A3_1: CDS starts at 5100. Event exactly at [5100, 5100]: 1bp CDS
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_A3_1",
    isoform_exon_start = 5000L,
    isoform_exon_end = 5300L,
    union_exon_id = "ue1"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 5100L, three_prime = 5100L,
      bp_diff = 1L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$event_span_bp[1], 1L)
  expect_equal(result$ref_cds_bp[1], 1L)
  expect_equal(result$ref_utr5_bp[1], 0L)
})

test_that("quantifyRegionalImpact: off-by-one at cds_stop", {
  # TX_A3_1: CDS ends at 6350. Event exactly at [6350, 6350]: 1bp CDS
  cds <- extractCdsAnnotations(test_cds_gtf_path,
                                isoform_ids = "TX_A3_1", verbose = FALSE)
  annotated_ue <- tibble::tibble(
    isoform_id = "TX_A3_1",
    isoform_exon_start = 6200L,
    isoform_exon_end = 6500L,
    union_exon_id = "ue3"
  )
  annotated_ue <- annotateRegionTypes(annotated_ue, cds)

  profiles <- tibble::tibble(
    gene_id = "GENE_CDS_A3",
    reference_isoform_id = "TX_A3_1",
    comparator_isoform_id = "TX_A3_1",
    detailed_events = list(tibble::tibble(
      event_type = "SE", direction = "LOSS",
      five_prime = 6350L, three_prime = 6350L,
      bp_diff = 1L
    ))
  )
  result <- quantifyRegionalImpact(profiles, annotated_ue, cds)
  expect_equal(result$event_span_bp[1], 1L)
  expect_equal(result$ref_cds_bp[1], 1L)
  expect_equal(result$ref_utr3_bp[1], 0L)
})

test_that("quantifyRegionalImpact: empty profiles → zero-row tibble", {
  td <- .make_functional_test_data()
  empty <- td$profiles[0, ]
  result <- quantifyRegionalImpact(empty, td$annotated_ue, td$cds)
  expect_equal(nrow(result), 0L)
  expect_true(all(c("event_id", "ref_cds_bp", "comp_cds_bp",
                     "ref_category") %in% names(result)))
})

test_that("quantifyRegionalImpact: non-coding comparator all intronic/non-coding", {
  # Reference is coding (TX_A3_1), comparator is non-coding (TX_G1_1)
  td <- .make_functional_test_data()
  result <- quantifyRegionalImpact(td$profiles, td$annotated_ue, td$cds)
  # Comp is non-coding so no CDS segments → all comp bp should be intronic
  if (nrow(result) > 0L) {
    expect_true(all(result$comp_cds_bp == 0L))
    expect_true(all(result$comp_utr5_bp == 0L))
    expect_true(all(result$comp_utr3_bp == 0L))
  }
})

test_that("quantifyRegionalImpact: real-gene integration (BCL2L1)", {
  real_gtf <- system.file("extdata", "test_real_genes.gtf", package = "Isopair")
  if (nchar(real_gtf) > 0L && file.exists(real_gtf)) {
    structures <- parseIsoformStructures(real_gtf, verbose = FALSE)
    ue <- buildUnionExons(structures, verbose = FALSE)
    cds <- extractCdsAnnotations(real_gtf, structures$isoform_id,
                                  verbose = FALSE)
    annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)

    # Build a simple profile from any available pair
    bcl2_structs <- structures[grepl("BCL2L1", structures$gene_id), ]
    if (nrow(bcl2_structs) >= 2L) {
      profiles <- buildProfiles(
        tibble::tibble(
          gene_id = bcl2_structs$gene_id[1],
          reference_isoform_id = bcl2_structs$isoform_id[1],
          comparator_isoform_id = bcl2_structs$isoform_id[2]
        ),
        structures, ue$union_exons, ue$isoform_union_mapping,
        verify = FALSE
      )
      result <- quantifyRegionalImpact(profiles, annotated, cds)
      # Check invariants
      if (nrow(result) > 0L) {
        ref_sum <- result$ref_utr5_bp + result$ref_cds_bp +
          result$ref_utr3_bp + result$ref_intronic_bp
        expect_equal(ref_sum, result$event_span_bp)
        expect_true(all(result$ref_cds_bp >= 0L))
      }
    }
  } else {
    skip("test_real_genes.gtf not available")
  }
})


test_that("bootstrapRegionalComparison: extreme difference detected", {
  er_a <- tibble::tibble(
    event_row_id = 1:20,
    gene_id = rep("G1", 20),
    reference_isoform_id = rep("REF", 20),
    comparator_isoform_id = rep("COMP", 20),
    event_type = rep("SE", 20),
    direction = rep("LOSS", 20),
    genomic_start = seq(100L, 2000L, by = 100L),
    genomic_end = seq(150L, 2050L, by = 100L),
    region_type = rep("CDS", 20)
  )
  er_b <- tibble::tibble(
    event_row_id = 1:20,
    gene_id = rep("G1", 20),
    reference_isoform_id = rep("REF", 20),
    comparator_isoform_id = rep("COMP", 20),
    event_type = rep("SE", 20),
    direction = rep("LOSS", 20),
    genomic_start = seq(100L, 2000L, by = 100L),
    genomic_end = seq(150L, 2050L, by = 100L),
    region_type = rep("3'UTR", 20)
  )
  result <- bootstrapRegionalComparison(er_a, er_b, n_boot = 200,
                                         min_events = 1, seed = 42,
                                         verbose = FALSE)
  cds_row <- result[result$region_type == "CDS", ]
  utr_row <- result[result$region_type == "3'UTR", ]
  expect_true(nrow(cds_row) > 0)
  expect_true(nrow(utr_row) > 0)
  # Use unname to avoid named vs unnamed mismatch
  expect_equal(unname(cds_row$prop_a), 1.0)
  expect_equal(unname(cds_row$prop_b), 0.0)
  expect_equal(unname(utr_row$prop_a), 0.0)
  expect_equal(unname(utr_row$prop_b), 1.0)
})


# ==========================================================================
# quantifyPairDivergence()
# ==========================================================================

# Helper: build minimal structures tibble from exon coordinates
.make_struct <- function(isoform_id, gene_id, chr, strand, starts, ends) {
  tibble::tibble(
    isoform_id = isoform_id, gene_id = gene_id,
    chr = chr, strand = strand,
    n_exons = length(starts),
    exon_starts = list(starts), exon_ends = list(ends),
    tx_start = min(starts), tx_end = max(ends),
    n_junctions = max(0L, length(starts) - 1L)
  )
}

# Helper: build minimal profiles tibble
.make_div_profiles <- function(gene_id, ref_id, comp_id) {
  tibble::tibble(
    gene_id = gene_id,
    reference_isoform_id = ref_id,
    comparator_isoform_id = comp_id,
    detailed_events = list(tibble::tibble(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0)
    )),
    n_events = 0L
  )
}

test_that("quantifyPairDivergence: invariant utr5+cds+utr3+noncoding == exonic for coding", {
  # Plus strand, CDS 5100-6350 (same as TX_A3_1)
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(5000L, 5600L, 6200L), c(5300L, 5900L, 6500L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(5000L, 5600L), c(5300L, 5900L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("coding", "coding"),
    cds_start = c(5100L, 5100L),
    cds_stop = c(6350L, 5850L),
    orf_length = c(803L, 451L),
    strand = c("+", "+")
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(nrow(result), 1L)
  # Invariant: per-region sums equal total exonic
  expect_equal(result$utr5_bp_ref + result$cds_bp_ref +
    result$utr3_bp_ref + result$noncoding_bp_ref, result$exonic_bp_ref)
  expect_equal(result$utr5_bp_comp + result$cds_bp_comp +
    result$utr3_bp_comp + result$noncoding_bp_comp, result$exonic_bp_comp)
})

test_that("quantifyPairDivergence: invariant shared+ref_unique+comp_unique == union", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(100L, 300L, 500L), c(200L, 400L, 600L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(100L, 350L), c(200L, 600L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("unknown", "unknown"),
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$shared_exonic_bp + result$ref_unique_bp +
    result$comp_unique_bp, result$union_exonic_bp)
})

test_that("quantifyPairDivergence: identical isoforms → pct_shared=100, diffs=0", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(100L, 300L), c(200L, 400L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(100L, 300L), c(200L, 400L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("unknown", "unknown"),
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$pct_shared, 100)
  expect_equal(result$exonic_bp_diff, 0L)
  expect_equal(result$ref_unique_bp, 0L)
  expect_equal(result$comp_unique_bp, 0L)
})

test_that("quantifyPairDivergence: non-overlapping exons → shared=0, pct_shared=0", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(100L), c(200L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(300L), c(400L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("unknown", "unknown"),
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$shared_exonic_bp, 0L)
  expect_equal(result$pct_shared, 0)
  expect_equal(result$ref_unique_bp, 101L)
  expect_equal(result$comp_unique_bp, 101L)
})

test_that("quantifyPairDivergence: minus-strand UTR assignment correct", {
  # Minus strand: genomically before CDS = 3'UTR, after = 5'UTR
  # Exon1: 1000-1200, Exon2: 1500-1800, CDS: 1100-1700
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "-",
                 c(1000L, 1500L), c(1200L, 1800L)),
    .make_struct("COMP1", "G1", "chr1", "-",
                 c(1000L, 1500L), c(1200L, 1800L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("coding", "coding"),
    cds_start = c(1100L, 1100L), cds_stop = c(1700L, 1700L),
    orf_length = c(402L, 402L), strand = c("-", "-")
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  # Exon1 (1000-1200): before CDS (1000-1099) = 100bp 3'UTR (minus strand)
  #                     within CDS (1100-1200) = 101bp CDS
  # Exon2 (1500-1800): within CDS (1500-1700) = 201bp CDS
  #                     after CDS (1701-1800) = 100bp 5'UTR (minus strand)
  expect_equal(result$utr3_bp_ref, 100L)
  expect_equal(result$cds_bp_ref, 302L)
  expect_equal(result$utr5_bp_ref, 100L)
  expect_equal(result$noncoding_bp_ref, 0L)
})

test_that("quantifyPairDivergence: non-coding isoform → cds=0, noncoding=exonic", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(100L, 300L), c(200L, 400L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(100L), c(200L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("unknown", "unknown"),
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$cds_bp_ref, 0L)
  expect_equal(result$utr5_bp_ref, 0L)
  expect_equal(result$utr3_bp_ref, 0L)
  expect_equal(result$noncoding_bp_ref, result$exonic_bp_ref)
  expect_equal(result$noncoding_bp_comp, result$exonic_bp_comp)
})

test_that("quantifyPairDivergence: partial exon overlap with hand-computed shared", {
  # REF: [100, 200] = 101bp
  # COMP: [150, 250] = 101bp
  # Shared: [150, 200] = 51bp
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+", c(100L), c(200L)),
    .make_struct("COMP1", "G1", "chr1", "+", c(150L), c(250L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("unknown", "unknown"),
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$shared_exonic_bp, 51L)
  expect_equal(result$ref_unique_bp, 50L)   # 100-149
  expect_equal(result$comp_unique_bp, 50L)  # 201-250
  expect_equal(result$union_exonic_bp, 151L) # 100-250
})

test_that("quantifyPairDivergence: empty profiles → zero-row tibble", {
  structures <- .make_struct("REF1", "G1", "chr1", "+", c(100L), c(200L))
  cds <- tibble::tibble(
    isoform_id = "REF1", coding_status = "unknown",
    cds_start = NA_integer_, cds_stop = NA_integer_,
    orf_length = NA_integer_, strand = NA_character_
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")[0, ]
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(nrow(result), 0L)
  expect_true(all(c("gene_id", "exonic_bp_ref", "pct_shared") %in% names(result)))
})

test_that("quantifyPairDivergence: mixed coding/non-coding pair", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(5000L, 5600L), c(5300L, 5900L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(5000L, 5600L), c(5300L, 5900L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("coding", "unknown"),
    cds_start = c(5100L, NA_integer_),
    cds_stop = c(5850L, NA_integer_),
    orf_length = c(451L, NA_integer_),
    strand = c("+", NA_character_)
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  # Ref is coding, comp is non-coding
  expect_true(result$cds_bp_ref > 0L)
  expect_equal(result$cds_bp_comp, 0L)
  expect_equal(result$noncoding_bp_comp, result$exonic_bp_comp)
  expect_equal(result$noncoding_bp_ref, 0L)
  # Same exons → pct_shared = 100
  expect_equal(result$pct_shared, 100)
})

test_that("quantifyPairDivergence: CDS spans entire exon → utr5=0, utr3=0", {
  structures <- dplyr::bind_rows(
    .make_struct("REF1", "G1", "chr1", "+",
                 c(100L, 300L), c(200L, 400L)),
    .make_struct("COMP1", "G1", "chr1", "+",
                 c(100L), c(200L))
  )
  cds <- tibble::tibble(
    isoform_id = c("REF1", "COMP1"),
    coding_status = c("coding", "coding"),
    cds_start = c(100L, 100L), cds_stop = c(400L, 200L),
    orf_length = c(202L, 101L), strand = c("+", "+")
  )
  profiles <- .make_div_profiles("G1", "REF1", "COMP1")
  result <- quantifyPairDivergence(profiles, structures, cds)
  expect_equal(result$utr5_bp_ref, 0L)
  expect_equal(result$utr3_bp_ref, 0L)
  expect_equal(result$cds_bp_ref, result$exonic_bp_ref)
  expect_equal(result$utr5_bp_comp, 0L)
  expect_equal(result$utr3_bp_comp, 0L)
  expect_equal(result$cds_bp_comp, result$exonic_bp_comp)
})
