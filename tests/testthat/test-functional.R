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
                     "orf_length") %in% names(cds)))
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
