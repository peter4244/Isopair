# Tests for classifyDomainEffects() and testDomainEnrichment()

# ============================================================================
# Helper: build minimal test data
# ============================================================================

.make_protein_events <- function() {
  tibble::tibble(
    gene_id = c("gene1", "gene1"),
    reference_isoform_id = c("tx1", "tx1"),
    comparator_isoform_id = c("tx2", "tx2"),
    event_type = c("SE", "SE"),
    direction = c("LOSS", "GAIN"),
    genomic_start = c(350L, 600L),
    genomic_end = c(450L, 700L),
    cds_nt_start = c(100L, 250L),
    cds_nt_end = c(200L, 350L),
    protein_start_aa = c(34L, 84L),
    protein_end_aa = c(67L, 117L),
    protein_start_aa_buffered = c(29L, 79L),
    protein_end_aa_buffered = c(72L, 122L),
    event_affects_frame = c(FALSE, FALSE)
  )
}

.make_domain_annotations <- function() {
  tibble::tibble(
    isoform_id = c("tx1", "tx1", "tx1", "tx2", "tx2"),
    domain_id = c("PF00091", "PF00200", "PF00300", "PF00091", "PF00400"),
    domain_name = c("Tubulin", "Kinase", "Zinc_finger", "Tubulin", "SH3"),
    domain_start = c(40L, 10L, 150L, 90L, 100L),
    domain_end = c(60L, 20L, 200L, 110L, 130L)
  )
}


# ============================================================================
# classifyDomainEffects()
# ============================================================================

test_that("classifyDomainEffects returns correct effects for LOSS events", {
  events <- .make_protein_events()[1L, ]  # LOSS event: AA 34-67
  domains <- .make_domain_annotations()

  result <- classifyDomainEffects(events, domains)

  # tx1 Tubulin (40-60) fully within 34-67 -> "lost"
  tubulin <- result[result$domain_id == "PF00091", ]
  expect_equal(nrow(tubulin), 1L)
  expect_equal(tubulin$effect, "lost")
  expect_equal(tubulin$isoform_source, "reference")
  expect_equal(tubulin$overlap_start, 40L)
  expect_equal(tubulin$overlap_end, 60L)

  # tx1 Kinase (10-20) does NOT overlap 34-67 -> excluded
  expect_equal(sum(result$domain_id == "PF00200"), 0L)

  # tx1 Zinc_finger (150-200) does NOT overlap 34-67 -> excluded
  expect_equal(sum(result$domain_id == "PF00300"), 0L)

  # tx2 domains should NOT be checked for LOSS events
  expect_equal(sum(result$isoform_source == "comparator"), 0L)
})

test_that("classifyDomainEffects classifies disrupted for partial LOSS overlap", {
  # LOSS event: AA 34-67, domain spans 50-80 (partially overlaps)
  events <- .make_protein_events()[1L, ]
  domains <- tibble::tibble(
    isoform_id = "tx1",
    domain_id = "PF00500",
    domain_name = "Partial_domain",
    domain_start = 50L,
    domain_end = 80L
  )

  result <- classifyDomainEffects(events, domains)
  expect_equal(nrow(result), 1L)
  expect_equal(result$effect, "disrupted")
  expect_equal(result$overlap_start, 50L)
  expect_equal(result$overlap_end, 67L)
})

test_that("classifyDomainEffects returns gained_intact for GAIN events", {
  events <- .make_protein_events()[2L, ]  # GAIN event: AA 84-117
  # tx2 SH3 domain: 100-130 partially overlaps -> gained_partial
  # tx2 Tubulin domain: 90-110 partially overlaps -> gained_partial
  domains <- .make_domain_annotations()

  result <- classifyDomainEffects(events, domains)

  # Only tx2 (comparator) domains should be checked for GAIN
  expect_true(all(result$isoform_source == "comparator"))

  # tx2 Tubulin (90-110): not entirely within 84-117 -> gained_partial?
  # Actually 90-110 IS within 84-117. 90 >= 84 and 110 <= 117. -> gained_intact
  tubulin <- result[result$domain_id == "PF00091", ]
  expect_equal(nrow(tubulin), 1L)
  expect_equal(tubulin$effect, "gained_intact")

  # tx2 SH3 (100-130): 130 > 117, partially overlaps -> gained_partial

  sh3 <- result[result$domain_id == "PF00400", ]
  expect_equal(nrow(sh3), 1L)
  expect_equal(sh3$effect, "gained_partial")
  expect_equal(sh3$overlap_start, 100L)
  expect_equal(sh3$overlap_end, 117L)
})

test_that("classifyDomainEffects handles non-overlapping domains", {
  events <- .make_protein_events()[1L, ]  # LOSS AA 34-67
  domains <- tibble::tibble(
    isoform_id = "tx1",
    domain_id = "PF99999",
    domain_name = "Far_away",
    domain_start = 200L,
    domain_end = 300L
  )

  result <- classifyDomainEffects(events, domains)
  expect_equal(nrow(result), 0L)
})

test_that("classifyDomainEffects returns zero-row tibble for empty inputs", {
  events <- .make_protein_events()[0L, ]
  domains <- .make_domain_annotations()

  result <- classifyDomainEffects(events, domains)
  expect_equal(nrow(result), 0L)
  expect_true("domain_id" %in% names(result))
  expect_true("effect" %in% names(result))
})

test_that("classifyDomainEffects handles NA domain positions", {
  events <- .make_protein_events()[1L, ]
  domains <- tibble::tibble(
    isoform_id = "tx1",
    domain_id = c("PF00001", "PF00002"),
    domain_name = c("Good", "Bad"),
    domain_start = c(40L, NA_integer_),
    domain_end = c(60L, NA_integer_)
  )

  result <- classifyDomainEffects(events, domains)
  # Only PF00001 should appear (PF00002 has NA positions)
  expect_equal(nrow(result), 1L)
  expect_equal(result$domain_id, "PF00001")
})

test_that("classifyDomainEffects handles NA protein coordinates in events", {
  events <- .make_protein_events()[1L, ]
  events$protein_start_aa <- NA_integer_

  domains <- .make_domain_annotations()
  result <- classifyDomainEffects(events, domains)
  expect_equal(nrow(result), 0L)
})

test_that("classifyDomainEffects preserves all protein_events columns", {
  events <- .make_protein_events()[1L, ]
  domains <- tibble::tibble(
    isoform_id = "tx1",
    domain_id = "PF00091",
    domain_name = "Tubulin",
    domain_start = 40L,
    domain_end = 60L
  )

  result <- classifyDomainEffects(events, domains)
  # Check that all original columns from protein_events are present
  for (col in names(events)) {
    expect_true(col %in% names(result),
                info = paste("Missing column:", col))
  }
})

test_that("classifyDomainEffects validates input columns", {
  events <- tibble::tibble(gene_id = "g1")
  domains <- .make_domain_annotations()

  expect_error(classifyDomainEffects(events, domains),
               "missing required columns")
})

test_that("classifyDomainEffects validates input types", {
  expect_error(classifyDomainEffects("not_df", data.frame()),
               "must be a data frame")
  expect_error(classifyDomainEffects(data.frame(), "not_df"),
               "must be a data frame")
})

test_that("classifyDomainEffects handles multiple events and domains", {
  events <- .make_protein_events()  # both LOSS and GAIN
  domains <- .make_domain_annotations()

  result <- classifyDomainEffects(events, domains)
  expect_true(nrow(result) > 0L)

  # Should have both reference and comparator sources
  expect_true("reference" %in% result$isoform_source)
  expect_true("comparator" %in% result$isoform_source)
})


# ============================================================================
# testDomainEnrichment()
# ============================================================================

test_that("testDomainEnrichment returns expected structure", {
  domain_effects <- tibble::tibble(
    domain_id = c("PF00091", "PF00091"),
    domain_name = c("Tubulin", "Tubulin"),
    reference_isoform_id = c("tx1", "tx2"),
    comparator_isoform_id = c("tx1c", "tx2c"),
    effect = c("lost", "disrupted")
  )
  domain_annotations <- tibble::tibble(
    isoform_id = c("tx1", "tx2", "tx3"),
    domain_id = c("PF00091", "PF00091", "PF00091"),
    domain_name = c("Tubulin", "Tubulin", "Tubulin"),
    domain_start = c(10L, 10L, 10L),
    domain_end = c(50L, 50L, 50L)
  )
  ref_ids <- c("tx1", "tx2", "tx3", "tx4", "tx5")

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  ref_ids, n_pairs = 5L)

  expect_true(is.data.frame(result))
  expected_cols <- c("domain_id", "domain_name", "n_ref_with_domain",
                     "n_affected", "n_pairs_total", "fisher_OR",
                     "fisher_p", "padj", "significant")
  for (col in expected_cols) {
    expect_true(col %in% names(result), info = paste("Missing:", col))
  }

  expect_equal(result$n_ref_with_domain, 3L)
  expect_equal(result$n_affected, 2L)
  expect_equal(result$n_pairs_total, 5L)
})

test_that("testDomainEnrichment returns empty tibble for empty domain_effects", {
  domain_effects <- tibble::tibble(
    domain_id = character(0),
    domain_name = character(0),
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    effect = character(0)
  )
  domain_annotations <- tibble::tibble(
    isoform_id = "tx1",
    domain_id = "PF00091",
    domain_name = "Tubulin",
    domain_start = 10L,
    domain_end = 50L
  )

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  "tx1", n_pairs = 10L)
  expect_equal(nrow(result), 0L)
})

test_that("testDomainEnrichment handles single domain", {
  domain_effects <- tibble::tibble(
    domain_id = "PF00091",
    domain_name = "Tubulin",
    reference_isoform_id = "tx1",
    comparator_isoform_id = "tx1c",
    effect = "lost"
  )
  domain_annotations <- tibble::tibble(
    isoform_id = c("tx1", "tx2"),
    domain_id = c("PF00091", "PF00091"),
    domain_name = c("Tubulin", "Tubulin"),
    domain_start = c(10L, 10L),
    domain_end = c(50L, 50L)
  )

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  c("tx1", "tx2", "tx3"), n_pairs = 3L)
  expect_equal(nrow(result), 1L)
  expect_equal(result$n_affected, 1L)
  expect_equal(result$n_ref_with_domain, 2L)
  expect_true(!is.na(result$fisher_p))
})

test_that("testDomainEnrichment filters by effects_to_test", {
  domain_effects <- tibble::tibble(
    domain_id = c("PF00091", "PF00091"),
    domain_name = c("Tubulin", "Tubulin"),
    reference_isoform_id = c("tx1", "tx2"),
    comparator_isoform_id = c("tx1c", "tx2c"),
    effect = c("lost", "gained_intact")  # only "lost" matches default
  )
  domain_annotations <- tibble::tibble(
    isoform_id = c("tx1", "tx2"),
    domain_id = c("PF00091", "PF00091"),
    domain_name = c("Tubulin", "Tubulin"),
    domain_start = c(10L, 10L),
    domain_end = c(50L, 50L)
  )

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  c("tx1", "tx2"), n_pairs = 2L)
  # Only "lost" should count as affected (default effects_to_test)
  expect_equal(result$n_affected, 1L)

  # Now test with all effects
  result2 <- testDomainEnrichment(domain_effects, domain_annotations,
                                   c("tx1", "tx2"), n_pairs = 2L,
                                   effects_to_test = c("lost",
                                     "gained_intact"))
  expect_equal(result2$n_affected, 2L)
})

test_that("testDomainEnrichment result is sorted by fisher_p", {
  domain_effects <- tibble::tibble(
    domain_id = c("PF00001", "PF00002", "PF00002"),
    domain_name = c("DomA", "DomB", "DomB"),
    reference_isoform_id = c("tx1", "tx2", "tx3"),
    comparator_isoform_id = c("tx1c", "tx2c", "tx3c"),
    effect = c("disrupted", "lost", "lost")
  )
  domain_annotations <- tibble::tibble(
    isoform_id = c("tx1", "tx2", "tx3", "tx4"),
    domain_id = c("PF00001", "PF00002", "PF00002", "PF00001"),
    domain_name = c("DomA", "DomB", "DomB", "DomA"),
    domain_start = c(10L, 10L, 10L, 10L),
    domain_end = c(50L, 50L, 50L, 50L)
  )

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  c("tx1", "tx2", "tx3", "tx4"),
                                  n_pairs = 4L)
  # Should be sorted by p-value
  expect_true(all(diff(result$fisher_p) >= 0, na.rm = TRUE))
})

test_that("testDomainEnrichment validates inputs", {
  expect_error(testDomainEnrichment("not_df", data.frame(), "tx1", 1L),
               "must be a data frame")

  de <- tibble::tibble(domain_id = "x", reference_isoform_id = "a",
                       comparator_isoform_id = "b", effect = "lost")
  da <- tibble::tibble(isoform_id = "a", domain_id = "x",
                       domain_name = "X")

  expect_error(testDomainEnrichment(de, da, character(0), 1L),
               "non-empty character vector")
  expect_error(testDomainEnrichment(de, da, "a", 0L),
               "positive integer")
})

test_that("testDomainEnrichment applies multiple testing correction", {
  domain_effects <- tibble::tibble(
    domain_id = c("PF00001", "PF00002"),
    domain_name = c("DomA", "DomB"),
    reference_isoform_id = c("tx1", "tx2"),
    comparator_isoform_id = c("tx1c", "tx2c"),
    effect = c("disrupted", "lost")
  )
  domain_annotations <- tibble::tibble(
    isoform_id = c("tx1", "tx2"),
    domain_id = c("PF00001", "PF00002"),
    domain_name = c("DomA", "DomB"),
    domain_start = c(10L, 10L),
    domain_end = c(50L, 50L)
  )

  result <- testDomainEnrichment(domain_effects, domain_annotations,
                                  c("tx1", "tx2"), n_pairs = 2L)
  # padj should be >= fisher_p (BH correction)
  expect_true(all(result$padj >= result$fisher_p, na.rm = TRUE))
})
