# Tests for pair-generation.R
#
# Inline test fixtures (no file dependencies) for speed.

# ============================================================================
# Test fixtures
# ============================================================================

# 3 genes, 8 isoforms
.make_test_gene_map <- function() {
  data.frame(
    isoform_id = c("TX_A1", "TX_A2", "TX_A3",
                    "TX_B1", "TX_B2",
                    "TX_C1"),
    gene_id = c("GENE_A", "GENE_A", "GENE_A",
                "GENE_B", "GENE_B",
                "GENE_C"),
    stringsAsFactors = FALSE
  )
}

# Expression matrix: A1 dominant in GENE_A, B1 dominant in GENE_B, C1 is sole
.make_test_expression <- function() {
  mat <- matrix(
    c(
      # S1    S2    S3    S4    S5    S6
      200,  220,  190,  210,  205,  215,   # TX_A1 (dominant)
       30,   25,   35,   28,   32,   27,   # TX_A2
        5,    8,    3,    6,    7,    4,   # TX_A3
      150,  160,  140,  155,  145,  158,   # TX_B1 (dominant)
       40,   35,   45,   38,   42,   37,   # TX_B2
       80,   85,   75,   82,   78,   83    # TX_C1 (sole)
    ),
    nrow = 6, ncol = 6, byrow = TRUE,
    dimnames = list(
      c("TX_A1", "TX_A2", "TX_A3", "TX_B1", "TX_B2", "TX_C1"),
      paste0("S", 1:6)
    )
  )
  mat
}


# ============================================================================
# identifyDominantIsoforms
# ============================================================================

test_that("identifyDominantIsoforms selects correct dominant", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()
  samples <- paste0("S", 1:6)

  dom <- identifyDominantIsoforms(expr_mat, gene_map, samples)

  expect_s3_class(dom, "tbl_df")
  # GENE_A dominant should be TX_A1

  dom_a <- dom[dom$gene_id == "GENE_A", ]
  expect_equal(dom_a$dominant_isoform_id, "TX_A1")
  # GENE_B dominant should be TX_B1
  dom_b <- dom[dom$gene_id == "GENE_B", ]
  expect_equal(dom_b$dominant_isoform_id, "TX_B1")
  # GENE_C has only one isoform
  dom_c <- dom[dom$gene_id == "GENE_C", ]
  expect_equal(dom_c$dominant_isoform_id, "TX_C1")
  expect_equal(dom_c$proportion, 1.0)
})

test_that("identifyDominantIsoforms respects threshold", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()
  samples <- paste0("S", 1:6)

  # TX_A1 proportion ~ 200/235 ~ 0.85 — should pass 0.5

  dom_50 <- identifyDominantIsoforms(expr_mat, gene_map, samples,
                                     threshold = 0.5)
  expect_true("GENE_A" %in% dom_50$gene_id)

  # With threshold 0.99 — only GENE_C (sole isoform, proportion=1.0) passes
  dom_99 <- identifyDominantIsoforms(expr_mat, gene_map, samples,
                                     threshold = 0.99)
  expect_equal(nrow(dom_99), 1L)
  expect_equal(dom_99$gene_id, "GENE_C")
})

test_that("identifyDominantIsoforms handles all-zero gene", {
  gene_map <- data.frame(
    isoform_id = c("TX1", "TX2"),
    gene_id = c("G1", "G1"),
    stringsAsFactors = FALSE
  )
  mat <- matrix(0, nrow = 2, ncol = 2,
                dimnames = list(c("TX1", "TX2"), c("S1", "S2")))

  dom <- identifyDominantIsoforms(mat, gene_map, c("S1", "S2"))
  # No dominant when all zero
  expect_equal(nrow(dom), 0L)
})

test_that("identifyDominantIsoforms handles single-sample", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()

  dom <- identifyDominantIsoforms(expr_mat, gene_map, "S1")
  expect_true(nrow(dom) > 0)
  expect_equal(dom$dominant_isoform_id[dom$gene_id == "GENE_A"], "TX_A1")
})

test_that("identifyDominantIsoforms validates inputs", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()

  expect_error(identifyDominantIsoforms("not_a_matrix", gene_map, "S1"),
               "must be a numeric matrix")
  expect_error(identifyDominantIsoforms(expr_mat, gene_map, "NONEXISTENT"),
               "Samples not found")
})


# ============================================================================
# generatePairsExpression
# ============================================================================

test_that("generatePairsExpression top_two creates 1 pair per gene", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()
  samples <- paste0("S", 1:6)

  pairs <- generatePairsExpression(expr_mat, gene_map, samples,
                                   method = "top_two")

  expect_s3_class(pairs, "tbl_df")
  # GENE_A: TX_A1 vs TX_A2 (top two by expression)
  pair_a <- pairs[pairs$gene_id == "GENE_A", ]
  expect_equal(nrow(pair_a), 1L)
  expect_equal(pair_a$reference_isoform_id, "TX_A1")
  expect_equal(pair_a$comparator_isoform_id, "TX_A2")

  # GENE_C has only 1 isoform — should NOT generate a pair
  expect_false("GENE_C" %in% pairs$gene_id)

  # Check metadata columns
  expect_true(all(pairs$source == "expression"))
  expect_true(all(is.na(pairs$direction)))
})

test_that("generatePairsExpression dominant_vs_rest creates N-1 pairs", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()
  samples <- paste0("S", 1:6)

  pairs <- generatePairsExpression(expr_mat, gene_map, samples,
                                   method = "dominant_vs_rest",
                                   dominant_threshold = 0.5)

  # GENE_A: TX_A1 (dominant) vs TX_A2 and TX_A3 = 2 pairs
  pairs_a <- pairs[pairs$gene_id == "GENE_A", ]
  expect_equal(nrow(pairs_a), 2L)
  expect_true(all(pairs_a$reference_isoform_id == "TX_A1"))
  expect_true(all(c("TX_A2", "TX_A3") %in% pairs_a$comparator_isoform_id))

  # GENE_B: TX_B1 vs TX_B2 = 1 pair
  pairs_b <- pairs[pairs$gene_id == "GENE_B", ]
  expect_equal(nrow(pairs_b), 1L)
})

test_that("generatePairsExpression min_expression filters comparators", {
  gene_map <- .make_test_gene_map()
  expr_mat <- .make_test_expression()
  samples <- paste0("S", 1:6)

  # TX_A3 mean ~ 5.5, set min_expression = 10 to exclude it
  pairs <- generatePairsExpression(expr_mat, gene_map, samples,
                                   method = "dominant_vs_rest",
                                   dominant_threshold = 0.5,
                                   min_expression = 10)

  pairs_a <- pairs[pairs$gene_id == "GENE_A", ]
  expect_equal(nrow(pairs_a), 1L)  # Only TX_A2 left
  expect_equal(pairs_a$comparator_isoform_id, "TX_A2")
})


# ============================================================================
# generatePairsDE
# ============================================================================

test_that("generatePairsDE pairs significant isoforms with dominant", {
  de_results <- data.frame(
    tx = c("TX_A1", "TX_A2", "TX_A3", "TX_B1", "TX_B2"),
    gn = c("GENE_A", "GENE_A", "GENE_A", "GENE_B", "GENE_B"),
    lfc = c(0.5, 2.5, -1.8, 0.1, 3.0),
    pval = c(0.8, 0.001, 0.01, 0.9, 0.002),
    stringsAsFactors = FALSE
  )

  dom <- data.frame(
    gene_id = c("GENE_A", "GENE_B"),
    dominant_isoform_id = c("TX_A1", "TX_B1"),
    mean_expression = c(200, 150),
    proportion = c(0.85, 0.79),
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDE(de_results, "tx", "gn", "lfc", "pval",
                           dominant_isoforms = dom)

  expect_s3_class(pairs, "tbl_df")
  # TX_A2 (up, lfc=2.5, p=0.001) paired with TX_A1
  # TX_A3 (down, lfc=-1.8, p=0.01) paired with TX_A1
  # TX_B2 (up, lfc=3.0, p=0.002) paired with TX_B1
  expect_equal(nrow(pairs), 3L)
  expect_true(all(pairs$source == "DE"))

  # Check direction
  pair_a2 <- pairs[pairs$comparator_isoform_id == "TX_A2", ]
  expect_equal(pair_a2$direction, "up")
  pair_a3 <- pairs[pairs$comparator_isoform_id == "TX_A3", ]
  expect_equal(pair_a3$direction, "down")
})

test_that("generatePairsDE skips self-pairs", {
  # Dominant isoform is itself significant
  de_results <- data.frame(
    tx = c("TX_A1", "TX_A2"),
    gn = c("GENE_A", "GENE_A"),
    lfc = c(2.0, 3.0),
    pval = c(0.01, 0.001),
    stringsAsFactors = FALSE
  )

  dom <- data.frame(
    gene_id = "GENE_A",
    dominant_isoform_id = "TX_A1",
    mean_expression = 200,
    proportion = 0.85,
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDE(de_results, "tx", "gn", "lfc", "pval",
                           dominant_isoforms = dom)

  # TX_A1 is dominant AND significant — should not pair with itself
  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$comparator_isoform_id, "TX_A2")
})

test_that("generatePairsDE logfc_threshold filters", {
  de_results <- data.frame(
    tx = c("TX_A2", "TX_A3"),
    gn = c("GENE_A", "GENE_A"),
    lfc = c(0.5, 2.0),
    pval = c(0.01, 0.01),
    stringsAsFactors = FALSE
  )

  dom <- data.frame(
    gene_id = "GENE_A",
    dominant_isoform_id = "TX_A1",
    mean_expression = 200,
    proportion = 0.85,
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDE(de_results, "tx", "gn", "lfc", "pval",
                           dominant_isoforms = dom,
                           logfc_threshold = 1.0)

  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$comparator_isoform_id, "TX_A3")
})

test_that("generatePairsDE validates column names", {
  de_results <- data.frame(x = 1, y = 2)
  dom <- data.frame(gene_id = "G", dominant_isoform_id = "T",
                    mean_expression = 1, proportion = 1)

  expect_error(
    generatePairsDE(de_results, "tx", "gn", "lfc", "pval", dom),
    "Columns not found"
  )
})


# ============================================================================
# generatePairsDU
# ============================================================================

test_that("generatePairsDU extreme mode pairs positive vs negative dPSI", {
  du_results <- data.frame(
    tx = c("TX_A1", "TX_A2", "TX_A3", "TX_B1", "TX_B2"),
    gn = c("GENE_A", "GENE_A", "GENE_A", "GENE_B", "GENE_B"),
    dpsi = c(0.3, -0.25, 0.1, 0.4, -0.35),
    pval = c(0.001, 0.01, 0.8, 0.002, 0.003),
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDU(du_results, "tx", "gn", "dpsi", "pval",
                           method = "extreme")

  expect_s3_class(pairs, "tbl_df")
  # GENE_A: TX_A1 (most positive sig) vs TX_A2 (most negative sig)
  # TX_A3 is not significant (p=0.8)
  pair_a <- pairs[pairs$gene_id == "GENE_A", ]
  expect_equal(nrow(pair_a), 1L)
  expect_equal(pair_a$reference_isoform_id, "TX_A1")
  expect_equal(pair_a$comparator_isoform_id, "TX_A2")

  expect_true(all(pairs$source == "DU"))
})

test_that("generatePairsDU extreme skips genes without both directions", {
  du_results <- data.frame(
    tx = c("TX_A1", "TX_A2"),
    gn = c("GENE_A", "GENE_A"),
    dpsi = c(0.3, 0.2),  # Both positive
    pval = c(0.01, 0.01),
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDU(du_results, "tx", "gn", "dpsi", "pval",
                           method = "extreme")
  expect_equal(nrow(pairs), 0L)
})

test_that("generatePairsDU dominant_vs_sig mode works", {
  du_results <- data.frame(
    tx = c("TX_A1", "TX_A2", "TX_A3"),
    gn = c("GENE_A", "GENE_A", "GENE_A"),
    dpsi = c(0.1, -0.3, 0.2),
    pval = c(0.8, 0.001, 0.01),
    stringsAsFactors = FALSE
  )

  dom <- data.frame(
    gene_id = "GENE_A",
    dominant_isoform_id = "TX_A1",
    mean_expression = 200,
    proportion = 0.85,
    stringsAsFactors = FALSE
  )

  pairs <- generatePairsDU(du_results, "tx", "gn", "dpsi", "pval",
                           dominant_isoforms = dom,
                           method = "dominant_vs_sig")

  expect_equal(nrow(pairs), 1L)
  expect_equal(pairs$reference_isoform_id, "TX_A1")
  # TX_A2 is most significant (p=0.001)
  expect_equal(pairs$comparator_isoform_id, "TX_A2")
  expect_equal(pairs$direction, "down")  # dpsi = -0.3
})

test_that("generatePairsDU dominant_vs_sig requires dominant_isoforms", {
  du <- data.frame(tx = "TX1", gn = "G1", dpsi = 0.3, pval = 0.01)
  expect_error(
    generatePairsDU(du, "tx", "gn", "dpsi", "pval",
                    method = "dominant_vs_sig"),
    "dominant_isoforms is required"
  )
})


# ============================================================================
# End-to-end: pair generation → buildProfiles (strand auto-lookup)
# ============================================================================

test_that("generatePairsExpression → buildProfiles works with strand auto-lookup", {
  # Use the existing test fixtures from helper-setup.R
  # Pick 2 genes with at least 2 isoforms
  genes_with_2 <- test_structures |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::filter(dplyr::n() >= 2L) |>
    dplyr::ungroup()

  if (nrow(genes_with_2) < 2) skip("Need at least 2 multi-isoform genes")

  # Use first 2 genes
  selected_genes <- unique(genes_with_2$gene_id)[1:2]
  sub_structures <- dplyr::filter(genes_with_2,
                                  .data$gene_id %in% selected_genes)

  # Build a mock expression matrix
  iso_ids <- sub_structures$isoform_id
  n_iso <- length(iso_ids)
  set.seed(123)
  mock_expr <- matrix(
    runif(n_iso * 4, 10, 200), nrow = n_iso, ncol = 4,
    dimnames = list(iso_ids, paste0("S", 1:4))
  )

  gene_map <- data.frame(
    isoform_id = sub_structures$isoform_id,
    gene_id = sub_structures$gene_id,
    stringsAsFactors = FALSE
  )

  # Generate pairs (no strand column in output)
  pairs <- generatePairsExpression(mock_expr, gene_map, paste0("S", 1:4),
                                   method = "top_two")

  expect_true(nrow(pairs) > 0)
  expect_false("strand" %in% names(pairs))

  # Build union exons for selected genes
  sub_ue <- buildUnionExons(sub_structures, verbose = FALSE)

  # buildProfiles should auto-lookup strand
  profiles <- buildProfiles(
    pairs, sub_structures,
    sub_ue$union_exons, sub_ue$isoform_union_mapping,
    verify = TRUE, verbose = FALSE
  )

  expect_true(nrow(profiles) > 0)
  expect_true(all(!is.na(profiles$reconstruction_status)))
})
