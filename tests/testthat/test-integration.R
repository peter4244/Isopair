# Integration tests: end-to-end workflows

test_that("full workflow: GTF -> profiles -> characterization -> comparison", {
  gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")

  # Parse
  structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
  expect_gt(nrow(structures), 0L)


  # Union exons
  ue <- buildUnionExons(structures, verbose = FALSE)
  expect_true(is.list(ue))
  expect_gt(nrow(ue$union_exons), 0L)
  expect_gt(nrow(ue$isoform_union_mapping), 0L)

  # Expression pairs
  expr_file <- system.file("extdata", "example_expression.csv",
    package = "Isopair")
  expr_df <- read.csv(expr_file)
  expr_mat <- as.matrix(expr_df[, -1])
  rownames(expr_mat) <- expr_df$isoform_id
  gene_map <- unique(as.data.frame(structures[, c("isoform_id", "gene_id")]))
  pairs <- generatePairsExpression(expr_mat, gene_map,
    colnames(expr_mat), method = "top_two")
  expect_gt(nrow(pairs), 0L)
  expect_true(all(c("gene_id", "reference_isoform_id",
    "comparator_isoform_id") %in% names(pairs)))

  # Build profiles
  profiles <- buildProfiles(pairs, structures,
    ue$union_exons, ue$isoform_union_mapping, verbose = FALSE)
  expect_gt(nrow(profiles), 0L)
  expect_true(all(profiles$reconstruction_status == "PASS"))

  # Summarize
  summ <- summarizeProfiles(profiles)
  expect_true(is.list(summ))
  expect_equal(summ$n_profiles, nrow(profiles))

  # Co-occurrence
  cooc <- testCooccurrence(profiles)
  expect_true(is.data.frame(cooc))
  expect_true("odds_ratio" %in% names(cooc))

  # CDS pipeline
  cds <- extractCdsAnnotations(gtf_path, structures$isoform_id,
    verbose = FALSE)
  expect_gt(nrow(cds), 0L)

  annotated <- annotateRegionTypes(ue$isoform_union_mapping, cds)
  expect_true("region_type" %in% names(annotated))

  er <- mapEventsToRegions(profiles, annotated, cds)
  expect_true(is.data.frame(er))

  # PTC
  ptc <- computePtcStatus(structures, cds)
  expect_true("has_ptc" %in% names(ptc))

  # Compare pair sets (self-comparison)
  result <- comparePairSets(profiles, profiles,
    min_overlap = 2L, verbose = FALSE)
  expect_true(is.list(result))
  expect_true(all(c("event_frequency", "event_count", "gain_loss",
    "cooccurrence") %in% names(result)))
})


test_that("DE pair generation workflow", {
  de_file <- system.file("extdata", "example_de.csv", package = "Isopair")
  de <- read.csv(de_file)

  expr_file <- system.file("extdata", "example_expression.csv",
    package = "Isopair")
  expr_df <- read.csv(expr_file)
  expr_mat <- as.matrix(expr_df[, -1])
  rownames(expr_mat) <- expr_df$isoform_id

  data(example_structures, envir = environment())
  gene_map <- unique(as.data.frame(
    example_structures[, c("isoform_id", "gene_id")]))

  dom <- identifyDominantIsoforms(expr_mat, gene_map, colnames(expr_mat))
  expect_true(is.data.frame(dom))
  expect_true("dominant_isoform_id" %in% names(dom))
  expect_gt(nrow(dom), 0L)

  pairs <- generatePairsDE(de, "isoform_id", "gene_id", "logFC",
    "adj_p_val", dom)
  expect_true(is.data.frame(pairs))
  expect_true(all(c("gene_id", "reference_isoform_id",
    "comparator_isoform_id", "source") %in% names(pairs)))
})


test_that("DU pair generation workflow", {
  du_file <- system.file("extdata", "example_du.csv", package = "Isopair")
  du <- read.csv(du_file)

  pairs <- generatePairsDU(du, "isoform_id", "gene_id", "dPSI",
    "adj_p_val", method = "extreme")
  expect_true(is.data.frame(pairs))
  expect_true(all(c("gene_id", "reference_isoform_id",
    "comparator_isoform_id", "source") %in% names(pairs)))
})
