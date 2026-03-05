# Tests for buildProfiles and new computed columns

test_that("buildProfiles produces correct output for a small subset", {
  # Use first 5 pairs from test data
  pairs_subset <- test_pairs_raw[1:min(5, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  profiles <- buildProfiles(pairs_df, test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = TRUE, verbose = FALSE)

  expect_gt(nrow(profiles), 0)

  # Check required columns exist
  expected_cols <- c("gene_id", "reference_isoform_id",
                     "comparator_isoform_id", "n_events", "n_gain", "n_loss",
                     "bp_gain", "bp_loss", "net_bp", "junction_similarity",
                     "pct_ref_length", "detailed_events",
                     "reconstruction_status")
  expect_true(all(expected_cols %in% names(profiles)))
})

test_that("gain/loss metrics are non-negative", {
  pairs_subset <- test_pairs_raw[1:min(10, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  profiles <- buildProfiles(pairs_df, test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = FALSE, verbose = FALSE)

  expect_true(all(profiles$n_gain >= 0))
  expect_true(all(profiles$n_loss >= 0))
  expect_true(all(profiles$bp_gain >= 0))
  expect_true(all(profiles$bp_loss >= 0))
})

test_that("junction_similarity is between 0 and 1 (or NA)", {
  pairs_subset <- test_pairs_raw[1:min(10, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  profiles <- buildProfiles(pairs_df, test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = FALSE, verbose = FALSE)

  non_na <- profiles$junction_similarity[!is.na(profiles$junction_similarity)]
  if (length(non_na) > 0) {
    expect_true(all(non_na >= 0 & non_na <= 1))
  }
})

test_that("summarizeProfiles returns expected structure", {
  pairs_subset <- test_pairs_raw[1:min(5, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  profiles <- buildProfiles(pairs_df, test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = FALSE, verbose = FALSE)

  summary <- summarizeProfiles(profiles)
  expect_true(is.list(summary))
  expect_true("n_profiles" %in% names(summary))
  expect_true("n_genes" %in% names(summary))
  expect_true("event_counts" %in% names(summary))
})
