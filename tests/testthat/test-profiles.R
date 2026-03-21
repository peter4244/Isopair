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

# ==========================================================================
# Checkpoint support
# ==========================================================================

test_that("buildProfiles: checkpoint_dir=NULL is backward compatible", {
  pairs_subset <- test_pairs_raw[1:min(3, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )
  profiles <- buildProfiles(pairs_df, test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = FALSE, verbose = FALSE,
                            checkpoint_dir = NULL)
  expect_gt(nrow(profiles), 0)
  expect_true("gene_id" %in% names(profiles))
})

test_that("buildProfiles: checkpoint creates and cleans up files", {
  ckpt_dir <- file.path(tempdir(), "isopair_ckpt_test_create")
  on.exit(unlink(ckpt_dir, recursive = TRUE), add = TRUE)

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
                            verify = FALSE, verbose = FALSE,
                            checkpoint_dir = ckpt_dir,
                            checkpoint_interval = 2L)
  expect_gt(nrow(profiles), 0)
  # After successful completion, checkpoint dir should be cleaned up
  expect_false(dir.exists(ckpt_dir))
})

test_that("buildProfiles: checkpoint resume produces same result", {
  ckpt_dir <- file.path(tempdir(), "isopair_ckpt_test_resume")
  on.exit(unlink(ckpt_dir, recursive = TRUE), add = TRUE)

  pairs_subset <- test_pairs_raw[1:min(5, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  # Run without checkpoint (baseline)
  profiles_baseline <- buildProfiles(pairs_df, test_structures,
                                      test_union_exons, test_isoform_union_mapping,
                                      verify = FALSE, verbose = FALSE)

  # Simulate a partial run by creating a checkpoint manually
  dir.create(ckpt_dir, recursive = TRUE)
  partial_profiles <- list()
  # Process first 2 pairs manually to create a checkpoint
  partial <- buildProfiles(pairs_df[1:2, ], test_structures,
                            test_union_exons, test_isoform_union_mapping,
                            verify = FALSE, verbose = FALSE)
  for (i in seq_len(nrow(partial))) {
    partial_profiles[[i]] <- partial[i, ]
  }

  ckpt_data <- list(
    last_index = 2L,
    n_pairs = nrow(pairs_df),
    first_pair = paste(pairs_df$reference_isoform_id[1],
                        pairs_df$comparator_isoform_id[1]),
    profiles = partial_profiles
  )
  saveRDS(ckpt_data, file.path(ckpt_dir, "checkpoint_000002.rds"))

  # Resume from checkpoint
  profiles_resumed <- buildProfiles(pairs_df, test_structures,
                                     test_union_exons, test_isoform_union_mapping,
                                     verify = FALSE, verbose = FALSE,
                                     checkpoint_dir = ckpt_dir,
                                     checkpoint_interval = 100L)

  expect_equal(nrow(profiles_resumed), nrow(profiles_baseline))
  expect_equal(profiles_resumed$gene_id, profiles_baseline$gene_id)
  expect_equal(profiles_resumed$n_events, profiles_baseline$n_events)
})

test_that("buildProfiles: checkpoint input mismatch warns and starts fresh", {
  ckpt_dir <- file.path(tempdir(), "isopair_ckpt_test_mismatch")
  on.exit(unlink(ckpt_dir, recursive = TRUE), add = TRUE)
  dir.create(ckpt_dir, recursive = TRUE)

  pairs_subset <- test_pairs_raw[1:min(5, nrow(test_pairs_raw)), ]
  pairs_df <- data.frame(
    gene_id = pairs_subset$gene_id,
    reference_isoform_id = pairs_subset$isoform_A,
    comparator_isoform_id = pairs_subset$isoform_B,
    strand = pairs_subset$strand,
    stringsAsFactors = FALSE
  )

  # Create a checkpoint with wrong n_pairs
  ckpt_data <- list(
    last_index = 2L,
    n_pairs = 999L,  # Wrong!
    first_pair = "WRONG WRONG",
    profiles = list()
  )
  saveRDS(ckpt_data, file.path(ckpt_dir, "checkpoint_000002.rds"))

  expect_warning(
    profiles <- buildProfiles(pairs_df, test_structures,
                               test_union_exons, test_isoform_union_mapping,
                               verify = FALSE, verbose = FALSE,
                               checkpoint_dir = ckpt_dir,
                               checkpoint_interval = 100L),
    "mismatch"
  )
  expect_gt(nrow(profiles), 0)
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
