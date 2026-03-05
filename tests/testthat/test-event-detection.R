# Tests for event detection and reconstruction (126 pairs)
#
# Ported from scripts/tests/run_tests.R
# Each pair: detect events -> reconstruct -> verify PASS

test_that("all 126 test pairs reconstruct correctly", {
  # Process each pair
  n_pairs <- nrow(test_pairs_raw)
  results <- character(n_pairs)
  reasons <- character(n_pairs)

  for (i in seq_len(n_pairs)) {
    pair <- test_pairs_raw[i, ]

    ref_exons <- test_exons[test_exons$gene_id == pair$gene_id &
                            test_exons$isoform_id == pair$isoform_A, ]
    comp_exons <- test_exons[test_exons$gene_id == pair$gene_id &
                             test_exons$isoform_id == pair$isoform_B, ]

    if (nrow(ref_exons) == 0 || nrow(comp_exons) == 0) {
      results[i] <- "ERROR"
      reasons[i] <- "Missing exons"
      next
    }

    gene_strand <- unique(ref_exons$strand)[1]

    result <- tryCatch({
      events <- detectEvents(
        ref_exons, comp_exons,
        pair$gene_id, pair$isoform_A, pair$isoform_B, gene_strand
      )

      if (nrow(events) == 0) {
        vr <- verifyReconstruction(ref_exons, comp_exons, gene_strand)
        list(status = vr$status, reason = vr$reason)
      } else {
        comp_for_recon <- comp_exons[, c("chr", "exon_start", "exon_end",
                                          "strand", "gene_id",
                                          "transcript_id")]
        reconstructed <- reconstructDominant(comp_for_recon, events)

        if (nrow(reconstructed) == 0) {
          list(status = "FAIL", reason = "Reconstruction produced 0 exons")
        } else {
          vr <- verifyReconstruction(ref_exons, reconstructed, gene_strand)
          list(status = vr$status, reason = vr$reason)
        }
      }
    }, error = function(e) {
      list(status = "ERROR", reason = paste0("Exception: ", e$message))
    })

    results[i] <- result$status
    reasons[i] <- result$reason
  }

  # Report failures
  failures <- which(results != "PASS")
  if (length(failures) > 0) {
    fail_msgs <- vapply(failures, function(i) {
      sprintf("  [%s] %s: %s vs %s -- %s",
              results[i], test_pairs_raw$gene_id[i],
              test_pairs_raw$isoform_A[i], test_pairs_raw$isoform_B[i],
              reasons[i])
    }, character(1))
    fail(paste(c(sprintf("%d/%d pairs failed:", length(failures), n_pairs),
                 fail_msgs), collapse = "\n"))
  }

  expect_equal(sum(results == "PASS"), n_pairs)
})

# Individual pair tests for key event types
test_that("basic SE detection works", {
  pair <- test_pairs_raw[test_pairs_raw$gene_id == "TEST_SE_Basic", ]
  if (nrow(pair) == 0) skip("TEST_SE_Basic not in pairs")

  ref_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                          test_exons$isoform_id == pair$isoform_A[1], ]
  comp_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                           test_exons$isoform_id == pair$isoform_B[1], ]

  events <- detectEvents(ref_exons, comp_exons,
                         pair$gene_id[1], pair$isoform_A[1],
                         pair$isoform_B[1], unique(ref_exons$strand)[1])

  expect_true("SE" %in% events$event_type)
})

test_that("basic IR detection works", {
  pair <- test_pairs_raw[test_pairs_raw$gene_id == "TEST_IR_Basic", ]
  if (nrow(pair) == 0) skip("TEST_IR_Basic not in pairs")

  ref_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                          test_exons$isoform_id == pair$isoform_A[1], ]
  comp_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                           test_exons$isoform_id == pair$isoform_B[1], ]

  events <- detectEvents(ref_exons, comp_exons,
                         pair$gene_id[1], pair$isoform_A[1],
                         pair$isoform_B[1], unique(ref_exons$strand)[1])

  ir_types <- c("IR", "IR_diff_5", "IR_diff_3", "IR_diff_5_3")
  expect_true(any(events$event_type %in% ir_types))
})

test_that("multi-event detection works", {
  pair <- test_pairs_raw[test_pairs_raw$gene_id == "TEST_MultiEvent", ]
  if (nrow(pair) == 0) skip("TEST_MultiEvent not in pairs")

  ref_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                          test_exons$isoform_id == pair$isoform_A[1], ]
  comp_exons <- test_exons[test_exons$gene_id == pair$gene_id[1] &
                           test_exons$isoform_id == pair$isoform_B[1], ]

  events <- detectEvents(ref_exons, comp_exons,
                         pair$gene_id[1], pair$isoform_A[1],
                         pair$isoform_B[1], unique(ref_exons$strand)[1])

  expect_gt(nrow(events), 1)
})
