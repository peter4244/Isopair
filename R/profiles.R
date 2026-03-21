#' Build Splicing Choice Profiles
#'
#' For each isoform pair, detects structural events, optionally reconstructs
#' and verifies the reference isoform, and builds a profile row summarizing
#' the structural differences.
#'
#' @param pairs A data frame with columns: gene_id, reference_isoform_id,
#'   comparator_isoform_id, and optionally strand. If strand is absent, it
#'   is looked up from `structures` using the reference isoform ID.
#' @param structures A tibble from [parseIsoformStructures()] with nested
#'   exon_starts and exon_ends.
#' @param union_exons Union exon tibble from [buildUnionExons()]$union_exons.
#' @param isoform_union_mapping Mapping tibble from
#'   [buildUnionExons()]$isoform_union_mapping.
#' @param verify Logical; if TRUE (default), run reconstruction verification.
#' @param verbose Logical; if TRUE, print progress messages.
#' @param checkpoint_dir Optional directory path. If provided, saves partial
#'   results every \code{checkpoint_interval} pairs and can resume from the
#'   last checkpoint on re-invocation. After successful completion, checkpoint
#'   files are removed.
#' @param checkpoint_interval Integer; save a checkpoint every this many pairs
#'   (default 1000). Only used when \code{checkpoint_dir} is not NULL.
#' @return A tibble of splicing choice profiles with one row per pair.
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' expr_file <- system.file("extdata", "example_expression.csv",
#'   package = "Isopair")
#' expr_df <- read.csv(expr_file)
#' expr_mat <- as.matrix(expr_df[, -1])
#' rownames(expr_mat) <- expr_df$isoform_id
#' gene_map <- unique(structures[, c("isoform_id", "gene_id")])
#' pairs <- generatePairsExpression(expr_mat, gene_map,
#'   colnames(expr_mat), method = "top_two")
#' profiles <- buildProfiles(pairs, structures, ue$union_exons,
#'   ue$isoform_union_mapping, verbose = FALSE)
#' @export
#' @importFrom dplyr filter select mutate left_join group_by summarise n
#'   bind_rows n_distinct arrange
#' @importFrom tibble tibble
#' @importFrom rlang .data
buildProfiles <- function(pairs, structures, union_exons,
                          isoform_union_mapping,
                          verify = TRUE, verbose = TRUE,
                          checkpoint_dir = NULL,
                          checkpoint_interval = 1000L) {

  if (verbose) message(sprintf("Building profiles for %d pairs...", nrow(pairs)))

  # Auto-lookup strand if not in pairs
  if (!"strand" %in% names(pairs)) {
    strand_idx <- match(pairs$reference_isoform_id, structures$isoform_id)
    if (any(is.na(strand_idx))) {
      missing <- pairs$reference_isoform_id[is.na(strand_idx)]
      stop(sprintf(
        "Cannot determine strand: reference isoform '%s' not found in structures.",
        missing[1L]
      ))
    }
    pairs$strand <- structures$strand[strand_idx]
    if (verbose) message("  Strand not in pairs; looked up from structures.")
  }

  # Expand structures to per-exon form
  expanded <- .expandStructures(structures)

  profiles <- list()
  start_idx <- 1L

  # Checkpoint: resume from previous run if checkpoint_dir is set
  use_checkpoint <- !is.null(checkpoint_dir)
  if (use_checkpoint) {
    if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)

    input_fingerprint <- paste(pairs$reference_isoform_id[1L],
                                pairs$comparator_isoform_id[1L])

    # Look for existing checkpoints
    ckpt_files <- sort(list.files(checkpoint_dir, pattern = "^checkpoint_.*\\.rds$",
                                   full.names = TRUE))
    if (length(ckpt_files) > 0L) {
      latest <- readRDS(ckpt_files[length(ckpt_files)])
      if (latest$n_pairs == nrow(pairs) &&
          latest$first_pair == input_fingerprint) {
        profiles <- latest$profiles
        start_idx <- latest$last_index + 1L
        if (verbose) {
          message(sprintf("  Resuming from checkpoint: %d/%d pairs already done",
                          latest$last_index, nrow(pairs)))
        }
      } else {
        warning("Checkpoint input mismatch; starting fresh.")
        # Clean up stale checkpoints
        file.remove(ckpt_files)
      }
    }
  }

  if (start_idx > nrow(pairs)) {
    # All pairs already processed from checkpoint
    result <- dplyr::bind_rows(profiles)
    if (use_checkpoint) .cleanCheckpoints(checkpoint_dir)
    if (verbose) {
      message(sprintf("  Total profiles: %d covering %d genes",
                      nrow(result), dplyr::n_distinct(result$gene_id)))
    }
    return(result)
  }

  for (idx in seq(start_idx, nrow(pairs))) {
    pair <- pairs[idx, ]
    gene <- pair$gene_id
    ref_id <- pair$reference_isoform_id
    comp_id <- pair$comparator_isoform_id
    gene_strand <- pair$strand

    # Get exon structures
    ref_structure <- expanded[expanded$isoform_id == ref_id, ]
    comp_structure <- expanded[expanded$isoform_id == comp_id, ]

    if (nrow(ref_structure) == 0L || nrow(comp_structure) == 0L) next

    ref_structure <- ref_structure[order(ref_structure$exon_number), ]
    comp_structure <- comp_structure[order(comp_structure$exon_number), ]

    # Prepare exons for detection (needs chr, transcript_id columns)
    ref_for_detect <- ref_structure
    names(ref_for_detect)[names(ref_for_detect) == "chr"] <- "chr"
    ref_for_detect$transcript_id <- ref_for_detect$isoform_id

    comp_for_detect <- comp_structure
    names(comp_for_detect)[names(comp_for_detect) == "chr"] <- "chr"
    comp_for_detect$transcript_id <- comp_for_detect$isoform_id

    # Detect events
    events <- detectEvents(
      ref_for_detect, comp_for_detect,
      gene, ref_id, comp_id, gene_strand
    )

    # Reconstruction verification
    recon_status <- NA_character_
    recon_reason <- NA_character_

    if (verify) {
      recon_result <- tryCatch({
        comp_for_recon <- comp_structure[, c("chr", "exon_start", "exon_end",
                                             "strand", "gene_id")]
        comp_for_recon$transcript_id <- comp_id

        if (nrow(events) == 0L) {
          ref_for_verify <- ref_structure[, c("chr", "exon_start", "exon_end",
                                              "strand", "gene_id")]
          ref_for_verify$transcript_id <- ref_id
          vr <- verifyReconstruction(ref_for_verify, comp_for_recon,
                                     gene_strand)
          list(status = vr$status, reason = vr$reason)
        } else {
          reconstructed <- reconstructReference(comp_for_recon, events)
          if (nrow(reconstructed) == 0L) {
            list(status = "FAIL",
                 reason = "Reconstruction produced 0 exons")
          } else {
            vr <- verifyReconstruction(ref_structure, reconstructed,
                                       gene_strand)
            list(status = vr$status, reason = vr$reason)
          }
        }
      }, error = function(e) {
        list(status = "ERROR", reason = paste0("Exception: ", e$message))
      })

      recon_status <- recon_result$status
      recon_reason <- recon_result$reason
    }

    # Event counts
    n_events <- nrow(events)
    if (n_events > 0L) {
      n_a5ss <- sum(events$event_type == "A5SS")
      n_a3ss <- sum(events$event_type == "A3SS")
      n_partial_ir <- sum(events$event_type %in%
                            c("Partial_IR_5", "Partial_IR_3"))
      n_ir <- sum(events$event_type == "IR")
      n_se <- sum(events$event_type == "SE")
      n_missing_internal <- sum(events$event_type == "Missing_Internal")
      n_ir_diff <- sum(grepl("^IR_diff", events$event_type))
      n_alt_tss <- sum(events$event_type == "Alt_TSS")
      n_alt_tes <- sum(events$event_type == "Alt_TES")
      tss_changed <- n_alt_tss > 0L
      tes_changed <- n_alt_tes > 0L
    } else {
      n_a5ss <- 0L; n_a3ss <- 0L; n_partial_ir <- 0L
      n_ir <- 0L; n_se <- 0L; n_missing_internal <- 0L
      n_ir_diff <- 0L; n_alt_tss <- 0L; n_alt_tes <- 0L
      tss_changed <- FALSE; tes_changed <- FALSE
    }

    # Gain/loss computation
    gain_loss <- .computeGainLoss(events)

    # Complexity metrics from compact structures
    ref_compact <- structures[structures$isoform_id == ref_id, ]
    comp_compact <- structures[structures$isoform_id == comp_id, ]

    n_exons_ref <- if (nrow(ref_compact) > 0L) ref_compact$n_exons[1L] else NA_integer_
    n_junctions_ref <- if (nrow(ref_compact) > 0L) ref_compact$n_junctions[1L] else NA_integer_
    length_ref <- if (nrow(ref_compact) > 0L) {
      ref_compact$tx_end[1L] - ref_compact$tx_start[1L] + 1L
    } else NA_integer_
    n_exons_comp <- if (nrow(comp_compact) > 0L) comp_compact$n_exons[1L] else NA_integer_
    n_junctions_comp <- if (nrow(comp_compact) > 0L) comp_compact$n_junctions[1L] else NA_integer_
    length_comp <- if (nrow(comp_compact) > 0L) {
      comp_compact$tx_end[1L] - comp_compact$tx_start[1L] + 1L
    } else NA_integer_

    # UE comparison
    ue_comp <- .computeUeComparison(ref_id, comp_id, isoform_union_mapping)

    # Junction similarity
    ref_junctions_all <- .computeJunctions(
      .orderExonsBiological(ref_structure, gene_strand))
    comp_junctions_all <- .computeJunctions(
      .orderExonsBiological(comp_structure, gene_strand))
    jxn_stats <- .computeSharedJunctions(ref_junctions_all, comp_junctions_all)

    # Percent reference length
    pct_ref_length <- if (!is.na(length_ref) && !is.na(length_comp) &&
                          length_ref > 0L) {
      100.0 * length_comp / length_ref
    } else {
      NA_real_
    }

    profiles[[idx]] <- tibble::tibble(
      gene_id = gene,
      reference_isoform_id = ref_id,
      comparator_isoform_id = comp_id,

      # UE comparison
      n_union_exons_total = ue_comp$n_union_exons_total,
      n_exons_shared = ue_comp$n_exons_shared,
      n_exons_ref_only = ue_comp$n_exons_ref_only,
      n_exons_comp_only = ue_comp$n_exons_comp_only,

      # Complexity metrics
      n_exons_ref = n_exons_ref,
      n_junctions_ref = n_junctions_ref,
      length_ref = length_ref,
      n_exons_comp = n_exons_comp,
      n_junctions_comp = n_junctions_comp,
      length_comp = length_comp,

      # TSS/TES changes
      tss_changed = tss_changed,
      tes_changed = tes_changed,

      # Event counts
      n_a5ss = n_a5ss,
      n_a3ss = n_a3ss,
      n_partial_ir = n_partial_ir,
      n_ir = n_ir,
      n_se = n_se,
      n_missing_internal = n_missing_internal,
      n_ir_diff = n_ir_diff,
      n_alt_tss = n_alt_tss,
      n_alt_tes = n_alt_tes,
      n_events = n_events,

      # Gain/loss metrics
      n_gain = gain_loss$n_gain,
      n_loss = gain_loss$n_loss,
      bp_gain = gain_loss$bp_gain,
      bp_loss = gain_loss$bp_loss,
      net_bp = gain_loss$net_bp,

      # Junction similarity (Jaccard)
      junction_similarity = jxn_stats$junction_similarity,

      # Relative length
      pct_ref_length = pct_ref_length,

      # Detailed events (nested)
      detailed_events = list(events),

      # Reconstruction
      reconstruction_status = recon_status,
      reconstruction_reason = recon_reason
    )

    if (verbose && idx %% 100L == 0L) {
      message(sprintf("  Processed %d/%d pairs", idx, nrow(pairs)))
    }

    # Save checkpoint
    if (use_checkpoint && idx %% checkpoint_interval == 0L) {
      ckpt_data <- list(
        last_index = idx,
        n_pairs = nrow(pairs),
        first_pair = input_fingerprint,
        profiles = profiles
      )
      tmp_path <- file.path(checkpoint_dir,
                             sprintf(".tmp_checkpoint_%06d.rds", idx))
      final_path <- file.path(checkpoint_dir,
                                sprintf("checkpoint_%06d.rds", idx))
      saveRDS(ckpt_data, tmp_path)
      file.rename(tmp_path, final_path)
      if (verbose) message(sprintf("  Checkpoint saved: %d/%d", idx, nrow(pairs)))
    }
  }

  result <- dplyr::bind_rows(profiles)

  # Clean up checkpoints on successful completion
  if (use_checkpoint) .cleanCheckpoints(checkpoint_dir)

  if (verbose) {
    message(sprintf("  Total profiles: %d covering %d genes",
                    nrow(result), dplyr::n_distinct(result$gene_id)))
  }

  result
}


#' Remove checkpoint files and directory
#' @keywords internal
.cleanCheckpoints <- function(checkpoint_dir) {
  ckpt_files <- list.files(checkpoint_dir, pattern = "\\.rds$",
                            full.names = TRUE)
  if (length(ckpt_files) > 0L) file.remove(ckpt_files)
  # Remove directory if empty
  remaining <- list.files(checkpoint_dir, all.files = TRUE, no.. = TRUE)
  if (length(remaining) == 0L) unlink(checkpoint_dir, recursive = TRUE)
}


#' Summarize Profiles
#'
#' Compute summary statistics across a set of splicing choice profiles.
#'
#' @param profiles A tibble from [buildProfiles()].
#' @return A named list of summary statistics.
#' @examples
#' data(example_profiles)
#' summarizeProfiles(example_profiles)
#' @export
summarizeProfiles <- function(profiles) {
  list(
    n_profiles = nrow(profiles),
    n_genes = dplyr::n_distinct(profiles$gene_id),
    total_events = sum(profiles$n_events),
    mean_events_per_profile = mean(profiles$n_events),
    event_counts = c(
      Alt_TSS = sum(profiles$n_alt_tss),
      Alt_TES = sum(profiles$n_alt_tes),
      A5SS = sum(profiles$n_a5ss),
      A3SS = sum(profiles$n_a3ss),
      SE = sum(profiles$n_se),
      Missing_Internal = sum(profiles$n_missing_internal),
      IR = sum(profiles$n_ir),
      IR_diff = sum(profiles$n_ir_diff),
      Partial_IR = sum(profiles$n_partial_ir)
    ),
    pct_tss_changed = 100 * mean(profiles$tss_changed, na.rm = TRUE),
    pct_tes_changed = 100 * mean(profiles$tes_changed, na.rm = TRUE),
    mean_junction_similarity = mean(profiles$junction_similarity, na.rm = TRUE)
  )
}


# ============================================================================
# Internal Helpers
# ============================================================================

#' Expand nested structures to per-exon rows
#' @keywords internal
.expandStructures <- function(structures) {
  rows <- list()
  for (i in seq_len(nrow(structures))) {
    s <- structures[i, ]
    starts <- s$exon_starts[[1]]
    ends <- s$exon_ends[[1]]
    rows[[i]] <- tibble::tibble(
      isoform_id = s$isoform_id,
      gene_id = s$gene_id,
      chr = s$chr,
      strand = s$strand,
      exon_number = seq_along(starts),
      exon_start = as.integer(starts),
      exon_end = as.integer(ends)
    )
  }
  dplyr::bind_rows(rows)
}

#' Compute gain/loss metrics from events
#' @keywords internal
.computeGainLoss <- function(events) {
  if (nrow(events) == 0L) {
    return(list(n_gain = 0L, n_loss = 0L, bp_gain = 0L, bp_loss = 0L,
                net_bp = 0L))
  }

  gain_events <- events[events$direction == "GAIN", ]
  loss_events <- events[events$direction == "LOSS", ]

  n_gain <- nrow(gain_events)
  n_loss <- nrow(loss_events)

  # Compute bp for each event
  bp_gain <- .sumEventBp(gain_events)
  bp_loss <- .sumEventBp(loss_events)
  net_bp <- bp_loss - bp_gain  # positive = reference is larger

  list(n_gain = n_gain, n_loss = n_loss, bp_gain = bp_gain,
       bp_loss = bp_loss, net_bp = net_bp)
}

#' Sum bp for a set of events
#' @keywords internal
.sumEventBp <- function(events) {
  if (nrow(events) == 0L) return(0L)
  total <- 0L
  for (i in seq_len(nrow(events))) {
    evt <- events[i, ]
    if (!is.na(evt$bp_diff)) {
      total <- total + evt$bp_diff
    } else {
      # IR events: compute from five_prime/three_prime
      total <- total + .computeIrBp(evt)
    }
  }
  total
}

#' Compute bp for IR events where bp_diff is NA
#' @keywords internal
.computeIrBp <- function(event) {
  abs(event$three_prime - event$five_prime) + 1L
}

#' Compute shared junctions between two isoforms (Jaccard)
#' @keywords internal
.computeSharedJunctions <- function(ref_junctions, comp_junctions) {
  n_ref <- length(ref_junctions)
  n_comp <- length(comp_junctions)

  if (n_ref == 0L && n_comp == 0L) {
    return(list(n_shared_junctions = 0L, junction_similarity = NA_real_))
  }

  shared <- intersect(ref_junctions, comp_junctions)
  n_shared <- length(shared)
  total_unique <- length(union(ref_junctions, comp_junctions))
  similarity <- if (total_unique > 0L) n_shared / total_unique else 0.0

  list(n_shared_junctions = n_shared, junction_similarity = similarity)
}

#' Compare union exon coverage between two isoforms
#' @keywords internal
.computeUeComparison <- function(ref_id, comp_id, isoform_union_mapping) {
  ref_ues <- unique(isoform_union_mapping$union_exon_id[
    isoform_union_mapping$isoform_id == ref_id])
  comp_ues <- unique(isoform_union_mapping$union_exon_id[
    isoform_union_mapping$isoform_id == comp_id])

  all_ues <- union(ref_ues, comp_ues)
  shared <- intersect(ref_ues, comp_ues)

  list(
    n_union_exons_total = length(all_ues),
    n_exons_shared = length(shared),
    n_exons_ref_only = length(setdiff(ref_ues, comp_ues)),
    n_exons_comp_only = length(setdiff(comp_ues, ref_ues))
  )
}
