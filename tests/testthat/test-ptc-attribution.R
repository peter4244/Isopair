# Tests for attributePtcEvents(), attribute3UtrSplice(), traceReferenceAtg()

# ==========================================================================
# Helper: build minimal test data
# ==========================================================================

# Plus-strand gene with 4 exons:
# Exon 1: [100,200], Exon 2: [300,400], Exon 3: [500,600], Exon 4: [700,800]
make_test_structures <- function() {
  tibble::tibble(
    isoform_id = c("ref1", "comp1", "comp2", "comp3", "comp4"),
    gene_id = "gene1",
    strand = "+",
    n_exons = c(4L, 4L, 3L, 4L, 4L),
    tx_start = c(100L, 100L, 100L, 100L, 100L),
    tx_end = c(800L, 800L, 600L, 800L, 800L),
    exon_starts = list(
      c(100L, 300L, 500L, 700L),  # ref1: 4 exons
      c(100L, 300L, 500L, 700L),  # comp1: same structure
      c(100L, 300L, 500L),         # comp2: missing exon 4
      c(100L, 300L, 500L, 700L),  # comp3: same structure
      c(100L, 300L, 500L, 700L)   # comp4: same structure
    ),
    exon_ends = list(
      c(200L, 400L, 600L, 800L),
      c(200L, 400L, 600L, 800L),
      c(200L, 400L, 600L),
      c(200L, 400L, 600L, 800L),
      c(200L, 400L, 600L, 800L)
    )
  )
}

make_test_cds <- function() {
  tibble::tibble(
    isoform_id = c("ref1", "comp1", "comp2", "comp3", "comp4"),
    coding_status = "coding",
    strand = "+",
    cds_start = c(150L, 150L, 150L, 350L, 150L),  # comp3: different CDS start
    cds_stop = c(750L, 750L, 550L, 750L, 750L),    # comp2: earlier stop
    orf_length = c(200L, 200L, 133L, 133L, 200L)
  )
}

# ==========================================================================
# genomicToTranscript()
# ==========================================================================

test_that("genomicToTranscript: plus strand basic", {
  # Exon 1: [100,200] = 101 bp → tx pos 1-101
  # Exon 2: [300,400] = 101 bp → tx pos 102-202
  expect_equal(genomicToTranscript(100, c(100, 300), c(200, 400), "+"), 1L)
  expect_equal(genomicToTranscript(200, c(100, 300), c(200, 400), "+"), 101L)
  expect_equal(genomicToTranscript(300, c(100, 300), c(200, 400), "+"), 102L)
  expect_equal(genomicToTranscript(400, c(100, 300), c(200, 400), "+"), 202L)
})

test_that("genomicToTranscript: minus strand basic", {
  # Minus strand: exons reversed, tx pos 1 = highest genomic
  # Exon [300,400] is first in transcript (highest coords)
  # Exon [100,200] is second
  expect_equal(genomicToTranscript(400, c(100, 300), c(200, 400), "-"), 1L)
  expect_equal(genomicToTranscript(300, c(100, 300), c(200, 400), "-"), 101L)
  expect_equal(genomicToTranscript(200, c(100, 300), c(200, 400), "-"), 102L)
  expect_equal(genomicToTranscript(100, c(100, 300), c(200, 400), "-"), 202L)
})

test_that("genomicToTranscript: intronic position returns NA", {
  expect_true(is.na(genomicToTranscript(250, c(100, 300), c(200, 400), "+")))
})

test_that("genomicToTranscript: NA input returns NA", {
  expect_true(is.na(genomicToTranscript(NA, c(100, 300), c(200, 400), "+")))
})

# ==========================================================================
# transcriptToGenomic()
# ==========================================================================

test_that("transcriptToGenomic: roundtrip with genomicToTranscript", {
  exon_s <- c(100L, 300L, 500L)
  exon_e <- c(200L, 400L, 600L)

  for (strand in c("+", "-")) {
    for (gpos in c(100, 150, 200, 300, 350, 400, 500, 600)) {
      tx <- genomicToTranscript(gpos, exon_s, exon_e, strand)
      back <- transcriptToGenomic(tx, exon_s, exon_e, strand)
      expect_equal(back, gpos,
        info = sprintf("roundtrip failed: gpos=%d, strand=%s", gpos, strand))
    }
  }
})

test_that("transcriptToGenomic: beyond transcript length returns NA", {
  expect_true(is.na(transcriptToGenomic(999, c(100), c(200), "+")))
})

# ==========================================================================
# getStrandAwareStop()
# ==========================================================================

test_that("getStrandAwareStop: plus strand uses cds_stop", {
  cds <- data.frame(
    isoform_id = "tx1", cds_start = 100, cds_stop = 400, strand = "+",
    coding_status = "coding")
  result <- getStrandAwareStop("tx1", cds)
  expect_equal(unname(result), 400)
})

test_that("getStrandAwareStop: minus strand uses cds_start", {
  cds <- data.frame(
    isoform_id = "tx1", cds_start = 100, cds_stop = 400, strand = "-",
    coding_status = "coding")
  result <- getStrandAwareStop("tx1", cds)
  expect_equal(unname(result), 100)
})

test_that("getStrandAwareStop: missing ID returns NA", {
  cds <- data.frame(
    isoform_id = "tx1", cds_start = 100, cds_stop = 400, strand = "+",
    coding_status = "coding")
  result <- getStrandAwareStop("tx_missing", cds)
  expect_true(is.na(result))
})

# ==========================================================================
# computeUtrLengths()
# ==========================================================================

test_that("computeUtrLengths: plus strand basic", {
  structures <- tibble::tibble(
    isoform_id = "tx1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 100L, tx_end = 500L,
    exon_starts = list(100L), exon_ends = list(500L))
  cds <- tibble::tibble(
    isoform_id = "tx1", coding_status = "coding", strand = "+",
    cds_start = 200L, cds_stop = 400L, orf_length = 67L)

  result <- computeUtrLengths(structures, cds)
  expect_equal(result$utr5_bp, 100L)  # 100-199
  expect_equal(result$utr3_bp, 100L)  # 401-500
  expect_equal(result$cds_bp, 201L)   # 200-400
  expect_equal(result$tx_length, 401L)
})

test_that("computeUtrLengths: minus strand", {
  structures <- tibble::tibble(
    isoform_id = "tx1", gene_id = "g1", strand = "-", n_exons = 1L,
    tx_start = 100L, tx_end = 500L,
    exon_starts = list(100L), exon_ends = list(500L))
  cds <- tibble::tibble(
    isoform_id = "tx1", coding_status = "coding", strand = "-",
    cds_start = 200L, cds_stop = 400L, orf_length = 67L)

  # Minus strand: 5'UTR is at high coords (401-500), 3'UTR at low (100-199)
  result <- computeUtrLengths(structures, cds)
  expect_equal(result$utr5_bp, 100L)  # 401-500 (upstream on - strand)
  expect_equal(result$utr3_bp, 100L)  # 100-199 (downstream on - strand)
})

test_that("computeUtrLengths: multi-exon with spliced UTR", {
  structures <- tibble::tibble(
    isoform_id = "tx1", gene_id = "g1", strand = "+", n_exons = 3L,
    tx_start = 100L, tx_end = 600L,
    exon_starts = list(c(100L, 300L, 500L)),
    exon_ends = list(c(200L, 400L, 600L)))
  cds <- tibble::tibble(
    isoform_id = "tx1", coding_status = "coding", strand = "+",
    cds_start = 150L, cds_stop = 550L, orf_length = 133L)

  result <- computeUtrLengths(structures, cds)
  # 5'UTR: exon1 100-149 = 50 bp
  expect_equal(result$utr5_bp, 50L)
  # 3'UTR: exon3 551-600 = 50 bp
  expect_equal(result$utr3_bp, 50L)
})

test_that("computeUtrLengths: non-coding excluded", {
  structures <- tibble::tibble(
    isoform_id = "tx1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 100L, tx_end = 500L,
    exon_starts = list(100L), exon_ends = list(500L))
  cds <- tibble::tibble(
    isoform_id = "tx1", coding_status = "unknown", strand = "+",
    cds_start = NA_integer_, cds_stop = NA_integer_, orf_length = NA_integer_)

  result <- computeUtrLengths(structures, cds)
  expect_equal(nrow(result), 0L)
})

# ==========================================================================
# extractAllEvents()
# ==========================================================================

test_that("extractAllEvents: basic concatenation", {
  profiles <- data.frame(
    comparator_isoform_id = c("c1", "c2"),
    reference_isoform_id = c("r1", "r2"))
  profiles$detailed_events <- list(
    data.frame(event_type = c("SE", "A3SS"), bp_diff = c(100L, 50L)),
    data.frame(event_type = "SE", bp_diff = 200L))

  result <- extractAllEvents(profiles)
  expect_equal(nrow(result), 3L)
  expect_equal(sum(result$event_type == "SE"), 2L)
})

test_that("extractAllEvents: column subsetting", {
  profiles <- data.frame(
    comparator_isoform_id = "c1",
    reference_isoform_id = "r1")
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", bp_diff = 100L, direction = "GAIN"))

  result <- extractAllEvents(profiles, cols = c("event_type", "direction"))
  expect_equal(ncol(result), 2L)
  expect_false("bp_diff" %in% names(result))
})

test_that("extractAllEvents: empty events handled", {
  profiles <- data.frame(
    comparator_isoform_id = c("c1", "c2"),
    reference_isoform_id = c("r1", "r2"))
  profiles$detailed_events <- list(NULL, data.frame(event_type = "SE", bp_diff = 100L))

  result <- extractAllEvents(profiles)
  expect_equal(nrow(result), 1L)
})

# ==========================================================================
# attributePtcEvents()
# ==========================================================================

test_that("attributePtcEvents: frameshift mechanism detected", {
  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  fw_events <- tibble::tibble(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1",
    event_type = "SE",
    genomic_start = 350L, genomic_end = 400L,
    is_frameshift = TRUE)

  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", five_prime = 350, three_prime = 400, bp_diff = 51))

  ptc_pos <- c(comp1 = 600L)
  strand_v <- c(comp1 = "+")
  is_fs <- c(comp1 = TRUE)

  result <- attributePtcEvents(pairs, fw_events, profiles, ptc_pos,
                                strand_vec = strand_v,
                                is_frameshift_vec = is_fs)

  expect_equal(nrow(result), 1L)
  expect_equal(result$mechanism, "Frameshift")
  expect_equal(result$ptc_causing_event, "SE")
  expect_equal(result$attribution, "direct")
})

test_that("attributePtcEvents: in-frame stop via coordinate containment", {
  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  fw_events <- tibble::tibble(
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    event_type = character(0),
    genomic_start = integer(0), genomic_end = integer(0),
    is_frameshift = logical(0))

  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", five_prime = 500, three_prime = 600, bp_diff = 99))

  ptc_pos <- c(comp1 = 550L)  # inside the SE event
  strand_v <- c(comp1 = "+")
  is_fs <- c(comp1 = FALSE)

  result <- attributePtcEvents(pairs, fw_events, profiles, ptc_pos,
                                strand_vec = strand_v,
                                is_frameshift_vec = is_fs)

  expect_equal(result$mechanism, "In-frame stop")
  expect_equal(result$ptc_causing_event, "SE")
})

test_that("attributePtcEvents: split-codon buffer catches near-boundary events", {
  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  fw_events <- tibble::tibble(
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    event_type = character(0),
    genomic_start = integer(0), genomic_end = integer(0),
    is_frameshift = logical(0))

  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = "Partial_IR_5", five_prime = 500, three_prime = 598, bp_diff = 99))

  # Stop at 600 — 2 bp outside the event boundary (598). Buffer should catch it.
  ptc_pos <- c(comp1 = 600L)
  strand_v <- c(comp1 = "+")
  is_fs <- c(comp1 = FALSE)

  result <- attributePtcEvents(pairs, fw_events, profiles, ptc_pos,
                                strand_vec = strand_v,
                                is_frameshift_vec = is_fs)

  expect_equal(result$mechanism, "In-frame stop")
  expect_equal(result$attribution, "direct")
})

test_that("attributePtcEvents: region filtering with atg_genomic_pos", {
  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  # Frameshift event at [800,850] — OUTSIDE the ATG(100)-to-stop(600) region
  fw_events <- tibble::tibble(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1",
    event_type = "SE",
    genomic_start = 800L, genomic_end = 850L,
    is_frameshift = TRUE)

  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", five_prime = 800, three_prime = 850, bp_diff = 51))

  ptc_pos <- c(comp1 = 600L)
  atg_pos <- c(comp1 = 100L)
  strand_v <- c(comp1 = "+")

  result <- attributePtcEvents(pairs, fw_events, profiles, ptc_pos,
                                atg_genomic_pos = atg_pos,
                                strand_vec = strand_v)

  # Should be unresolved — the frameshift is outside the ATG-to-stop region
  expect_equal(result$mechanism, "Frameshift")
  expect_equal(result$attribution, "unresolved")
})

test_that("attributePtcEvents: composite key handles duplicate comparators", {
  # Same comparator in two pairs with different references
  pairs <- data.frame(
    comparator_isoform_id = c("comp1", "comp1"),
    reference_isoform_id = c("ref1", "ref2"))

  fw_events <- tibble::tibble(
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    event_type = character(0),
    genomic_start = integer(0), genomic_end = integer(0),
    is_frameshift = logical(0))

  profiles <- data.frame(
    reference_isoform_id = c("ref1", "ref2"),
    comparator_isoform_id = c("comp1", "comp1"))
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", five_prime = 500, three_prime = 600, bp_diff = 99),
    data.frame(event_type = "A3SS", five_prime = 500, three_prime = 600, bp_diff = 30))

  ptc_pos <- c(comp1 = 550L)
  strand_v <- c(comp1 = "+")
  is_fs <- c(comp1 = FALSE)

  result <- attributePtcEvents(pairs, fw_events, profiles, ptc_pos,
                                strand_vec = strand_v,
                                is_frameshift_vec = is_fs)

  expect_equal(nrow(result), 2L)
  # First pair should get SE (from ref1::comp1 events)
  expect_equal(result$ptc_causing_event[1], "SE")
  # Second pair should get A3SS (from ref2::comp1 events)
  expect_equal(result$ptc_causing_event[2], "A3SS")
})

# ==========================================================================
# attribute3UtrSplice()
# ==========================================================================

test_that("attribute3UtrSplice: nearest downstream event found", {
  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = c("SE", "Alt_TES"),
               five_prime = c(100, 700), three_prime = c(200, 800)))

  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  result <- attribute3UtrSplice(
    pairs, profiles,
    stop_genomic_pos = c(comp1 = 600L),
    strand_vec = c(comp1 = "+"))

  expect_equal(result$ptc_causing_event, "Alt_TES")
  expect_equal(result$attribution, "direct")
})

test_that("attribute3UtrSplice: no downstream event → unresolved", {
  profiles <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")
  profiles$detailed_events <- list(
    data.frame(event_type = "SE", five_prime = 100, three_prime = 200))

  pairs <- data.frame(
    comparator_isoform_id = "comp1",
    reference_isoform_id = "ref1")

  result <- attribute3UtrSplice(
    pairs, profiles,
    stop_genomic_pos = c(comp1 = 600L),
    strand_vec = c(comp1 = "+"))

  expect_equal(result$attribution, "unresolved")
})

# ==========================================================================
# traceReferenceAtg()
# ==========================================================================

test_that("traceReferenceAtg: basic effectively_ptc case", {
  structures <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    gene_id = "g1", strand = "+", n_exons = 3L,
    tx_start = 1L, tx_end = 300L,
    exon_starts = list(c(1L, 101L, 201L), c(1L, 101L, 201L)),
    exon_ends = list(c(100L, 200L, 300L), c(100L, 200L, 300L)))

  cds <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    coding_status = "coding", strand = "+",
    cds_start = c(10L, 50L),   # comp has different CDS start
    cds_stop = c(280L, 280L),
    orf_length = c(90L, 77L))

  # Build sequences: ref has ATG at position 10 (tx pos 10)
  # Comp has ATG at position 10 too (same genomic pos)
  # Need sequences where ref ATG traces to a premature stop in comp
  ref_seq <- paste0(
    paste(rep("N", 9), collapse = ""),   # 9 nt before ATG
    "ATG",                                # ATG at pos 10
    paste(rep("NNN", 80), collapse = ""), # 240 nt of codons (no stops)
    "TAA",                                # stop at pos 253
    paste(rep("N", 47), collapse = ""))   # padding to 300

  # Comp: same ATG but a stop codon earlier
  comp_seq <- paste0(
    paste(rep("N", 9), collapse = ""),
    "ATG",
    paste(rep("NNN", 20), collapse = ""),  # 60 nt of codons
    "TAG",                                  # premature stop at pos 72
    paste(rep("N", 228), collapse = ""))   # padding to 300

  sequences <- c(ref1 = ref_seq, comp1 = comp_seq)

  pairs <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")

  result <- traceReferenceAtg(pairs, structures, cds, sequences)

  expect_equal(nrow(result), 1L)
  expect_true(result$ref_atg_exonic_in_comp)
  expect_equal(result$category, "effectively_ptc")
  expect_true(result$comp_orf_length < result$ref_orf_length)
  expect_true(result$n_downstream_ejc > 0L)
})

test_that("traceReferenceAtg: ref_atg_lost when ATG not exonic", {
  structures <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    gene_id = "g1", strand = "+", n_exons = c(2L, 1L),
    tx_start = c(1L, 201L), tx_end = c(400L, 400L),
    exon_starts = list(c(1L, 201L), 201L),  # comp missing first exon
    exon_ends = list(c(100L, 400L), 400L))

  cds <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    coding_status = "coding", strand = "+",
    cds_start = c(50L, 250L), cds_stop = c(350L, 350L),
    orf_length = c(100L, 33L))

  pairs <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")

  # Ref ATG at genomic 50, which is in exon [1,100] of ref but NOT in comp
  result <- traceReferenceAtg(pairs, structures, cds,
                               sequences = c(ref1 = "A", comp1 = "A"))

  expect_equal(result$category, "ref_atg_lost")
  expect_false(result$ref_atg_exonic_in_comp)
})

test_that("traceReferenceAtg: resolve_alt_start recovers alt-start PTC categories", {
  # Comparator is missing ref's first exon (ref ATG position not exonic).
  # Default resolve_alt_start = FALSE -> ref_atg_lost (unchanged behavior).
  # resolve_alt_start = TRUE -> scan comp for first viable ATG and classify.
  structures <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    gene_id = "g1", strand = "+", n_exons = c(2L, 2L),
    tx_start = c(1L, 201L), tx_end = c(600L, 600L),
    exon_starts = list(c(1L, 201L), c(201L, 401L)),
    exon_ends   = list(c(100L, 600L), c(300L, 600L))
  )
  cds <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    coding_status = c("coding", "non-coding"), strand = "+",
    cds_start = c(50L, NA_integer_), cds_stop = c(350L, NA_integer_),
    orf_length = c(100L, NA_integer_))

  # Build a comparator sequence that has a downstream junction after
  # the in-frame stop so the category should be alt_start_effectively_ptc.
  # Comparator transcript length = (300-201+1) + (600-401+1) = 100 + 200 = 300.
  # Put an ATG at tx pos 20 and a stop at tx pos 29 (inclusive), so the ORF
  # is 9 nt (below default min 30 -> should be rejected).
  # Put a second ATG at tx pos 50 followed by plenty of non-stop codons
  # and an in-frame TAA at tx pos 98. First-exon length is 100, so the
  # stop_end = 100 is within the first exon and there's a junction at 101.
  # Threshold = 50 means junctions > stop_end + 49 count as downstream.
  # Junction 101 > (100 + 49) = 149? No. So need a later stop.
  # Simpler: construct a deterministic sequence with controlled positions.
  len <- 300L
  s <- rep("A", len)
  # short alt ATG at pos 10 with immediate stop at 13 -> 3-nt ORF, rejected
  s[10:12] <- c("A","T","G"); s[13:15] <- c("T","A","A")
  # viable ATG at pos 50, stop at pos 80 (30 nt ORF, just meets default min)
  s[50:52] <- c("A","T","G"); s[80:82] <- c("T","A","A")
  # fill interior codons with serines (TCN) to avoid premature stops
  fill_pos <- seq(53L, 77L, by = 3L)
  for (p in fill_pos) s[p:(p + 2L)] <- c("T","C","A")
  comp_seq <- paste0(s, collapse = "")

  sequences <- c(ref1 = paste0(rep("A", 400L), collapse = ""),
                 comp1 = comp_seq)

  pairs <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")

  # Default: ref_atg_lost, as before.
  r_default <- traceReferenceAtg(pairs, structures, cds, sequences)
  expect_equal(r_default$category, "ref_atg_lost")
  expect_true(is.na(r_default$alt_start_tx_pos))

  # With resolution: either alt_start_effectively_ptc or
  # alt_start_no_downstream_ejc depending on EJC geometry. Here the first
  # exon ends at tx pos 100 (so junction position in tx space is 100), the
  # stop ends at tx pos 82, and ejc_threshold = 50 -> junction must be
  # > 82 + 49 = 131 to count as downstream. Junction 100 does NOT count
  # -> expect alt_start_no_downstream_ejc.
  r_resolved <- traceReferenceAtg(pairs, structures, cds, sequences,
                                  resolve_alt_start = TRUE)
  expect_equal(r_resolved$category, "alt_start_no_downstream_ejc")
  expect_equal(r_resolved$alt_start_tx_pos, 50L)
  expect_equal(r_resolved$alt_start_orf_length, 30L)

  # Raising ejc_threshold doesn't matter here since there's no
  # downstream junction close enough. But if we move the alt-start stop
  # earlier so the junction IS downstream, we should get the PTC category.
  # Construct a version where stop is at pos 40 (first-exon junction at
  # 100 becomes > 40 + 49 = 89 -> downstream).
  s2 <- s
  s2[40:42] <- c("T","A","A")
  # make sure no earlier in-frame stop between ATG(50) and old stop(80)
  sequences2 <- sequences
  sequences2["comp1"] <- paste0(s2, collapse = "")
  # Move viable-ORF ATG earlier so stop at 40 is its in-frame stop.
  s3 <- rep("A", len)
  s3[10:12] <- c("A","T","G"); s3[40:42] <- c("T","A","A")
  fill_pos2 <- seq(13L, 37L, by = 3L)
  for (p in fill_pos2) s3[p:(p + 2L)] <- c("T","C","A")
  sequences3 <- c(ref1 = sequences["ref1"], comp1 = paste0(s3, collapse = ""))
  names(sequences3)[1] <- "ref1"
  # 30 nt ORF (10..39). Stop_end = 42. Junction at 100 > 42 + 49 = 91 -> yes.
  r_ptc <- traceReferenceAtg(pairs, structures, cds, sequences3,
                             resolve_alt_start = TRUE)
  expect_equal(r_ptc$category, "alt_start_effectively_ptc")
  expect_equal(r_ptc$alt_start_tx_pos, 10L)
  expect_equal(r_ptc$alt_start_orf_length, 30L)
  expect_gt(r_ptc$n_downstream_ejc, 0L)
})

test_that("traceReferenceAtg: resolve_alt_start returns ref_atg_lost_no_viable_start when no ATG", {
  structures <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    gene_id = "g1", strand = "+", n_exons = c(2L, 1L),
    tx_start = c(1L, 201L), tx_end = c(400L, 400L),
    exon_starts = list(c(1L, 201L), 201L),
    exon_ends   = list(c(100L, 400L), 400L))
  cds <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    coding_status = c("coding", "non-coding"), strand = "+",
    cds_start = c(50L, NA_integer_), cds_stop = c(350L, NA_integer_),
    orf_length = c(100L, NA_integer_))
  pairs <- data.frame(reference_isoform_id = "ref1",
                      comparator_isoform_id = "comp1")
  sequences <- c(ref1 = paste0(rep("A", 400L), collapse = ""),
                 comp1 = paste0(rep("A", 200L), collapse = ""))
  r <- traceReferenceAtg(pairs, structures, cds, sequences,
                         resolve_alt_start = TRUE)
  expect_equal(r$category, "ref_atg_lost_no_viable_start")
  expect_true(is.na(r$alt_start_tx_pos))
})

test_that("traceReferenceAtg: no_ref_cds for non-coding reference", {
  structures <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = 300L,
    exon_starts = list(1L, 1L), exon_ends = list(300L, 300L))

  cds <- tibble::tibble(
    isoform_id = c("ref1", "comp1"),
    coding_status = c("unknown", "coding"), strand = "+",
    cds_start = c(NA_integer_, 50L), cds_stop = c(NA_integer_, 250L),
    orf_length = c(NA_integer_, 67L))

  pairs <- data.frame(
    reference_isoform_id = "ref1",
    comparator_isoform_id = "comp1")

  result <- traceReferenceAtg(pairs, structures, cds,
                               sequences = c(ref1 = "A", comp1 = "A"))

  expect_equal(result$category, "no_ref_cds")
})
