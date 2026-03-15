# Tests for mapSpliceToProtein() and internal helpers

# ============================================================================
# Helper: build minimal test data
# ============================================================================

# + strand gene: 3 exons
# Exon 1: 100-250, Exon 2: 400-550, Exon 3: 700-900
# CDS: 200-800
# Exon 1 CDS: 200-250 = 51 nt
# Exon 2 CDS: 400-550 = 151 nt
# Exon 3 CDS: 700-800 = 101 nt
# Total CDS: 303 nt = 101 AA

.make_plus_structures <- function() {
  tibble::tibble(
    isoform_id = c("tx_plus", "tx_plus_comp"),
    gene_id = c("gene_plus", "gene_plus"),
    chr = c("chr1", "chr1"),
    strand = c("+", "+"),
    n_exons = c(3L, 2L),
    exon_starts = list(c(100L, 400L, 700L), c(100L, 700L)),
    exon_ends = list(c(250L, 550L, 900L), c(250L, 900L)),
    tx_start = c(100L, 100L),
    tx_end = c(900L, 900L),
    n_junctions = c(2L, 1L)
  )
}

.make_plus_cds <- function() {
  tibble::tibble(
    isoform_id = "tx_plus",
    coding_status = "coding",
    cds_start = 200L,
    cds_stop = 800L,
    orf_length = 303L,
    strand = "+"
  )
}

# - strand gene: 3 exons (genomic ascending, but transcript reads right to left)
# Exon 1: 1000-1150, Exon 2: 1300-1450, Exon 3: 1600-1800
# CDS: 1100-1700
# In transcript order (5'->3'): Exon 3 first, then Exon 2, then Exon 1
# Exon 3 CDS (transcript first): 1600-1700 = 101 nt
# Exon 2 CDS: 1300-1450 = 151 nt
# Exon 1 CDS: 1100-1150 = 51 nt
# Total CDS: 303 nt = 101 AA

.make_minus_structures <- function() {
  tibble::tibble(
    isoform_id = c("tx_minus", "tx_minus_comp"),
    gene_id = c("gene_minus", "gene_minus"),
    chr = c("chr1", "chr1"),
    strand = c("-", "-"),
    n_exons = c(3L, 2L),
    exon_starts = list(c(1000L, 1300L, 1600L), c(1000L, 1600L)),
    exon_ends = list(c(1150L, 1450L, 1800L), c(1150L, 1800L)),
    tx_start = c(1000L, 1000L),
    tx_end = c(1800L, 1800L),
    n_junctions = c(2L, 1L)
  )
}

.make_minus_cds <- function() {
  tibble::tibble(
    isoform_id = "tx_minus",
    coding_status = "coding",
    cds_start = 1100L,
    cds_stop = 1700L,
    orf_length = 303L,
    strand = "-"
  )
}


# ============================================================================
# .genomicToCdsPosition() internal helper
# ============================================================================

test_that(".genomicToCdsPosition: + strand at CDS start", {
  # CDS starts at 200, exon 1 starts at 100
  pos <- Isopair:::.genomicToCdsPosition(
    200L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_equal(pos, 1L)
})

test_that(".genomicToCdsPosition: + strand middle of first CDS exon", {
  # Position 210 in exon 1 (CDS starts at 200)
  # CDS-relative: 210 - 200 + 1 = 11
  pos <- Isopair:::.genomicToCdsPosition(
    210L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_equal(pos, 11L)
})

test_that(".genomicToCdsPosition: + strand start of second CDS exon", {
  # Exon 1 CDS: 200-250 = 51 nt
  # Position 400 is start of exon 2, first CDS base in exon 2
  # CDS-relative: 51 + 1 = 52
  pos <- Isopair:::.genomicToCdsPosition(
    400L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_equal(pos, 52L)
})

test_that(".genomicToCdsPosition: + strand at CDS end", {
  # Exon 1 CDS: 200-250 = 51 nt
  # Exon 2 CDS: 400-550 = 151 nt
  # Exon 3 CDS: 700-800 = 101 nt
  # Total = 303
  pos <- Isopair:::.genomicToCdsPosition(
    800L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_equal(pos, 303L)
})

test_that(".genomicToCdsPosition: + strand position in UTR returns NA", {
  # Position 150 is in 5'UTR (before CDS start at 200)
  pos <- Isopair:::.genomicToCdsPosition(
    150L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_true(is.na(pos))
})

test_that(".genomicToCdsPosition: + strand intronic position returns NA", {
  # Position 300 is between exons (intron: 251-399)
  pos <- Isopair:::.genomicToCdsPosition(
    300L, c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L, "+"
  )
  expect_true(is.na(pos))
})

test_that(".genomicToCdsPosition: - strand at translation start (cds_stop)", {
  # On minus strand, translation starts at cds_stop (1700), which is in exon 3
  # CDS-relative position = 1
  pos <- Isopair:::.genomicToCdsPosition(
    1700L, c(1000L, 1300L, 1600L), c(1150L, 1450L, 1800L), 1100L, 1700L, "-"
  )
  expect_equal(pos, 1L)
})

test_that(".genomicToCdsPosition: - strand middle of first transcript exon", {
  # On - strand, exon 3 (1600-1800) is the 5'-most in transcript order
  # CDS portion of exon 3: 1600-1700 = 101 nt
  # Position 1650: from 1700 going toward 1600
  # CDS-relative: 1700 - 1650 + 1 = 51
  pos <- Isopair:::.genomicToCdsPosition(
    1650L, c(1000L, 1300L, 1600L), c(1150L, 1450L, 1800L), 1100L, 1700L, "-"
  )
  expect_equal(pos, 51L)
})

test_that(".genomicToCdsPosition: - strand start of second transcript exon", {
  # Exon 3 CDS: 1600-1700 = 101 nt (transcript positions 1-101)
  # Exon 2 starts at transcript position 102. Its 5' end (in transcript)
  # is genomic 1450 (highest coordinate in exon 2).
  pos <- Isopair:::.genomicToCdsPosition(
    1450L, c(1000L, 1300L, 1600L), c(1150L, 1450L, 1800L), 1100L, 1700L, "-"
  )
  expect_equal(pos, 102L)
})

test_that(".genomicToCdsPosition: - strand at translation stop (cds_start)", {
  # Translation stop at genomic 1100, which is in exon 1
  # Total CDS = 303
  pos <- Isopair:::.genomicToCdsPosition(
    1100L, c(1000L, 1300L, 1600L), c(1150L, 1450L, 1800L), 1100L, 1700L, "-"
  )
  expect_equal(pos, 303L)
})


# ============================================================================
# .computeCdsLength() internal helper
# ============================================================================

test_that(".computeCdsLength: 3-exon + strand", {
  len <- Isopair:::.computeCdsLength(
    c(100L, 400L, 700L), c(250L, 550L, 900L), 200L, 800L
  )
  # Exon 1: 200-250 = 51, Exon 2: 400-550 = 151, Exon 3: 700-800 = 101
  expect_equal(len, 303L)
})

test_that(".computeCdsLength: single exon fully within CDS", {
  len <- Isopair:::.computeCdsLength(500L, 800L, 400L, 900L)
  expect_equal(len, 301L)
})


# ============================================================================
# Test 1: Simple + strand case with one SE event in the CDS
# ============================================================================

test_that("mapSpliceToProtein: + strand SE event in CDS", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # SE event: comparator lost exon 2 (400-550), entirely within CDS
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 400L, three_prime = 550L,
      bp_diff = 151L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)

  expect_equal(nrow(result), 1L)
  expect_equal(result$gene_id, "gene_plus")
  expect_equal(result$event_type, "SE")
  expect_equal(result$genomic_start, 400L)
  expect_equal(result$genomic_end, 550L)

  # CDS-relative: exon 1 CDS = 51 nt, so exon 2 CDS starts at nt 52
  # Exon 2 CDS ends at nt 202 (51 + 151)
  expect_equal(result$cds_nt_start, 52L)
  expect_equal(result$cds_nt_end, 202L)

  # AA positions: ceil(52/3)=18, ceil(202/3)=68
  expect_equal(result$protein_start_aa, 18L)
  expect_equal(result$protein_end_aa, 68L)

  # 151 bp is not divisible by 3 -> affects frame
  expect_true(result$event_affects_frame)
})


# ============================================================================
# Test 2: - strand case
# ============================================================================

test_that("mapSpliceToProtein: - strand SE event in CDS", {
  structures <- .make_minus_structures()
  cds_metadata <- .make_minus_cds()

  # SE event: comparator lost exon 2 (1300-1450), entirely within CDS
  profiles <- tibble::tibble(
    gene_id = "gene_minus",
    reference_isoform_id = "tx_minus",
    comparator_isoform_id = "tx_minus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 1450L, three_prime = 1300L,
      bp_diff = 151L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)

  expect_equal(nrow(result), 1L)
  expect_equal(result$genomic_start, 1300L)
  expect_equal(result$genomic_end, 1450L)

  # On - strand: exon 3 CDS (1600-1700) = 101 nt (transcript positions 1-101)
  # Exon 2 (1300-1450) starts at transcript position 102
  # CDS position of 1450 = 102 (5' end of exon 2 in transcript)
  # CDS position of 1300 = 252 (3' end of exon 2 in transcript)
  expect_equal(result$cds_nt_start, 102L)
  expect_equal(result$cds_nt_end, 252L)

  # AA: ceil(102/3) = 34, ceil(252/3) = 84
  expect_equal(result$protein_start_aa, 34L)
  expect_equal(result$protein_end_aa, 84L)

  expect_true(result$event_affects_frame)
})


# ============================================================================
# Test 3: Event entirely in UTR (should be excluded)
# ============================================================================

test_that("mapSpliceToProtein: event in 5'UTR excluded", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # Event entirely in 5'UTR (100-190, CDS starts at 200)
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "A5SS", direction = "LOSS",
      five_prime = 100L, three_prime = 190L,
      bp_diff = 91L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 0L)
})

test_that("mapSpliceToProtein: event in 3'UTR excluded", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # Event entirely in 3'UTR (850-900, CDS ends at 800)
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "A3SS", direction = "GAIN",
      five_prime = 850L, three_prime = 900L,
      bp_diff = 51L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 0L)
})


# ============================================================================
# Test 4: Multiple events, some in CDS, some not
# ============================================================================

test_that("mapSpliceToProtein: mixed CDS and UTR events", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = c("A5SS", "SE", "A3SS"),
      direction = c("GAIN", "LOSS", "LOSS"),
      five_prime = c(100L, 400L, 850L),
      three_prime = c(190L, 550L, 900L),
      bp_diff = c(91L, 151L, 51L),
      stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)

  # Only the SE event (400-550) overlaps CDS; the other two are in UTRs
  expect_equal(nrow(result), 1L)
  expect_equal(result$event_type, "SE")
})


# ============================================================================
# Test 5: Buffer parameter verification
# ============================================================================

test_that("mapSpliceToProtein: buffer_aa extends boundaries", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 400L, three_prime = 550L,
      bp_diff = 151L, stringsAsFactors = FALSE
    ))
  )

  # Default buffer = 5 AA
  result_default <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(result_default$protein_start_aa_buffered,
               result_default$protein_start_aa - 5L)
  expect_equal(result_default$protein_end_aa_buffered,
               result_default$protein_end_aa + 5L)

  # Buffer = 0
  result_no_buf <- mapSpliceToProtein(profiles, cds_metadata, structures,
                                       buffer_aa = 0L)
  expect_equal(result_no_buf$protein_start_aa_buffered,
               result_no_buf$protein_start_aa)
  expect_equal(result_no_buf$protein_end_aa_buffered,
               result_no_buf$protein_end_aa)

  # Large buffer clamps to valid range
  result_big <- mapSpliceToProtein(profiles, cds_metadata, structures,
                                    buffer_aa = 500L)
  expect_equal(result_big$protein_start_aa_buffered, 1L)
  # max_aa = ceil(303/3) = 101
  expect_equal(result_big$protein_end_aa_buffered, 101L)
})


# ============================================================================
# Test 6: Event partially overlapping CDS (clipped to CDS boundaries)
# ============================================================================

test_that("mapSpliceToProtein: event spanning UTR/CDS boundary is clipped", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # Event spans from 5'UTR into CDS: 150-230 (CDS starts at 200)
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "A5SS", direction = "LOSS",
      five_prime = 150L, three_prime = 230L,
      bp_diff = 81L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 1L)

  # Clipped to CDS: 200-230
  # CDS-relative: 1 to 31
  expect_equal(result$cds_nt_start, 1L)
  expect_equal(result$cds_nt_end, 31L)

  # AA: ceil(1/3) = 1, ceil(31/3) = 11
  expect_equal(result$protein_start_aa, 1L)
  expect_equal(result$protein_end_aa, 11L)
})


# ============================================================================
# Test 7: Non-coding reference isoform (skip pair)
# ============================================================================

test_that("mapSpliceToProtein: non-coding reference skipped", {
  structures <- .make_plus_structures()
  cds_metadata <- tibble::tibble(
    isoform_id = "tx_plus",
    coding_status = "unknown",
    cds_start = NA_integer_,
    cds_stop = NA_integer_,
    orf_length = NA_integer_,
    strand = NA_character_
  )

  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 400L, three_prime = 550L,
      bp_diff = 151L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 0L)
})


# ============================================================================
# Test 8: Single-exon gene
# ============================================================================

test_that("mapSpliceToProtein: single-exon gene works", {
  structures <- tibble::tibble(
    isoform_id = c("tx_single", "tx_single_comp"),
    gene_id = c("gene_single", "gene_single"),
    chr = c("chr1", "chr1"),
    strand = c("+", "+"),
    n_exons = c(1L, 1L),
    exon_starts = list(100L, 100L),
    exon_ends = list(600L, 550L),
    tx_start = c(100L, 100L),
    tx_end = c(600L, 550L),
    n_junctions = c(0L, 0L)
  )

  cds_metadata <- tibble::tibble(
    isoform_id = "tx_single",
    coding_status = "coding",
    cds_start = 200L,
    cds_stop = 500L,
    orf_length = 301L,
    strand = "+"
  )

  # A3SS event at the 3' end of the exon, partially in CDS
  profiles <- tibble::tibble(
    gene_id = "gene_single",
    reference_isoform_id = "tx_single",
    comparator_isoform_id = "tx_single_comp",
    detailed_events = list(data.frame(
      event_type = "A3SS", direction = "LOSS",
      five_prime = 480L, three_prime = 600L,
      bp_diff = 121L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 1L)

  # Clipped to CDS: 480-500
  # CDS-relative: 480 - 200 + 1 = 281, 500 - 200 + 1 = 301
  expect_equal(result$cds_nt_start, 281L)
  expect_equal(result$cds_nt_end, 301L)
})


# ============================================================================
# Test 9: event_affects_frame logic
# ============================================================================

test_that("mapSpliceToProtein: in-frame event has event_affects_frame=FALSE", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # bp_diff = 150, which IS divisible by 3
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 400L, three_prime = 550L,
      bp_diff = 150L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_false(result$event_affects_frame)
})


# ============================================================================
# Test 10: Input validation
# ============================================================================

test_that("mapSpliceToProtein: validates required columns", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()
  bad_profiles <- tibble::tibble(gene_id = "x")

  expect_error(
    mapSpliceToProtein(bad_profiles, cds_metadata, structures),
    "missing required columns"
  )
})

test_that("mapSpliceToProtein: validates buffer_aa", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "SE", direction = "LOSS",
      five_prime = 400L, three_prime = 550L,
      bp_diff = 151L, stringsAsFactors = FALSE
    ))
  )

  expect_error(
    mapSpliceToProtein(profiles, cds_metadata, structures, buffer_aa = -1L),
    "non-negative"
  )
})


# ============================================================================
# Test 11: Empty inputs
# ============================================================================

test_that("mapSpliceToProtein: empty profiles returns empty tibble", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  profiles <- tibble::tibble(
    gene_id = character(0),
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0),
    detailed_events = list()
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 0L)
  expect_true(all(c("gene_id", "protein_start_aa", "event_affects_frame") %in%
                    names(result)))
})

test_that("mapSpliceToProtein: no events returns empty tibble", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = character(0), direction = character(0),
      five_prime = integer(0), three_prime = integer(0),
      bp_diff = integer(0), stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 0L)
})


# ============================================================================
# Test 12: Event at exact CDS boundary
# ============================================================================

test_that("mapSpliceToProtein: event at exact CDS start boundary", {
  structures <- .make_plus_structures()
  cds_metadata <- .make_plus_cds()

  # Event exactly at CDS start: 200-210
  profiles <- tibble::tibble(
    gene_id = "gene_plus",
    reference_isoform_id = "tx_plus",
    comparator_isoform_id = "tx_plus_comp",
    detailed_events = list(data.frame(
      event_type = "A5SS", direction = "LOSS",
      five_prime = 200L, three_prime = 210L,
      bp_diff = 11L, stringsAsFactors = FALSE
    ))
  )

  result <- mapSpliceToProtein(profiles, cds_metadata, structures)
  expect_equal(nrow(result), 1L)
  expect_equal(result$cds_nt_start, 1L)
  expect_equal(result$cds_nt_end, 11L)
})
