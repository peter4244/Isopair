# ==============================================================================
# test-uorf.R — Tests for 5'UTR feature scanning, uORF detection,
#                comparison, and attribution
# ==============================================================================

# --- Helpers for constructing test data ---------------------------------------

# Build a minimal structures tibble
make_structures <- function(isoform_id, gene_id, chr = "chr1", strand = "+",
                            exon_starts, exon_ends) {
  tibble::tibble(
    isoform_id = isoform_id,
    gene_id = gene_id,
    chr = chr,
    strand = strand,
    n_exons = length(exon_starts),
    exon_starts = list(as.integer(exon_starts)),
    exon_ends = list(as.integer(exon_ends))
  )
}

# Build a minimal cds_metadata tibble
make_cds <- function(isoform_id, coding_status = "coding",
                     cds_start, cds_stop, orf_length = NA_integer_,
                     strand = "+") {
  tibble::tibble(
    isoform_id = isoform_id,
    coding_status = coding_status,
    cds_start = as.integer(cds_start),
    cds_stop = as.integer(cds_stop),
    orf_length = as.integer(orf_length),
    strand = strand
  )
}


# ==============================================================================
# .extract5UtrCoords()
# ==============================================================================

test_that(".extract5UtrCoords works on plus strand single exon", {
  # Exon 100-300, CDS starts at 150
  # 5'UTR = 100-149
  coords <- Isopair:::.extract5UtrCoords("+", 100L, 300L, 150L, 280L)
  expect_equal(nrow(coords), 1)
  expect_equal(coords$start, 100L)
  expect_equal(coords$end, 149L)
})

test_that(".extract5UtrCoords works on plus strand multi exon", {
  # Exon1: 100-200, Exon2: 400-600, CDS starts at 450
  # 5'UTR: all of exon1 (100-200) + part of exon2 (400-449)
  coords <- Isopair:::.extract5UtrCoords("+", c(100L, 400L), c(200L, 600L),
                                          400L, 580L)
  expect_equal(nrow(coords), 1)
  # Exon1 is entirely before cds_start=400, so 100-200
  # Exon2 starts at 400 = cds_start, no partial 5'UTR in exon2
  expect_equal(coords$start, 100L)
  expect_equal(coords$end, 200L)
})

test_that(".extract5UtrCoords handles CDS starting mid-exon", {
  # Exon1: 100-200, Exon2: 400-600, CDS starts at 450
  coords <- Isopair:::.extract5UtrCoords("+", c(100L, 400L), c(200L, 600L),
                                          450L, 580L)
  expect_equal(nrow(coords), 2)
  expect_equal(coords$start, c(100L, 400L))
  expect_equal(coords$end, c(200L, 449L))
})

test_that(".extract5UtrCoords returns empty for no 5'UTR (+ strand)", {
  # CDS starts at exon start
  coords <- Isopair:::.extract5UtrCoords("+", 100L, 300L, 100L, 280L)
  expect_equal(nrow(coords), 0)
})

test_that(".extract5UtrCoords works on minus strand single exon", {
  # Exon 100-300, CDS: cds_start=120, cds_stop=250
  # Minus strand: translation starts at cds_stop=250
  # 5'UTR = exonic bases > 250 → 251-300
  coords <- Isopair:::.extract5UtrCoords("-", 100L, 300L, 120L, 250L)
  expect_equal(nrow(coords), 1)
  expect_equal(coords$start, 251L)
  expect_equal(coords$end, 300L)
})

test_that(".extract5UtrCoords works on minus strand multi exon", {
  # Exon1: 100-200, Exon2: 400-600
  # CDS: cds_start=100, cds_stop=450
  # Minus strand: translation starts at cds_stop=450
  # 5'UTR = exonic bases > 450 → 451-600 (part of exon2)
  coords <- Isopair:::.extract5UtrCoords("-", c(100L, 400L), c(200L, 600L),
                                          100L, 450L)
  expect_equal(nrow(coords), 1)
  expect_equal(coords$start, 451L)
  expect_equal(coords$end, 600L)
})

test_that(".extract5UtrCoords minus strand with full exon in 5'UTR", {
  # Exon1: 100-200, Exon2: 400-500, Exon3: 700-900
  # CDS: cds_start=100, cds_stop=450
  # Minus strand: translation starts at cds_stop=450
  # 5'UTR = exon3 entirely (700-900) + part of exon2 (451-500)
  coords <- Isopair:::.extract5UtrCoords("-", c(100L, 400L, 700L),
                                          c(200L, 500L, 900L),
                                          100L, 450L)
  expect_equal(nrow(coords), 2)
  # Sorted by genomic position
  expect_equal(coords$start, c(451L, 700L))
  expect_equal(coords$end, c(500L, 900L))
})

test_that(".extract5UtrCoords returns empty for no 5'UTR (- strand)", {
  # Minus strand, CDS goes to end of last exon
  coords <- Isopair:::.extract5UtrCoords("-", 100L, 300L, 100L, 300L)
  expect_equal(nrow(coords), 0)
})


# ==============================================================================
# .compute5UtrLength()
# ==============================================================================

test_that(".compute5UtrLength gives correct length plus strand", {
  # Exon 100-300, CDS starts at 150: 5'UTR = 150-100 = 50 bp
  len <- Isopair:::.compute5UtrLength("+", 100L, 300L, 150L, 280L)
  expect_equal(len, 50L)
})

test_that(".compute5UtrLength gives correct length multi-exon plus strand", {
  # Exon1: 100-200 (101bp), Exon2: 400-600, CDS starts at 450
  # 5'UTR = 101 + (449-400+1) = 101 + 50 = 151
  len <- Isopair:::.compute5UtrLength("+", c(100L, 400L), c(200L, 600L),
                                       450L, 580L)
  expect_equal(len, 151L)
})

test_that(".compute5UtrLength gives 0 for no 5'UTR", {
  len <- Isopair:::.compute5UtrLength("+", 100L, 300L, 100L, 280L)
  expect_equal(len, 0L)
})

test_that(".compute5UtrLength minus strand", {
  # Exon 100-300, CDS cds_start=120, cds_stop=250
  # 5'UTR = 300-250 = 50 bp
  len <- Isopair:::.compute5UtrLength("-", 100L, 300L, 120L, 250L)
  expect_equal(len, 50L)
})


# ==============================================================================
# .transcriptToGenomic()
# ==============================================================================

test_that(".transcriptToGenomic works on plus strand", {
  # Exon1: 100-109 (10bp), Exon2: 200-209 (10bp)
  # Transcript pos 1 → genomic 100, pos 10 → 109, pos 11 → 200, pos 20 → 209
  expect_equal(Isopair:::.transcriptToGenomic(1L, "+", c(100L, 200L),
                                               c(109L, 209L)), 100L)
  expect_equal(Isopair:::.transcriptToGenomic(10L, "+", c(100L, 200L),
                                                c(109L, 209L)), 109L)
  expect_equal(Isopair:::.transcriptToGenomic(11L, "+", c(100L, 200L),
                                                c(109L, 209L)), 200L)
  expect_equal(Isopair:::.transcriptToGenomic(20L, "+", c(100L, 200L),
                                                c(109L, 209L)), 209L)
})

test_that(".transcriptToGenomic returns NA beyond transcript", {
  expect_true(is.na(
    Isopair:::.transcriptToGenomic(21L, "+", c(100L, 200L), c(109L, 209L))
  ))
})

test_that(".transcriptToGenomic works on minus strand", {
  # Exon1: 100-109 (10bp), Exon2: 200-209 (10bp)
  # Minus strand: transcript 5' → 3' is genomic high → low
  # Transcript pos 1 → genomic 209, pos 10 → 200, pos 11 → 109, pos 20 → 100
  expect_equal(Isopair:::.transcriptToGenomic(1L, "-", c(100L, 200L),
                                               c(109L, 209L)), 209L)
  expect_equal(Isopair:::.transcriptToGenomic(10L, "-", c(100L, 200L),
                                                c(109L, 209L)), 200L)
  expect_equal(Isopair:::.transcriptToGenomic(11L, "-", c(100L, 200L),
                                                c(109L, 209L)), 109L)
  expect_equal(Isopair:::.transcriptToGenomic(20L, "-", c(100L, 200L),
                                                c(109L, 209L)), 100L)
})

test_that(".transcriptToGenomic returns NA for invalid input", {
  expect_true(is.na(Isopair:::.transcriptToGenomic(NA, "+", 100L, 200L)))
  expect_true(is.na(Isopair:::.transcriptToGenomic(0L, "+", 100L, 200L)))
})


# ==============================================================================
# .scoreKozak()
# ==============================================================================

test_that(".scoreKozak scores strong Kozak correctly", {
  # Strong: A/G at -3, G at +4
  # Sequence: ...A??ATG G...   pos of A in ATG = 4
  #                   1234567
  kz <- Isopair:::.scoreKozak("ACCATGG", 4L)
  expect_equal(kz$score, 2L)  # A at -3, G at +4
})

test_that(".scoreKozak scores weak Kozak correctly", {
  # Neither A/G at -3 nor G at +4
  kz <- Isopair:::.scoreKozak("TCCATGA", 4L)
  expect_equal(kz$score, 0L)
})

test_that(".scoreKozak returns NA at edge of sequence", {
  kz <- Isopair:::.scoreKozak("ATG", 1L)
  expect_true(is.na(kz$score))
})


# ==============================================================================
# scan5UtrFeatures()
# ==============================================================================

test_that("scan5UtrFeatures computes basic features correctly", {
  # Single-exon isoform, + strand
  # Exon: 100-130 (31 bp)
  # CDS: 110-130 (ATG at pos 110)
  # 5'UTR: 100-109 (10 bp, transcript positions 1-10)
  # Sequence: 10 bp 5'UTR + 21 bp CDS
  # Put ATG in 5'UTR at position 2, with in-frame stop at position 11 (in CDS)
  #           1234567890 1234567890 1
  utr5_seq <- "AATGTAAACT"  # ATG at pos 2, TAA at pos 5 → ORF = pos 2-7 (6 nt = 2 codons, too short)
  cds_seq  <- "ATGCCCCCCCCCCCCCCCCTGA"  # 21 bp, starts at transcript pos 11

  # Let's construct more carefully:
  # 5'UTR = 10bp, need ATG that produces a qualifying ORF (>= 9 nt)
  # Put ATG at pos 1, find stop at pos 1+9=10+ → need stop in CDS region
  #          1234567890 12345678901234567890 1
  full_seq <- "ATGAAAAAA" # 9bp so far
  # Need stop codon in-frame starting at position 10, 13, 16, ...
  # Position 10 is the start of CDS. Let's make the CDS:
  # pos 10 = start of CDS (last bp of 5'UTR is pos 9... wait.
  # Let me recalculate.
  # 5'UTR = positions 1-10. CDS starts at position 11.
  # ATG at pos 1 in 5'UTR. Walk in-frame: codons at 4, 7, 10, 13
  # Codon at pos 10-12 would be in CDS territory.

  # Let's just build a simple test case:
  # 5'UTR (12 bp): "AAATGAAATGA" + one more = "AAATGAAATGAA" (12bp)
  # Has ATG at pos 3 and another at pos 8
  # CDS (9 bp): "ATGCCCTGA"
  # Full: "AAATGAAATGAAATGCCCTGA" (21 bp total)
  # 5'UTR = 12 bp, CDS starts at pos 13

  # ATG at pos 3: walk in-frame → codons 6("AAA"), 9("TGA")... wait let me count:
  # pos 3-5 = "ATG", pos 6-8 = "AAA", pos 9-11 = "TGA"... no, let me recount
  # "AAATGAAATGAAATGCCCTGA"
  #  123456789012345678901
  # ATG at pos 3: next codons at 6,9,12,15,18
  # pos 6-8 = "AAA", pos 9-11 = "TGA" → not right, let me look at actual chars
  # A A A T G A A A T G A A A T G C C C T G A
  # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
  # ATG at pos 3? No, pos 3 = 'A'. ATG at pos 4: A=4,T=5,G=6
  # Actually let me use a cleaner approach

  full_seq <- "CCATGCCCCCCCCCTAACCATGCCCTGA"
  #             123456789012345678901234567890
  # C C A T G C C C C C C C C C T A A C C A T G C C C T G A
  # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8

  # 5'UTR = 18 bp (positions 1-18)
  # CDS starts at position 19 (pos 19-21 = "CCA"... hmm, need ATG at 19)

  # Simplify:
  full_seq <- "CCATGCCCCCCTAACCCCATGCCCTGA"
  #             1234567890123456789012345678
  # 5'UTR = 18 bp, CDS from position 19
  # But CDS must start with ATG at position 19
  # Let me just set utr5 to 15 bp:
  # full_seq = 15 bp UTR + "ATGCCCTGA" CDS (9bp) = 24 bp total
  full_seq <- "CCATGCCCCCCTAAGATGCCCCTATGA"
  #             123456789012345678901234567
  # C C A T G C C C C C C T A A G A T G C C C C T A T G A
  # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
  # Hmm this is getting messy. Let me use a programmatic approach.

  # Just make a clean example:
  # 5'UTR = 15bp, CDS = 12bp, total = 27bp
  utr5 <- "CCCATGAAACCCTAA"  # 15 bp, has ATG at pos 4
  cds  <- "ATGAAAAAATGA"      # 12 bp (ATG + 3 codons including stop)
  full_seq <- paste0(utr5, cds)
  # ATG at pos 4 in 5'UTR
  # Walk in-frame from ATG at 4: codons at 7("AAA"), 10("CCC"), 13("TAA") → STOP!
  # uORF = pos 4-15 = 12 nt (4 codons), entirely within 5'UTR
  # CDS ATG at pos 16 → frame = (16 - 4) %% 3 = 12 %% 3 = 0 → in-frame

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds <- make_cds("TX1", cds_start = 115L, cds_stop = 126L, strand = "+")
  # exon: 100-126, cds_start=115 means 5'UTR = 100-114 = 15 bp ✓

  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds, seqs, verbose = FALSE)
  expect_equal(nrow(result), 1)
  expect_equal(result$utr5_length, 15L)
  expect_equal(result$n_atg, 1L)  # one ATG in 5'UTR
  expect_equal(result$n_orfs, 1L)
  expect_equal(result$n_orfs_inframe, 1L)
  expect_equal(result$n_orfs_outframe, 0L)
  expect_equal(result$n_orfs_overlapping, 0L)
  expect_equal(result$longest_orf_nt, 12L)  # pos 4-15
  expect_true(result$atg_validated)
})

test_that("scan5UtrFeatures handles zero-length 5'UTR", {
  full_seq <- "ATGAAATGA"  # 9 bp CDS, no 5'UTR
  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L, exon_ends = 108L)
  cds <- make_cds("TX1", cds_start = 100L, cds_stop = 108L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds, seqs, verbose = FALSE)
  expect_equal(result$utr5_length, 0L)
  expect_equal(result$n_atg, 0L)
  expect_equal(result$n_orfs, 0L)
})

test_that("scan5UtrFeatures handles overlapping uORF", {
  # 5'UTR = 12 bp, CDS = 15 bp, total = 27 bp
  utr5 <- "CCCATGAAAAAA"  # 12 bp, ATG at pos 4
  cds  <- "ATGAAAAAAAAATGA"  # 15 bp
  full_seq <- paste0(utr5, cds)
  # ATG at pos 4: walk in-frame → codons at 7("AAA"), 10("AAA"), 13("ATG")... no stop in utr
  # Continue past CDS start (pos 13): 13("ATG"), 16("AAA"), 19("AAA"), 22("AAT"), 25("GA")...
  # Wait, let me recount:
  # C C C A T G A A A A A A A T G A A A A A A A A A T G A
  # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
  # ATG at 4: next codons at 7("AAA"), 10("AAA"), 13("ATG"), 16("AAA"), 19("AAA"), 22("AAT"), 25("GA")→ too short
  # No stop found → no uORF

  # Let me put an explicit stop:
  utr5 <- "CCCATGAAAAAA"  # 12 bp, ATG at pos 4
  cds  <- "ATGCCCTGACCCTGA"  # 15 bp
  full_seq <- paste0(utr5, cds)
  # C C C A T G A A A A A A A T G C C C T G A C C C T G A
  # 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
  # ATG at 4: codons at 7("AAA"), 10("AAA"), 13("ATG"), 16("CCC"), 19("TGA") → STOP at 19-21
  # uORF = pos 4-21 = 18 nt (6 codons)
  # CDS starts at pos 13, stop is at 21 → overlapping! (21 > 12)
  # Frame: (13 - 4) %% 3 = 9 %% 3 = 0 → in-frame

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 112L, cds_stop = 126L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(result$n_orfs, 1L)
  expect_equal(result$n_orfs_overlapping, 1L)
  expect_equal(result$longest_orf_nt, 18L)
})

test_that("scan5UtrFeatures handles multiple isoforms", {
  full1 <- "CCCATGAAACCCTAAATGCCCTGA"  # 15bp UTR + 9bp CDS
  full2 <- "ATGAAATGA"  # no UTR

  structures <- dplyr::bind_rows(
    make_structures("TX1", "GENE1", strand = "+", exon_starts = 100L, exon_ends = 123L),
    make_structures("TX2", "GENE2", strand = "+", exon_starts = 200L, exon_ends = 208L)
  )
  cds_meta <- dplyr::bind_rows(
    make_cds("TX1", cds_start = 115L, cds_stop = 123L, strand = "+"),
    make_cds("TX2", cds_start = 200L, cds_stop = 208L, strand = "+")
  )
  seqs <- c(TX1 = full1, TX2 = full2)

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(nrow(result), 2)
  expect_equal(result$isoform_id, c("TX1", "TX2"))
})

test_that("scan5UtrFeatures returns empty for no common IDs", {
  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L, exon_ends = 200L)
  cds_meta <- make_cds("TX2", cds_start = 100L, cds_stop = 200L, strand = "+")
  seqs <- c(TX3 = "ATGCCCTGA")

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(nrow(result), 0)
})

test_that("scan5UtrFeatures works on minus strand", {
  # Exon: 100-126 (27 bp)
  # Minus strand: CDS cds_start=100, cds_stop=114
  # Translation starts at cds_stop=114 (5' end)
  # 5'UTR = exonic bases > 114 → 115-126 = 12 bp

  # Transcript sequence is read high-to-low on minus strand:
  # Transcript pos 1 = genomic 126, pos 12 = genomic 115 (5'UTR)
  # Transcript pos 13 = genomic 114 (CDS start = ATG)

  utr5 <- "CCCATGAAAAAA"  # 12 bp 5'UTR
  cds  <- "ATGCCCCCCCCCTGA"  # 15 bp CDS
  full_seq <- paste0(utr5, cds)

  structures <- make_structures("TX1", "GENE1", strand = "-",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 100L, cds_stop = 114L, strand = "-")
  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(result$utr5_length, 12L)
  expect_equal(result$strand, "-")
})

test_that("scan5UtrFeatures computes frame-stratified counts", {
  # 5'UTR = 18 bp. CDS starts at pos 19.
  # ATG at pos 1: frame = (19-1) %% 3 = 0 → frame0
  # ATG at pos 5: frame = (19-5) %% 3 = 14 %% 3 = 2 → frame2
  # ATG at pos 10: frame = (19-10) %% 3 = 9 %% 3 = 0 → frame0
  utr5 <- "ATGCATGCCCCATGCCCCC"  # 19 chars... need exactly 18
  # Actually let me be precise:
  # pos: 1234567890123456789
  #      ATGCCATGCCCATGCCCCC
  # Hmm, need 18 bp. Let me just construct it:
  utr5 <- "ATGCCATGCCCATGCCCC"  # 18 bp
  # ATG at pos 1, 6, 12
  # Check: A(1)T(2)G(3)C(4)C(5)A(6)T(7)G(8)C(9)C(10)C(11)A(12)T(13)G(14)C(15)C(16)C(17)C(18)
  # ATG at pos 1, pos 6, pos 12
  # Frames: (19-1)%%3=0, (19-6)%%3=1, (19-12)%%3=1

  cds <- "ATGCCCTGA"  # 9 bp
  full_seq <- paste0(utr5, cds)

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 118L, cds_stop = 126L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(result$n_atg, 3L)
  expect_equal(result$n_atg_frame0, 1L)  # pos 1
  expect_equal(result$n_atg_frame1, 2L)  # pos 6, pos 12
  expect_equal(result$n_atg_frame2, 0L)
})

test_that("scan5UtrFeatures computes density features", {
  utr5 <- "CCCATGAAACCCTAA"  # 15 bp, 1 ATG, 1 TAA stop, 1 ORF
  cds  <- "ATGCCCTGA"         # 9 bp
  full_seq <- paste0(utr5, cds)

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 115L, cds_stop = 123L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- scan5UtrFeatures(structures, cds_meta, seqs, verbose = FALSE)
  expect_equal(result$atg_density, 1 / 15 * 100, tolerance = 1e-6)
  expect_true(result$stop_density > 0)
})


# ==============================================================================
# detectUorfs()
# ==============================================================================

test_that("detectUorfs finds uORFs correctly", {
  utr5 <- "CCCATGAAACCCTAA"  # 15 bp, ATG at pos 4
  cds  <- "ATGCCCTGA"
  full_seq <- paste0(utr5, cds)
  # ATG at pos 4: codons at 7("AAA"), 10("CCC"), 13("TAA") → STOP
  # uORF = pos 4-15 = 12 nt

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 115L, cds_stop = 123L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", kozak = TRUE, verbose = FALSE)

  expect_equal(nrow(result), 1)
  expect_equal(result$isoform_id, "TX1")
  expect_equal(result$atg_transcript_pos, 4L)
  expect_equal(result$stop_transcript_pos, 15L)
  expect_equal(result$uorf_length_nt, 12L)
  expect_false(result$is_overlapping)
  expect_equal(result$frame_relative_to_cds, 0L)  # (16-4)%%3=0
  expect_equal(result$stop_codon, "TAA")
  expect_equal(result$uorf_sequence, "ATGAAACCCTAA")
})

test_that("detectUorfs returns genomic coordinates", {
  utr5 <- "CCCATGAAACCCTAA"  # 15 bp
  cds  <- "ATGCCCTGA"
  full_seq <- paste0(utr5, cds)

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 115L, cds_stop = 123L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", kozak = FALSE, verbose = FALSE)

  # Transcript pos 4 → genomic 103 (100 + 4 - 1)
  # Transcript pos 15 → genomic 114 (100 + 15 - 1)
  expect_equal(result$atg_genomic_pos, 103L)
  expect_equal(result$stop_genomic_pos, 114L)
})

test_that("detectUorfs handles no uORFs", {
  full_seq <- "CCCCCCCCCCCCCCCATGCCCTGA"  # 15 bp UTR with no ATG
  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 115L, cds_stop = 123L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", verbose = FALSE)
  expect_equal(nrow(result), 0)
})

test_that("detectUorfs handles overlapping uORF", {
  utr5 <- "CCCATGAAAAAA"  # 12 bp, ATG at pos 4
  cds  <- "ATGCCCTGACCCTGA"  # 15 bp
  full_seq <- paste0(utr5, cds)
  # ATG at 4: codons 7("AAA"),10("AAA"),13("ATG"),16("CCC"),19("TGA") → stop at 19-21
  # uORF = 4-21 = 18 nt, overlapping (21 > 12)

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 112L, cds_stop = 126L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", verbose = FALSE)
  expect_equal(nrow(result), 1)
  expect_true(result$is_overlapping)
})

test_that("detectUorfs filters by minimum length", {
  # ATG at pos 4, stop TAA at pos 7 → ORF = 6 nt (2 codons) → too short (min=3)
  utr5 <- "CCCATGTAAAAAAAAAA"  # 17 bp
  cds  <- "ATGCCCTGA"
  full_seq <- paste0(utr5, cds)

  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L,
                                exon_ends = 100L + nchar(full_seq) - 1L)
  cds_meta <- make_cds("TX1", cds_start = 117L, cds_stop = 125L, strand = "+")
  seqs <- setNames(full_seq, "TX1")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", min_length_codons = 3L,
                         verbose = FALSE)
  expect_equal(nrow(result), 0)
})

test_that("detectUorfs returns empty tibble schema on empty input", {
  structures <- make_structures("TX1", "GENE1", strand = "+",
                                exon_starts = 100L, exon_ends = 200L)
  cds_meta <- make_cds("TX2", cds_start = 100L, cds_stop = 200L, strand = "+")
  seqs <- c(TX3 = "ATGCCCTGA")

  result <- detectUorfs(structures, cds_meta, seqs,
                         method = "simple", verbose = FALSE)
  expect_equal(nrow(result), 0)
  expect_true(all(c("isoform_id", "gene_id", "uorf_id", "atg_genomic_pos",
                     "uorf_length_nt", "is_overlapping",
                     "frame_relative_to_cds") %in% names(result)))
})


# ==============================================================================
# compareUorfs()
# ==============================================================================

test_that("compareUorfs detects shared, gained, and lost uORFs", {
  # Two isoforms with uORFs at specific genomic positions
  uorf_table <- tibble::tibble(
    isoform_id = c("REF", "REF", "COMP", "COMP"),
    gene_id = "GENE1",
    uorf_id = c("uORF_1", "uORF_2", "uORF_1", "uORF_2"),
    atg_transcript_pos = c(4L, 10L, 4L, 20L),
    stop_transcript_pos = c(15L, 21L, 15L, 31L),
    atg_genomic_pos = c(103L, 109L, 103L, 119L),
    stop_genomic_pos = c(114L, 120L, 114L, 130L),
    uorf_length_nt = c(12L, 12L, 12L, 12L),
    is_overlapping = FALSE,
    frame_relative_to_cds = 0L,
    kozak_score = NA_integer_,
    kozak_context = NA_character_,
    uorf_sequence = "ATGAAACCCTAA",
    stop_codon = "TAA"
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1",
    reference_isoform_id = "REF",
    comparator_isoform_id = "COMP"
  )

  result <- compareUorfs(profiles, uorf_table)

  expect_equal(result$pair_summary$n_shared, 1L)   # ATG at 103
  expect_equal(result$pair_summary$n_lost, 1L)     # ATG at 109 (ref only)
  expect_equal(result$pair_summary$n_gained, 1L)   # ATG at 119 (comp only)

  detail <- result$uorf_detail
  expect_equal(sum(detail$status == "shared"), 1)
  expect_equal(sum(detail$status == "lost"), 1)
  expect_equal(sum(detail$status == "gained"), 1)
})

test_that("compareUorfs detects stop_changed in shared uORFs", {
  uorf_table <- tibble::tibble(
    isoform_id = c("REF", "COMP"),
    gene_id = "GENE1",
    uorf_id = "uORF_1",
    atg_transcript_pos = 4L,
    stop_transcript_pos = c(15L, 21L),
    atg_genomic_pos = 103L,
    stop_genomic_pos = c(114L, 120L),  # Different stops!
    uorf_length_nt = c(12L, 18L),
    is_overlapping = FALSE,
    frame_relative_to_cds = 0L,
    kozak_score = NA_integer_,
    kozak_context = NA_character_,
    uorf_sequence = "ATGAAACCCTAA",
    stop_codon = "TAA"
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1",
    reference_isoform_id = "REF",
    comparator_isoform_id = "COMP"
  )

  result <- compareUorfs(profiles, uorf_table)
  expect_equal(result$pair_summary$n_shared, 1L)
  expect_equal(result$pair_summary$n_shared_diff_stop, 1L)
  expect_true(result$uorf_detail$stop_changed[1])
})

test_that("compareUorfs handles pairs with no uORFs", {
  uorf_table <- tibble::tibble(
    isoform_id = character(0), gene_id = character(0),
    uorf_id = character(0),
    atg_transcript_pos = integer(0), stop_transcript_pos = integer(0),
    atg_genomic_pos = integer(0), stop_genomic_pos = integer(0),
    uorf_length_nt = integer(0), is_overlapping = logical(0),
    frame_relative_to_cds = integer(0),
    kozak_score = integer(0), kozak_context = character(0),
    uorf_sequence = character(0), stop_codon = character(0)
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1",
    reference_isoform_id = "REF",
    comparator_isoform_id = "COMP"
  )

  result <- compareUorfs(profiles, uorf_table)
  expect_equal(result$pair_summary$n_shared, 0L)
  expect_equal(result$pair_summary$n_gained, 0L)
  expect_equal(result$pair_summary$n_lost, 0L)
})

test_that("compareUorfs returns empty for empty profiles", {
  uorf_table <- tibble::tibble(
    isoform_id = character(0), gene_id = character(0),
    uorf_id = character(0),
    atg_transcript_pos = integer(0), stop_transcript_pos = integer(0),
    atg_genomic_pos = integer(0), stop_genomic_pos = integer(0),
    uorf_length_nt = integer(0), is_overlapping = logical(0),
    frame_relative_to_cds = integer(0),
    kozak_score = integer(0), kozak_context = character(0),
    uorf_sequence = character(0), stop_codon = character(0)
  )
  profiles <- tibble::tibble(
    gene_id = character(0),
    reference_isoform_id = character(0),
    comparator_isoform_id = character(0)
  )

  result <- compareUorfs(profiles, uorf_table)
  expect_equal(nrow(result$pair_summary), 0)
  expect_equal(nrow(result$uorf_detail), 0)
})


# ==============================================================================
# attributeUorfsToEvents()
# ==============================================================================

test_that("attributeUorfsToEvents handles frameshift attribution", {
  uorf_comparison <- list(
    pair_summary = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      n_ref_uorfs = 1L, n_comp_uorfs = 1L,
      n_shared = 1L, n_gained = 0L, n_lost = 0L,
      n_shared_same_stop = 0L, n_shared_diff_stop = 1L
    ),
    uorf_detail = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      atg_genomic_pos = 103L,
      status = "shared",
      ref_stop_genomic_pos = 114L, comp_stop_genomic_pos = 120L,
      ref_uorf_length_nt = 12L, comp_uorf_length_nt = 18L,
      ref_is_overlapping = FALSE, comp_is_overlapping = FALSE,
      ref_frame_to_cds = 0L, comp_frame_to_cds = 0L,
      stop_changed = TRUE
    )
  )

  # Event between ATG (103) and stop (114/120) with non-3 bp_diff
  events <- data.frame(
    event_type = "alternative_5prime",
    direction = "LOSS",
    five_prime = 108L,
    three_prime = 112L,
    bp_diff = 5L,  # Not divisible by 3 → frameshift
    stringsAsFactors = FALSE
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1", reference_isoform_id = "REF",
    comparator_isoform_id = "COMP",
    detailed_events = list(events)
  )

  cds_meta <- make_cds("REF", cds_start = 115L, cds_stop = 130L, strand = "+")

  result <- attributeUorfsToEvents(profiles, uorf_comparison, cds_meta)
  expect_equal(nrow(result), 1)
  expect_equal(result$attribution_type, "frameshift")
  expect_equal(result$event_bp_diff, 5L)
})

test_that("attributeUorfsToEvents handles containment for gained uORF", {
  uorf_comparison <- list(
    pair_summary = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      n_ref_uorfs = 0L, n_comp_uorfs = 1L,
      n_shared = 0L, n_gained = 1L, n_lost = 0L,
      n_shared_same_stop = 0L, n_shared_diff_stop = 0L
    ),
    uorf_detail = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      atg_genomic_pos = 110L,
      status = "gained",
      ref_stop_genomic_pos = NA_integer_,
      comp_stop_genomic_pos = 120L,
      ref_uorf_length_nt = NA_integer_,
      comp_uorf_length_nt = 12L,
      ref_is_overlapping = NA,
      comp_is_overlapping = FALSE,
      ref_frame_to_cds = NA_integer_,
      comp_frame_to_cds = 0L,
      stop_changed = NA
    )
  )

  events <- data.frame(
    event_type = "skipped_exon",
    direction = "GAIN",
    five_prime = 105L,
    three_prime = 115L,
    bp_diff = 11L,  # Contains ATG at 110
    stringsAsFactors = FALSE
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1", reference_isoform_id = "REF",
    comparator_isoform_id = "COMP",
    detailed_events = list(events)
  )

  cds_meta <- make_cds("REF", cds_start = 125L, cds_stop = 140L, strand = "+")

  result <- attributeUorfsToEvents(profiles, uorf_comparison, cds_meta)
  expect_equal(nrow(result), 1)
  expect_equal(result$attribution_type, "containment")
})

test_that("attributeUorfsToEvents handles inframe_splice", {
  uorf_comparison <- list(
    pair_summary = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      n_ref_uorfs = 1L, n_comp_uorfs = 1L,
      n_shared = 1L, n_gained = 0L, n_lost = 0L,
      n_shared_same_stop = 0L, n_shared_diff_stop = 1L
    ),
    uorf_detail = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      atg_genomic_pos = 103L,
      status = "shared",
      ref_stop_genomic_pos = 114L, comp_stop_genomic_pos = 120L,
      ref_uorf_length_nt = 12L, comp_uorf_length_nt = 18L,
      ref_is_overlapping = FALSE, comp_is_overlapping = FALSE,
      ref_frame_to_cds = 0L, comp_frame_to_cds = 0L,
      stop_changed = TRUE
    )
  )

  events <- data.frame(
    event_type = "alternative_5prime",
    direction = "LOSS",
    five_prime = 108L,
    three_prime = 112L,
    bp_diff = 6L,  # Divisible by 3 → in-frame
    stringsAsFactors = FALSE
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1", reference_isoform_id = "REF",
    comparator_isoform_id = "COMP",
    detailed_events = list(events)
  )

  cds_meta <- make_cds("REF", cds_start = 125L, cds_stop = 140L, strand = "+")

  result <- attributeUorfsToEvents(profiles, uorf_comparison, cds_meta)
  expect_equal(nrow(result), 1)
  expect_equal(result$attribution_type, "inframe_splice")
})

test_that("attributeUorfsToEvents returns unresolved when no events match", {
  uorf_comparison <- list(
    pair_summary = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      n_ref_uorfs = 0L, n_comp_uorfs = 1L,
      n_shared = 0L, n_gained = 1L, n_lost = 0L,
      n_shared_same_stop = 0L, n_shared_diff_stop = 0L
    ),
    uorf_detail = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      atg_genomic_pos = 110L,
      status = "gained",
      ref_stop_genomic_pos = NA_integer_,
      comp_stop_genomic_pos = 120L,
      ref_uorf_length_nt = NA_integer_,
      comp_uorf_length_nt = 12L,
      ref_is_overlapping = NA,
      comp_is_overlapping = FALSE,
      ref_frame_to_cds = NA_integer_,
      comp_frame_to_cds = 0L,
      stop_changed = NA
    )
  )

  # Event far from the ATG
  events <- data.frame(
    event_type = "skipped_exon",
    direction = "GAIN",
    five_prime = 200L,
    three_prime = 210L,
    bp_diff = 11L,
    stringsAsFactors = FALSE
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1", reference_isoform_id = "REF",
    comparator_isoform_id = "COMP",
    detailed_events = list(events)
  )

  cds_meta <- make_cds("REF", cds_start = 125L, cds_stop = 140L, strand = "+")

  result <- attributeUorfsToEvents(profiles, uorf_comparison, cds_meta)
  expect_equal(result$attribution_type, "unresolved")
})

test_that("attributeUorfsToEvents returns empty for no actionable uORFs", {
  uorf_comparison <- list(
    pair_summary = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      n_ref_uorfs = 1L, n_comp_uorfs = 1L,
      n_shared = 1L, n_gained = 0L, n_lost = 0L,
      n_shared_same_stop = 1L, n_shared_diff_stop = 0L
    ),
    uorf_detail = tibble::tibble(
      gene_id = "GENE1", reference_isoform_id = "REF",
      comparator_isoform_id = "COMP",
      atg_genomic_pos = 103L,
      status = "shared",
      ref_stop_genomic_pos = 114L, comp_stop_genomic_pos = 114L,
      ref_uorf_length_nt = 12L, comp_uorf_length_nt = 12L,
      ref_is_overlapping = FALSE, comp_is_overlapping = FALSE,
      ref_frame_to_cds = 0L, comp_frame_to_cds = 0L,
      stop_changed = FALSE  # Same stop → not actionable
    )
  )

  profiles <- tibble::tibble(
    gene_id = "GENE1", reference_isoform_id = "REF",
    comparator_isoform_id = "COMP",
    detailed_events = list(data.frame(
      event_type = "alt_5", direction = "LOSS",
      five_prime = 108L, three_prime = 112L, bp_diff = 5L,
      stringsAsFactors = FALSE
    ))
  )

  cds_meta <- make_cds("REF", cds_start = 125L, cds_stop = 140L, strand = "+")

  result <- attributeUorfsToEvents(profiles, uorf_comparison, cds_meta)
  expect_equal(nrow(result), 0)
})
