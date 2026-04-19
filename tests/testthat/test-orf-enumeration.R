test_that("empiricalKozakThreshold returns a data-driven cutoff", {
  # Two coding isoforms with annotated CDS ATGs in different Kozak contexts.
  seq_strong <- paste0(
    paste0(rep("A", 18), collapse = ""),
    "GCCACC", "ATG", "GC",                     # strong Kozak ATG at pos 25
    paste0(rep(c("T","C","A"), 20), collapse = ""),
    "TAA"
  )
  seq_medium <- paste0(
    paste0(rep("A", 18), collapse = ""),
    "AAAGCC", "ATG", "AA",                     # weaker Kozak at pos 25
    paste0(rep(c("T","C","A"), 20), collapse = ""),
    "TAA"
  )
  len_s <- nchar(seq_strong); len_m <- nchar(seq_medium)
  structures <- tibble::tibble(
    isoform_id = c("iso_strong", "iso_medium"),
    gene_id    = "g1", strand = "+", n_exons = 1L,
    tx_start   = 1L, tx_end = c(len_s, len_m),
    exon_starts = list(1L, 1L),
    exon_ends   = list(len_s, len_m)
  )
  cds <- tibble::tibble(
    isoform_id = c("iso_strong", "iso_medium"),
    coding_status = "coding", strand = "+",
    cds_start = c(25L, 25L), cds_stop = c(100L, 100L),
    orf_length = c(75L, 75L)
  )
  seqs <- c(iso_strong = seq_strong, iso_medium = seq_medium)

  emp <- empiricalKozakThreshold(structures, cds, seqs, quantile = 0.05)
  expect_equal(emp$n_used, 2L)
  expect_true(is.finite(emp$threshold))
  expect_true(all(is.finite(emp$scores)))
  expect_length(emp$scores, 2L)
  expect_true(emp$threshold >= min(emp$scores) &&
              emp$threshold <= max(emp$scores))
  # Strong Kozak context should score higher than the medium one.
  expect_gt(emp$scores[1], emp$scores[2])
})

test_that("enumerateOrfs filters by Kozak by default", {
  # Construct a transcript with two ATGs: one in strong Kozak context
  # (GCCACCATGG*) and one in weak context (TTTTTT ATG*). The strong one
  # should pass the default threshold (0); the weak one should not.
  #
  # Layout (1-based):
  #   1..10    padding
  #   11..18   strong context lead-in  (positions -7..-1 before ATG @ 25)
  #   19..20   "GG"  (+1..+2 after ATG at pos 19-21 actually — recompute)
  #
  # Simpler: build it position-by-position.
  len <- 200L
  s <- rep("A", len)
  # Strong Kozak ATG at position 25:
  #   pos 19..24 = "GCCACC", ATG at 25..27, then +4 = G, +5 = C
  s[19:27] <- c("G","C","C","A","C","C","A","T","G")
  s[28:29] <- c("G","C")
  # Seed a stop codon in frame for the strong ATG: ATG at 25 -> stop at 58+
  s[58:60] <- c("T","A","A")
  for (p in seq(28L, 55L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  # Re-place the + context after touching +4/+5
  s[28:29] <- c("G","C")
  # Weak Kozak ATG at position 100:
  #   pos 94..99 = "TTTTTT", ATG at 100..102, then +4 = T, +5 = T
  s[94:102] <- c("T","T","T","T","T","T","A","T","G")
  s[103:104] <- c("T","T")
  # Stop in frame for weak ATG at ~130
  s[130:132] <- c("T","A","A")
  for (p in seq(103L, 127L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  s[103:104] <- c("T","T")  # preserve weak +4/+5

  seq_str <- paste0(s, collapse = "")
  structures <- tibble::tibble(
    isoform_id = "iso1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = len,
    exon_starts = list(1L), exon_ends = list(len)
  )
  cds <- tibble::tibble(isoform_id = "iso1", coding_status = "non-coding",
                       strand = "+", cds_start = NA_integer_,
                       cds_stop = NA_integer_, orf_length = NA_integer_)
  sequences <- c(iso1 = seq_str)

  # Default (kozak_filter = TRUE): only the strong-Kozak ATG at pos 25
  default_out <- enumerateOrfs(structures, cds, sequences)
  expect_true(25L %in% default_out$atg_tx_pos)
  expect_false(100L %in% default_out$atg_tx_pos)
  expect_true(all(default_out$kozak_score >= 0))

  # Disable filter -> weak ATG also appears
  unfiltered <- enumerateOrfs(structures, cds, sequences,
                               kozak_filter = FALSE)
  expect_true(25L %in% unfiltered$atg_tx_pos)
  expect_true(100L %in% unfiltered$atg_tx_pos)
  # Stricter threshold -> both filtered out (weak definitely; strong
  # depends on exact score; tighten to a level above its score)
  strict <- enumerateOrfs(structures, cds, sequences,
                          kozak_threshold = 5)
  expect_equal(nrow(strict), 0L)
})

test_that("enumerateOrfs emits per-ORF rows for a two-ORF transcript", {
  # Single-isoform test: construct a 300-nt transcript with two viable ORFs.
  #   ATG at pos 10, stop at pos 40   -> 30-nt ORF, stop inside exon 1
  #     (exon 1 is 1-100, junction at 100; stop_end=42; 42+49 = 91 < 100
  #      so junction 100 is > 91 => downstream -> effectively_ptc)
  #   ATG at pos 150 (exon 2), stop at pos 210 -> 60-nt ORF in last exon
  #     (exon 2 is 101-300; last exon; no downstream junction -> no_downstream_ejc)
  len <- 300L
  s <- rep("A", len)
  s[10:12] <- c("A","T","G"); s[40:42] <- c("T","A","A")
  for (p in seq(13L, 37L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  s[150:152] <- c("A","T","G"); s[210:212] <- c("T","A","A")
  for (p in seq(153L, 207L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  seq_str <- paste0(s, collapse = "")

  structures <- tibble::tibble(
    isoform_id = "iso1", gene_id = "g1", strand = "+", n_exons = 2L,
    tx_start = 1L, tx_end = 400L,
    exon_starts = list(c(1L, 201L)), exon_ends = list(c(100L, 400L))
  )
  cds <- tibble::tibble(
    isoform_id = "iso1", coding_status = "coding", strand = "+",
    cds_start = 50L, cds_stop = 250L, orf_length = 67L
  )
  sequences <- c(iso1 = seq_str)

  out <- enumerateOrfs(structures, cds, sequences, kozak_filter = FALSE)
  expect_equal(nrow(out), 2L)
  expect_true(all(out$isoform_id == "iso1"))
  expect_setequal(out$atg_tx_pos, c(10L, 150L))

  r1 <- out[out$atg_tx_pos == 10L, ]
  expect_equal(r1$category, "effectively_ptc")
  expect_equal(r1$orf_length, 30L)
  expect_gt(r1$n_downstream_ejc, 0L)

  r2 <- out[out$atg_tx_pos == 150L, ]
  expect_equal(r2$category, "no_downstream_ejc")
  expect_equal(r2$orf_length, 60L)
  expect_equal(r2$n_downstream_ejc, 0L)
})

test_that("enumerateOrfs filters ORFs below min_orf_nt", {
  len <- 50L
  s <- rep("A", len)
  s[5:7]  <- c("A","T","G"); s[11:13] <- c("T","A","A")  # 6-nt ORF (below 30)
  s[20:22] <- c("A","T","G"); s[50:50] <- "A"            # runs off
  seq_str <- paste0(s, collapse = "")
  structures <- tibble::tibble(
    isoform_id = "iso1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = 50L,
    exon_starts = list(1L), exon_ends = list(50L)
  )
  cds <- tibble::tibble(isoform_id = "iso1", coding_status = "non-coding",
                       strand = "+", cds_start = NA_integer_,
                       cds_stop = NA_integer_, orf_length = NA_integer_)
  sequences <- c(iso1 = seq_str)
  out <- enumerateOrfs(structures, cds, sequences, include_no_stop = FALSE,
                       kozak_filter = FALSE)
  # The 6-nt ORF is below 30 and runoff is excluded -> 0 rows
  expect_equal(nrow(out), 0L)
})

test_that("enumerateOrfs emits 'no_stop_in_frame' when requested", {
  seq_str <- paste0(c(rep("C", 9), "A", "T", "G",
                      rep(c("T", "C", "A"), 20L),
                      rep("C", 9)), collapse = "")  # starts have ATG at 10, no stops
  structures <- tibble::tibble(
    isoform_id = "iso1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = nchar(seq_str),
    exon_starts = list(1L), exon_ends = list(nchar(seq_str))
  )
  cds <- tibble::tibble(isoform_id = "iso1", coding_status = "non-coding",
                       strand = "+", cds_start = NA_integer_,
                       cds_stop = NA_integer_, orf_length = NA_integer_)
  out <- enumerateOrfs(structures, cds, c(iso1 = seq_str),
                       include_no_stop = TRUE, kozak_filter = FALSE)
  expect_true(any(out$category == "no_stop_in_frame"))
  expect_true(all(is.na(out$orf_length[out$category == "no_stop_in_frame"])))
  expect_true(all(is.na(out$n_downstream_ejc[out$category == "no_stop_in_frame"])))
})

test_that("enumerateOrfs flags the annotated-CDS ATG via is_annotated_cds", {
  len <- 200L
  s <- rep("A", len)
  s[50:52] <- c("A","T","G")     # ATG at tx pos 50 (weak Kozak context)
  s[100:102] <- c("T","A","A")   # stop at 100 -> 50-nt ORF
  for (p in seq(53L, 97L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  seq_str <- paste0(s, collapse = "")
  # Single-exon isoform: tx pos 50 maps to genomic 50 on + strand.
  structures <- tibble::tibble(
    isoform_id = "iso1", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = len,
    exon_starts = list(1L), exon_ends = list(len)
  )
  cds <- tibble::tibble(
    isoform_id = "iso1", coding_status = "coding", strand = "+",
    cds_start = 50L, cds_stop = 100L, orf_length = 50L
  )
  # Kozak context here is weak; disable filter so the test focuses on the
  # is_annotated_cds flag rather than Kozak gating.
  out <- enumerateOrfs(structures, cds, c(iso1 = seq_str),
                       kozak_filter = FALSE)
  expect_true(any(out$is_annotated_cds))
  # Exactly one row should be the annotated CDS
  expect_equal(sum(out$is_annotated_cds), 1L)
  expect_equal(out$atg_tx_pos[out$is_annotated_cds], 50L)
})
