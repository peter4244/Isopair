test_that("selectPrimaryOrf prefers the dominant isoform's ATG when exonic", {
  # Construct two isoforms. The dominant's CDS ATG at tx pos 30 is exonic
  # in both isoforms. Each has a second strong-Kozak ATG upstream (at 10)
  # that selectPrimaryOrf should NOT prefer when step 1 succeeds.
  len <- 200L
  s <- rep("A", len)
  # ATG at 10 (different Kozak context) with stop at 40 (30 nt ORF)
  s[10:12] <- c("A","T","G"); s[40:42] <- c("T","A","A")
  for (p in seq(13L, 37L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  # ATG at 30 with stop at 150 -- "dominant" ATG, 120 nt ORF
  s[30:32] <- c("A","T","G"); s[150:152] <- c("T","A","A")
  for (p in seq(33L, 147L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  seq_str <- paste0(s, collapse = "")

  # Both isoforms have the same exon structure so the dominant ATG at
  # genomic pos 30 is exonic in both. Kozak filter off for the test.
  structures <- tibble::tibble(
    isoform_id = c("dom", "other"),
    gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = len,
    exon_starts = list(1L, 1L), exon_ends = list(len, len)
  )
  cds <- tibble::tibble(
    isoform_id = c("dom", "other"),
    coding_status = "coding", strand = "+",
    cds_start = c(30L, 30L), cds_stop = c(150L, 150L),
    orf_length = c(120L, 120L))
  sequences <- c(dom = seq_str, other = seq_str)

  orfs <- enumerateOrfs(structures, cds, sequences, kozak_filter = FALSE)
  prim <- selectPrimaryOrf(orfs, structures, cds,
                           dominant_isoform_id = "dom")
  expect_equal(nrow(prim), 2L)
  expect_equal(prim$primary_atg_tx_pos[prim$isoform_id == "dom"],   30L)
  expect_equal(prim$primary_atg_tx_pos[prim$isoform_id == "other"], 30L)
  expect_true(all(prim$primary_source == "dominant"))
  expect_true(all(prim$primary_orf_length == 120L))
})

test_that("selectPrimaryOrf falls back to 5'-most ATG when dominant is not exonic", {
  # Target isoform is missing the exon containing the dominant's ATG at
  # genomic pos 50 -> step 1 fails, step 2 picks the 5'-most enumerated
  # ATG in the target.
  len <- 300L
  s <- rep("A", len)
  # first ATG at 120 (within the target's exon), stop at 150 (30 nt)
  s[120:122] <- c("A","T","G"); s[150:152] <- c("T","A","A")
  for (p in seq(123L, 147L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  # second ATG at 200, stop at 230 (30 nt)
  s[200:202] <- c("A","T","G"); s[230:232] <- c("T","A","A")
  for (p in seq(203L, 227L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  seq_str <- paste0(s, collapse = "")

  structures <- tibble::tibble(
    isoform_id = c("dom", "target"),
    gene_id = "g1", strand = "+",
    n_exons = c(2L, 1L),
    tx_start = c(1L, 101L), tx_end = c(300L, 300L),
    exon_starts = list(c(1L, 101L), 101L),
    exon_ends   = list(c(100L, 300L), 300L))
  cds <- tibble::tibble(
    isoform_id = c("dom", "target"),
    coding_status = c("coding", "non-coding"), strand = "+",
    cds_start = c(50L, NA_integer_), cds_stop = c(250L, NA_integer_),
    orf_length = c(200L, NA_integer_))
  sequences <- c(dom = seq_str, target = seq_str)

  orfs <- enumerateOrfs(structures, cds, sequences, kozak_filter = FALSE)
  prim <- selectPrimaryOrf(orfs, structures, cds,
                           dominant_isoform_id = "dom")
  t_row <- prim[prim$isoform_id == "target", ]
  expect_equal(t_row$primary_source, "first_hq")
  expect_equal(t_row$primary_atg_tx_pos, 120L)    # 5'-most ATG in target
})

test_that("selectPrimaryOrf with NULL dominant_isoform_id uses first_hq everywhere", {
  len <- 100L
  s <- rep("A", len)
  s[20:22] <- c("A","T","G"); s[50:52] <- c("T","A","A")
  for (p in seq(23L, 47L, by = 3L)) s[p:(p + 2L)] <- c("T","C","A")
  seq_str <- paste0(s, collapse = "")
  structures <- tibble::tibble(
    isoform_id = "iso", gene_id = "g1", strand = "+", n_exons = 1L,
    tx_start = 1L, tx_end = len,
    exon_starts = list(1L), exon_ends = list(len)
  )
  cds <- tibble::tibble(isoform_id = "iso", coding_status = "non-coding",
                       strand = "+", cds_start = NA_integer_,
                       cds_stop = NA_integer_, orf_length = NA_integer_)
  orfs <- enumerateOrfs(structures, cds, c(iso = seq_str), kozak_filter = FALSE)
  prim <- selectPrimaryOrf(orfs, structures, cds, dominant_isoform_id = NULL)
  expect_equal(prim$primary_source, "first_hq")
  expect_equal(prim$primary_atg_tx_pos, 20L)
})

test_that("defaultKozakThreshold loads values from inst/extdata", {
  t5 <- defaultKozakThreshold(0.05)
  expect_true(is.numeric(t5))
  expect_true(is.finite(t5))
  # MANE 5th-percentile on GENCODE v49 should be around -1.25
  expect_lt(t5, 0)
  expect_gt(t5, -3)
})
