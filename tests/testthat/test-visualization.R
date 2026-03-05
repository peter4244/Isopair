# Smoke tests for plotIsoformPair()

test_that("plotIsoformPair: returns ggplot object", {
  skip_if_not_installed("ggplot2")

  # Pick the first test pair
  pair <- test_pairs_raw[1, ]
  ref_exons <- test_exons[test_exons$isoform_id == pair$isoform_A, ]
  comp_exons <- test_exons[test_exons$isoform_id == pair$isoform_B, ]
  gene_strand <- unique(ref_exons$strand)[1]

  ref_df <- data.frame(
    exon_start = ref_exons$exon_start, exon_end = ref_exons$exon_end
  )
  comp_df <- data.frame(
    exon_start = comp_exons$exon_start, exon_end = comp_exons$exon_end
  )

  events <- detectEvents(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand
  )

  p <- plotIsoformPair(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    events = events,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand
  )

  expect_s3_class(p, "ggplot")
})

test_that("plotIsoformPair: works with reconstructed track", {
  skip_if_not_installed("ggplot2")

  pair <- test_pairs_raw[1, ]
  ref_exons <- test_exons[test_exons$isoform_id == pair$isoform_A, ]
  comp_exons <- test_exons[test_exons$isoform_id == pair$isoform_B, ]
  gene_strand <- unique(ref_exons$strand)[1]

  ref_df <- data.frame(
    exon_start = ref_exons$exon_start, exon_end = ref_exons$exon_end
  )
  comp_df <- data.frame(
    exon_start = comp_exons$exon_start, exon_end = comp_exons$exon_end
  )

  events <- detectEvents(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand
  )

  recon <- reconstructDominant(comp_df, events)

  p <- plotIsoformPair(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    events = events,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand,
    reconstructed_exons = recon
  )

  expect_s3_class(p, "ggplot")
})

test_that("plotIsoformPair: works with CDS metadata", {
  skip_if_not_installed("ggplot2")

  pair <- test_pairs_raw[1, ]
  ref_exons <- test_exons[test_exons$isoform_id == pair$isoform_A, ]
  comp_exons <- test_exons[test_exons$isoform_id == pair$isoform_B, ]
  gene_strand <- unique(ref_exons$strand)[1]

  ref_df <- data.frame(
    exon_start = ref_exons$exon_start, exon_end = ref_exons$exon_end
  )
  comp_df <- data.frame(
    exon_start = comp_exons$exon_start, exon_end = comp_exons$exon_end
  )

  events <- detectEvents(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand
  )

  cds <- extractCdsAnnotations(test_gtf_path, verbose = FALSE)

  p <- plotIsoformPair(
    reference_exons = ref_df,
    comparator_exons = comp_df,
    events = events,
    gene_id = pair$gene_id,
    reference_id = pair$isoform_A,
    comparator_id = pair$isoform_B,
    strand = gene_strand,
    cds_metadata = cds
  )

  expect_s3_class(p, "ggplot")
})

test_that("plotIsoformPair: works with no events", {
  skip_if_not_installed("ggplot2")

  exons <- data.frame(
    exon_start = c(1000L, 2000L, 3000L),
    exon_end = c(1500L, 2500L, 3500L)
  )

  events <- detectEvents(
    reference_exons = exons,
    comparator_exons = exons,
    gene_id = "GENE_IDENT",
    reference_id = "REF_IDENT",
    comparator_id = "COMP_IDENT",
    strand = "+"
  )

  p <- plotIsoformPair(
    reference_exons = exons,
    comparator_exons = exons,
    events = events,
    gene_id = "GENE_IDENT",
    reference_id = "REF_IDENT",
    comparator_id = "COMP_IDENT",
    strand = "+"
  )

  expect_s3_class(p, "ggplot")
})
