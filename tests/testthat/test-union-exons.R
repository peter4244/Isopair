# Tests for union exon construction

test_that("union exons have no overlaps within a gene", {
  genes <- unique(test_union_exons$gene_id)
  # Check a subset of genes
  check_genes <- head(genes, min(50, length(genes)))

  for (gene in check_genes) {
    gene_ues <- test_union_exons[test_union_exons$gene_id == gene, ]
    gene_ues <- gene_ues[order(gene_ues$start), ]

    if (nrow(gene_ues) < 2) next

    for (i in seq_len(nrow(gene_ues) - 1)) {
      expect_true(
        gene_ues$end[i] < gene_ues$start[i + 1],
        info = sprintf("Gene %s: UE %d [%d-%d] overlaps UE %d [%d-%d]",
                        gene, i, gene_ues$start[i], gene_ues$end[i],
                        i + 1, gene_ues$start[i + 1], gene_ues$end[i + 1])
      )
    }
  }
})

test_that("all isoform exon boundaries land on UE boundaries", {
  genes <- unique(test_union_exons$gene_id)
  check_genes <- head(genes, min(20, length(genes)))

  for (gene in check_genes) {
    gene_ues <- test_union_exons[test_union_exons$gene_id == gene, ]
    gene_isos <- test_structures[test_structures$gene_id == gene, ]

    ue_starts <- gene_ues$start
    ue_ends <- gene_ues$end

    for (j in seq_len(nrow(gene_isos))) {
      iso <- gene_isos[j, ]
      starts <- iso$exon_starts[[1]]
      ends <- iso$exon_ends[[1]]

      for (k in seq_along(starts)) {
        expect_true(
          starts[k] %in% ue_starts,
          info = sprintf("Gene %s, isoform %s: exon start %d not a UE boundary",
                          gene, iso$isoform_id, starts[k])
        )
        expect_true(
          ends[k] %in% ue_ends,
          info = sprintf("Gene %s, isoform %s: exon end %d not a UE boundary",
                          gene, iso$isoform_id, ends[k])
        )
      }
    }
  }
})

test_that("buildUnionExons returns correct structure", {
  result <- buildUnionExons(test_structures, verbose = FALSE)
  expect_named(result, c("union_exons", "isoform_union_mapping"))
  expect_true(all(c("gene_id", "union_exon_id", "chr", "start", "end",
                     "strand", "length") %in% names(result$union_exons)))
  expect_true(all(c("gene_id", "isoform_id", "union_exon_id") %in%
                    names(result$isoform_union_mapping)))
})
