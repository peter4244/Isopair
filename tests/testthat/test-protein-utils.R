# Tests for parseHmmscanOutput(), runHmmscan(), digestProtein(),
# queryPeptideAtlas()

# ============================================================================
# Helper: create mock domtblout file
# ============================================================================

.make_mock_domtblout <- function(tmp_path = tempfile(fileext = ".domtblout")) {
  # Format: 22 fixed columns + description
  # Columns: target_name accession tlen query_name query_accession qlen
  #   E-value score bias dom# ndom c-Evalue i-Evalue dom_score dom_bias
  #   hmm_from hmm_to ali_from ali_to env_from env_to acc description...
  lines <- c(
    "#                                                                            --- full sequence --- -------------- this domain ----------   hmm coord   ali coord   env coord",
    "# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias   from    to  from    to  from    to  acc description of target",
    "#--- -------------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------",
    "Tubulin              PF00091.32    431 ENST00000001         -    450   2.3e-55  189.4   0.1   1   2   3.5e-30   1.2e-28  170.2   0.0     5   195    10   210     8   212  0.95 Tubulin GTPase domain",
    "Kinase               PF00069.10    300 ENST00000001         -    450   1.1e-40  150.3   0.2   1   1   1.0e-20   5.5e-18  145.0   0.1     1   280    50   340    48   342  0.90 Protein kinase domain",
    "DUF9999              PF99999.1     100 ENST00000002         -    300   0.5      10.0    0.0   1   1   0.1       0.5       8.0    0.0     1    90    10    95     8    97  0.60 Domain of unknown function"
  )
  writeLines(lines, tmp_path)
  tmp_path
}


# ============================================================================
# parseHmmscanOutput()
# ============================================================================

test_that("parseHmmscanOutput parses domtblout correctly", {
  tmp <- .make_mock_domtblout()
  on.exit(unlink(tmp))

  result <- parseHmmscanOutput(tmp)

  # Should have 2 rows (DUF9999 has i-evalue 0.5 > default 1e-5)
  expect_equal(nrow(result), 2L)

  # Check column names
  expect_true(all(c("isoform_id", "domain_id", "domain_name", "domain_start",
                     "domain_end", "domain_score", "domain_evalue",
                     "description") %in% names(result)))

  # Check first row (Tubulin)
  tub <- result[result$domain_name == "Tubulin", ]
  expect_equal(nrow(tub), 1L)
  expect_equal(tub$isoform_id, "ENST00000001")
  expect_equal(tub$domain_id, "PF00091")  # version stripped
  expect_equal(tub$domain_start, 8L)      # env_from
  expect_equal(tub$domain_end, 212L)      # env_to
  expect_equal(tub$domain_score, 170.2)
  expect_true(grepl("Tubulin GTPase domain", tub$description))

  # Check second row (Kinase)
  kin <- result[result$domain_name == "Kinase", ]
  expect_equal(kin$domain_id, "PF00069")
  expect_equal(kin$domain_start, 48L)
  expect_equal(kin$domain_end, 342L)
})

test_that("parseHmmscanOutput respects evalue_threshold", {
  tmp <- .make_mock_domtblout()
  on.exit(unlink(tmp))

  # With a permissive threshold, all 3 should pass

  result_all <- parseHmmscanOutput(tmp, evalue_threshold = 1.0)
  expect_equal(nrow(result_all), 3L)

  # With strict threshold, only Tubulin passes
  result_strict <- parseHmmscanOutput(tmp, evalue_threshold = 1e-20)
  expect_equal(nrow(result_strict), 1L)
  expect_equal(result_strict$domain_name, "Tubulin")
})

test_that("parseHmmscanOutput returns empty tibble for comment-only file", {
  tmp <- tempfile(fileext = ".domtblout")
  writeLines(c("# comment line 1", "# comment line 2"), tmp)
  on.exit(unlink(tmp))

  result <- parseHmmscanOutput(tmp)
  expect_equal(nrow(result), 0L)
  expect_true(all(c("isoform_id", "domain_id", "domain_name") %in%
                    names(result)))
})

test_that("parseHmmscanOutput validates inputs", {
  expect_error(parseHmmscanOutput(42), "single character string")
  expect_error(parseHmmscanOutput("/no/such/file.txt"), "File not found")
  tmp <- .make_mock_domtblout()
  on.exit(unlink(tmp))
  expect_error(parseHmmscanOutput(tmp, evalue_threshold = -1),
               "positive number")
})


# ============================================================================
# runHmmscan()
# ============================================================================

test_that("runHmmscan errors gracefully when hmmscan not found", {
  # Mock Sys.which to return empty string
  mockr_available <- requireNamespace("testthat", quietly = TRUE)

  # Use withCallingHandlers/local to test without mockr
  # Create a temporary fake HMM database file
  fake_db <- tempfile(fileext = ".hmm")
  writeLines("fake", fake_db)
  on.exit(unlink(fake_db))

  seqs <- c(iso1 = "MKTLLVKPFLAG")

  # The real test: if hmmscan is not on PATH, the function should error
  # We test this by checking the error message pattern
  orig_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = "")
  on.exit(Sys.setenv(PATH = orig_path), add = TRUE)

  expect_error(
    runHmmscan(seqs, fake_db),
    "hmmscan.*not installed|not found"
  )
})

test_that("runHmmscan validates inputs", {
  expect_error(runHmmscan(character(0), "/fake/db"),
               "non-empty character vector")
  expect_error(runHmmscan(c(a = "MKT"), "/no/such/db.hmm"),
               "HMM database not found")

  fake_db <- tempfile(fileext = ".hmm")
  writeLines("fake", fake_db)
  on.exit(unlink(fake_db))

  expect_error(runHmmscan(c("MKT"), fake_db),  # unnamed
               "named character vector")
  expect_error(runHmmscan(c(a = "MKT"), fake_db, n_cores = 0L),
               "positive integer")
})


# ============================================================================
# digestProtein()
# ============================================================================

test_that("digestProtein performs trypsin digestion correctly", {
  # Simple sequence: MKTLLVK|PFLAGR|AATYR
  # Trypsin: cut after K/R, not before P
  # K at position 2: next char is T (not P) -> cut
  # K at position 7: next char is P -> NO cut (KP rule)
  # R at position 13: next char is A -> cut
  # R at position 18: end of sequence
  # So peptides: MK (1-2), TLLVKPFLAGR (3-13), AATYR (14-18)
  seqs <- c(iso1 = "MKTLLVKPFLAGRAATYR")
  result <- digestProtein(seqs, min_length = 1L, max_length = 100L)

  expect_equal(nrow(result), 3L)
  expect_equal(result$peptide, c("MK", "TLLVKPFLAGR", "AATYR"))
  expect_equal(result$start_aa, c(1L, 3L, 14L))
  expect_equal(result$end_aa, c(2L, 13L, 18L))
})

test_that("digestProtein handles stop codon markers", {
  seqs <- c(iso1 = "MKTLLVKPFLAGRAATYR*")
  result <- digestProtein(seqs, min_length = 1L, max_length = 100L)

  # Same result as without star
  expect_equal(nrow(result), 3L)
  expect_equal(result$peptide, c("MK", "TLLVKPFLAGR", "AATYR"))

  # Multiple stop codons
  seqs2 <- c(iso1 = "MKTLLVK***")
  result2 <- digestProtein(seqs2, min_length = 1L, max_length = 100L)
  expect_equal(nrow(result2), 2L)
  expect_equal(result2$peptide, c("MK", "TLLVK"))
})

test_that("digestProtein filters by length", {
  seqs <- c(iso1 = "MKTLLVKPFLAGRAATYR")
  # Default min_length = 7, max_length = 30
  result <- digestProtein(seqs)

  # MK (2aa) and AATYR (5aa) should be excluded, TLLVKPFLAGR (11aa) kept
  expect_equal(nrow(result), 1L)
  expect_equal(result$peptide, "TLLVKPFLAGR")
  expect_equal(result$start_aa, 3L)
  expect_equal(result$end_aa, 13L)
})

test_that("digestProtein handles multiple isoforms", {
  seqs <- c(iso1 = "MKTLLVKPFLAGR", iso2 = "AATYRLIVEGK")
  result <- digestProtein(seqs, min_length = 1L, max_length = 100L)

  iso1_rows <- result[result$isoform_id == "iso1", ]
  iso2_rows <- result[result$isoform_id == "iso2", ]

  expect_true(nrow(iso1_rows) > 0L)
  expect_true(nrow(iso2_rows) > 0L)
})

test_that("digestProtein handles sequence with no cleavage sites", {
  seqs <- c(iso1 = "MTTLLVVPFLAG")
  result <- digestProtein(seqs, min_length = 1L, max_length = 100L)

  expect_equal(nrow(result), 1L)
  expect_equal(result$peptide, "MTTLLVVPFLAG")
  expect_equal(result$start_aa, 1L)
  expect_equal(result$end_aa, 12L)
})

test_that("digestProtein validates inputs", {
  expect_error(digestProtein(character(0)), "non-empty")
  expect_error(digestProtein(c("MKT")), "named")
  expect_error(digestProtein(c(a = "MKT"), min_length = 0L),
               "positive integer")
  expect_error(digestProtein(c(a = "MKT"), min_length = 10L, max_length = 5L),
               "max_length")
})


# ============================================================================
# queryPeptideAtlas()
# ============================================================================

test_that("queryPeptideAtlas performs correct lookups", {
  atlas <- tempfile(fileext = ".txt")
  writeLines(c("MKTLLVK", "PFLAGR", "AATYR", "LIVEGK"), atlas)
  on.exit(unlink(atlas))

  result <- queryPeptideAtlas(c("MKTLLVK", "UNKNOWN", "AATYR"), atlas)

  expect_equal(nrow(result), 3L)
  expect_equal(result$peptide, c("MKTLLVK", "UNKNOWN", "AATYR"))
  expect_equal(result$found, c(TRUE, FALSE, TRUE))
})

test_that("queryPeptideAtlas handles empty atlas file", {
  atlas <- tempfile(fileext = ".txt")
  writeLines(character(0), atlas)
  on.exit(unlink(atlas))

  result <- queryPeptideAtlas(c("MKTLLVK"), atlas)
  expect_equal(nrow(result), 1L)
  expect_false(result$found)
})

test_that("queryPeptideAtlas validates inputs", {
  expect_error(queryPeptideAtlas(character(0), "/fake"),
               "non-empty character vector")
  expect_error(queryPeptideAtlas(c("MKT"), "/no/such/file.txt"),
               "File not found")
  expect_error(queryPeptideAtlas(c("MKT"), 42),
               "single character string")
})
