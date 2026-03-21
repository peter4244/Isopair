# Isopair 0.99.2

## Bug fixes

* `annotateRegionTypes()` now correctly assigns strand-aware region labels.
  Previously, 5'UTR/3'UTR and contains_orf_start/contains_orf_stop were
  based on genomic coordinate position without considering strand, causing
  all directional labels to be swapped for minus-strand genes (~50% of genes).
* `extractCdsAnnotations()` now returns a `strand` column ("+", "-", or NA
  for unknown isoforms). This column is required by `annotateRegionTypes()`.

# Isopair 0.99.1

* New `deduplicateStructures()` for merging and deduplicating isoform structures
  from multiple annotation sources (e.g., GENCODE + PacBio).
* New `bootstrapRegionalComparison()` for permutation-based testing of
  event-type x region-type proportion differences between two profile sets.
* `buildProfiles()` gains `checkpoint_dir` and `checkpoint_interval` parameters
  for intermediate saves and resumable processing of large pair sets.

# Isopair 0.99.0

* Initial Bioconductor submission.
* 12-type splicing event detection engine with reconstruction verification.
* Three pair generation methods: expression, differential expression, differential usage.
* Single-set characterization: co-occurrence, spatial, regional, PTC analysis.
* Cross-set comparison framework with self-documenting results.
* Random-effects meta-analysis via metafor.
