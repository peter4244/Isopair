# Isopair 0.99.4

## New features

* New exported function `scoreKozakPWM()` computes continuous Kozak initiation
  strength scores using a log-odds position weight matrix (PWM) across 8
  positions (-6 to -1, +4, +5). Accepts custom weight matrices (e.g., Noderer
  et al. 2014 FACS-seq values). Vectorized for efficient scoring of many ATGs.
* `scan5UtrFeatures()` gains a `best_kozak_pwm_score` column — the highest
  PWM-based Kozak score among all uATGs (continuous, alongside the existing
  0-2 categorical `best_kozak_score`).
* `detectUorfs()` gains a `kozak_pwm_score` column per uORF (continuous,
  alongside the existing 0-2 categorical `kozak_score`).

# Isopair 0.99.3

## Bug fixes

* Fixed version bump for internal changes.

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
