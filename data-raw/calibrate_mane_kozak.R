# Calibrate a Kozak-score threshold from MANE Select annotated CDS starts.
#
# Produces inst/extdata/kozak_mane_calibration.rds — a named list with:
#   threshold_q01, threshold_q05, threshold_q10, threshold_q25, threshold_q50
#   scores          : numeric vector of all per-MANE-isoform Kozak scores
#   n_mane_input    : number of MANE Select transcript IDs requested
#   n_mane_scored   : number successfully scored (after structure / CDS /
#                     sequence availability filters)
#   gencode_version : "GENCODE_v49"
#   computed_at     : timestamp
#
# Run once when MANE changes (annually) or when the default PWM changes.
# The 5th-percentile threshold is the canonical default used by
# enumerateOrfs() and selectPrimaryOrf() via defaultKozakThreshold().

suppressPackageStartupMessages({
  library(Isopair); library(Biostrings); library(dplyr); library(cli)
})

# ---- inputs -----------------------------------------------------------------
GTF    <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.sorted.gtf.gz"
GTF_RAW <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.primary_assembly.annotation.gtf.gz"
FASTA   <- "/Users/petecastaldi/claude_projects/nmd/reference_files/gencode.v49.transcripts.fa.gz"

# ---- 1. Enumerate MANE Select transcript IDs --------------------------------
cli::cli_h1("Calibrating Kozak threshold from MANE Select")

cli::cli_alert_info("Reading MANE Select IDs from the unsorted GENCODE GTF...")
# zcat | grep is fast for this ~93 MB file.
ids_raw <- system(sprintf(
  "zcat < %s | awk -F'\\t' '$3==\"transcript\" && /MANE_Select/' | grep -oE 'transcript_id \"[^\"]+\"' | sort -u",
  shQuote(GTF_RAW)
), intern = TRUE)
mane_ids <- sub('transcript_id \"([^\"]+)\"', '\\1', ids_raw)
cli::cli_alert_success("MANE Select transcripts: {length(mane_ids)}")

# ---- 2. Parse structures + CDS ---------------------------------------------
cli::cli_alert_info("Parsing structures for {length(mane_ids)} MANE transcripts (a few minutes)...")
t0 <- Sys.time()
structs <- parseIsoformStructures(GTF, isoform_ids = mane_ids, verbose = FALSE)
cds     <- extractCdsAnnotations(GTF, isoform_ids = mane_ids, verbose = FALSE)
cli::cli_alert_success("Parsed {nrow(structs)} structures, \\
                       {sum(cds$coding_status == 'coding')} coding CDS records \\
                       in {sprintf('%.0fs', as.numeric(Sys.time() - t0, units = 'secs'))}.")

# ---- 3. Load transcript sequences ------------------------------------------
cli::cli_alert_info("Loading transcript sequences from GENCODE FASTA...")
dss <- readDNAStringSet(FASTA)
keys <- sub("\\|.*$", "", sub("\\s.*$", "", names(dss)))
names(dss) <- keys
kept  <- intersect(mane_ids, keys)
seqs  <- toupper(as.character(dss[kept]))
names(seqs) <- kept
cli::cli_alert_success("Sequences available for {length(seqs)} / {length(mane_ids)} MANE transcripts.")

# ---- 4. Score Kozak at each annotated CDS start ----------------------------
cli::cli_alert_info("Scoring Kozak context at annotated CDS starts...")
t0 <- Sys.time()
emp <- empiricalKozakThreshold(structs, cds, seqs, quantile = 0.05)
cli::cli_alert_success("n_used = {emp$n_used}   \\
                       5%% threshold = {sprintf('%.4f', emp$threshold)}   \\
                       wall time = {sprintf('%.0fs', as.numeric(Sys.time() - t0, units = 'secs'))}")

# ---- 5. Report the distribution --------------------------------------------
cli::cli_h2("Score distribution (annotated MANE CDS starts)")
print(summary(emp$scores))
cat("\n")
for (q in c(0.01, 0.05, 0.10, 0.25, 0.50, 0.90)) {
  cat(sprintf("  quantile %4.2f : %.4f\n", q,
              stats::quantile(emp$scores, probs = q, na.rm = TRUE)))
}
cat("\n")
cat(sprintf("  mean  : %.4f\n", mean(emp$scores, na.rm = TRUE)))
cat(sprintf("  sd    : %.4f\n", sd(emp$scores, na.rm = TRUE)))
cat(sprintf("  fraction above 0 : %.3f\n",
            mean(emp$scores > 0, na.rm = TRUE)))

# ---- 6. Serialize to inst/extdata ------------------------------------------
out_path <- file.path(rprojroot::find_package_root_file(),
                     "inst", "extdata", "kozak_mane_calibration.rds")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

to_save <- list(
  threshold_q01 = as.numeric(stats::quantile(emp$scores, 0.01, na.rm = TRUE)),
  threshold_q05 = as.numeric(stats::quantile(emp$scores, 0.05, na.rm = TRUE)),
  threshold_q10 = as.numeric(stats::quantile(emp$scores, 0.10, na.rm = TRUE)),
  threshold_q25 = as.numeric(stats::quantile(emp$scores, 0.25, na.rm = TRUE)),
  threshold_q50 = as.numeric(stats::quantile(emp$scores, 0.50, na.rm = TRUE)),
  scores          = emp$scores,
  n_mane_input    = length(mane_ids),
  n_mane_scored   = emp$n_used,
  gencode_version = "GENCODE_v49",
  computed_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
saveRDS(to_save, out_path)
cli::cli_alert_success("Wrote {.path {out_path}}")
