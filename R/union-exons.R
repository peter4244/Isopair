#' Build Atomic Union Exon Models
#'
#' Creates atomic (non-overlapping) union exon segments by splitting at all
#' unique exon boundaries within each gene. Each atomic segment represents the
#' smallest unit that is either fully included or fully excluded in any given
#' isoform.
#'
#' @param structures A tibble from [parseIsoformStructures()] with nested
#'   exon_starts and exon_ends list columns.
#' @param verbose Logical; if TRUE, print progress messages.
#' @return A named list with two elements:
#'   \describe{
#'     \item{union_exons}{Tibble with columns: gene_id, union_exon_id,
#'       union_exon_number, chr, start, end, strand, length.}
#'     \item{isoform_union_mapping}{Tibble with columns: gene_id, isoform_id,
#'       union_exon_id, exon_number_in_isoform, isoform_exon_start,
#'       isoform_exon_end.}
#'   }
#' @examples
#' gtf_path <- system.file("extdata", "example_small.gtf", package = "Isopair")
#' structures <- parseIsoformStructures(gtf_path, verbose = FALSE)
#' ue <- buildUnionExons(structures, verbose = FALSE)
#' head(ue$union_exons)
#' @export
#' @importFrom dplyr n_distinct
#' @importFrom tibble tibble
#' @importFrom rlang .data
buildUnionExons <- function(structures, verbose = TRUE) {
  if (verbose) message("Building atomic union exon models...")

  gene_list <- split(structures, structures$gene_id)
  genes <- names(gene_list)

  if (verbose) message(sprintf("  %d genes to process", length(genes)))

  all_ue <- list()
  all_map <- list()

  for (gene in genes) {
    gene_isoforms <- gene_list[[gene]]
    if (nrow(gene_isoforms) == 0L) next

    chr_val <- gene_isoforms$chr[1L]
    strand_val <- gene_isoforms$strand[1L]

    # Expand list columns to vectors
    n_exons_per_iso <- lengths(gene_isoforms$exon_starts)
    all_exon_starts <- unlist(gene_isoforms$exon_starts)
    all_exon_ends <- unlist(gene_isoforms$exon_ends)
    all_iso_ids <- rep(gene_isoforms$isoform_id, n_exons_per_iso)
    all_exon_nums <- unlist(lapply(n_exons_per_iso, seq_len))

    if (length(all_exon_starts) == 0L) next

    # Unique sorted boundaries
    all_boundaries <- sort(unique(c(all_exon_starts, all_exon_ends)))
    if (length(all_boundaries) < 2L) next

    # Classify boundaries
    starts_set <- unique(all_exon_starts)
    ends_set <- unique(all_exon_ends)
    is_start <- all_boundaries %in% starts_set
    is_end <- all_boundaries %in% ends_set

    # Build segments
    n_bounds <- length(all_boundaries)
    n_both <- sum(is_start & is_end)
    max_segs <- n_bounds - 1L + n_both
    seg_s <- integer(max_segs)
    seg_e <- integer(max_segs)
    n_seg <- 0L

    current_start <- all_boundaries[1L]

    if (n_bounds > 2L) {
      for (k in 2L:(n_bounds - 1L)) {
        b <- all_boundaries[k]
        if (is_start[k] && is_end[k]) {
          n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b - 1L
          n_seg <- n_seg + 1L; seg_s[n_seg] <- b; seg_e[n_seg] <- b
          current_start <- b + 1L
        } else if (is_start[k]) {
          n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b - 1L
          current_start <- b
        } else {
          n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- b
          current_start <- b + 1L
        }
      }
    }

    # Final segment
    n_seg <- n_seg + 1L; seg_s[n_seg] <- current_start; seg_e[n_seg] <- all_boundaries[n_bounds]

    # Trim and filter phantoms
    seg_s <- seg_s[seq_len(n_seg)]
    seg_e <- seg_e[seq_len(n_seg)]
    valid <- seg_s <= seg_e
    seg_s <- seg_s[valid]
    seg_e <- seg_e[valid]

    if (length(seg_s) == 0L) next

    # Coverage check: each segment must be covered by at least one exon
    n_segs <- length(seg_s)
    covered <- logical(n_segs)
    for (j in seq_len(n_segs)) {
      covered[j] <- any(all_exon_starts <= seg_s[j] & all_exon_ends >= seg_e[j])
    }
    seg_s <- seg_s[covered]
    seg_e <- seg_e[covered]

    if (length(seg_s) == 0L) next

    # Create union exon records
    n_ue <- length(seg_s)
    ue_ids <- paste0(gene, "_UE", seq_len(n_ue))

    all_ue[[gene]] <- tibble::tibble(
      gene_id = gene,
      union_exon_id = ue_ids,
      union_exon_number = seq_len(n_ue),
      chr = chr_val,
      start = seg_s,
      end = seg_e,
      strand = strand_val,
      length = seg_e - seg_s + 1L
    )

    # Isoform-to-UE mapping via containment
    gene_maps <- list()
    unique_isos <- unique(all_iso_ids)

    for (iso_id in unique_isos) {
      mask <- all_iso_ids == iso_id
      ex_starts <- all_exon_starts[mask]
      ex_ends <- all_exon_ends[mask]
      ex_nums <- all_exon_nums[mask]

      for (e_idx in seq_along(ex_starts)) {
        contained_mask <- seg_s >= ex_starts[e_idx] & seg_e <= ex_ends[e_idx]
        if (any(contained_mask)) {
          which_contained <- which(contained_mask)
          gene_maps[[length(gene_maps) + 1L]] <- tibble::tibble(
            gene_id = gene,
            isoform_id = iso_id,
            union_exon_id = ue_ids[which_contained],
            exon_number_in_isoform = ex_nums[e_idx],
            isoform_exon_start = ex_starts[e_idx],
            isoform_exon_end = ex_ends[e_idx]
          )
        }
      }
    }

    if (length(gene_maps) > 0L) {
      all_map[[gene]] <- dplyr::bind_rows(gene_maps)
    }
  }

  union_exons <- dplyr::bind_rows(all_ue)
  isoform_union_mapping <- dplyr::bind_rows(all_map)

  if (verbose) {
    message(sprintf("  Created %d atomic union exons across %d genes",
                    nrow(union_exons), dplyr::n_distinct(union_exons$gene_id)))
    message(sprintf("  Created %d isoform-UE mappings",
                    nrow(isoform_union_mapping)))
  }

  list(
    union_exons = union_exons,
    isoform_union_mapping = isoform_union_mapping
  )
}
