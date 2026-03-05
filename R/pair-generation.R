#' Identify Dominant Isoforms Per Gene
#'
#' For each gene, identifies the isoform with the highest mean expression
#' across the specified samples. Optionally requires the dominant isoform to
#' exceed a proportion threshold (proportion of gene total).
#'
#' @param expression_matrix Numeric matrix with isoform IDs as rownames and
#'   sample IDs as column names.
#' @param gene_map A data frame with columns `isoform_id` and `gene_id`.
#' @param samples Character vector of column names to use from the expression
#'   matrix.
#' @param threshold Numeric; minimum proportion of gene total required for
#'   an isoform to be called dominant. Default 0 means highest-expressed
#'   regardless of proportion.
#' @return A tibble with columns: gene_id, dominant_isoform_id,
#'   mean_expression, proportion.
#' @examples
#' # Inline expression matrix for 4 isoforms across 2 genes
#' expr_mat <- matrix(
#'   c(100, 20, 80, 10, 110, 25, 75, 15),
#'   nrow = 4, dimnames = list(
#'     c("tx1", "tx2", "tx3", "tx4"), c("S1", "S2")
#'   )
#' )
#' gene_map <- data.frame(
#'   isoform_id = c("tx1", "tx2", "tx3", "tx4"),
#'   gene_id = c("g1", "g1", "g2", "g2")
#' )
#' dom <- identifyDominantIsoforms(expr_mat, gene_map, c("S1", "S2"))
#' dom
#' @export
#' @importFrom dplyr filter select mutate group_by summarise arrange
#'   row_number ungroup left_join bind_rows n
#' @importFrom tibble tibble
#' @importFrom rlang .data
identifyDominantIsoforms <- function(expression_matrix, gene_map,
                                     samples, threshold = 0) {

  if (!is.matrix(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix with rownames.")
  }
  if (is.null(rownames(expression_matrix))) {
    stop("expression_matrix must have rownames (isoform IDs).")
  }
  if (!all(c("isoform_id", "gene_id") %in% names(gene_map))) {
    stop("gene_map must have columns 'isoform_id' and 'gene_id'.")
  }
  if (!all(samples %in% colnames(expression_matrix))) {
    missing <- setdiff(samples, colnames(expression_matrix))
    stop(sprintf("Samples not found in expression_matrix: %s",
                 paste(missing, collapse = ", ")))
  }

  # Compute mean expression per isoform across selected samples
  sub_mat <- expression_matrix[, samples, drop = FALSE]
  mean_expr <- rowMeans(sub_mat)

  # Build isoform-level table
  iso_df <- tibble::tibble(
    isoform_id = rownames(expression_matrix),
    mean_expression = as.numeric(mean_expr)
  )

  # Join with gene map (only keep isoforms in both)
  iso_df <- dplyr::left_join(iso_df, gene_map, by = "isoform_id")
  iso_df <- dplyr::filter(iso_df, !is.na(.data$gene_id))

  # Compute gene totals and proportions
  gene_totals <- iso_df |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(gene_total = sum(.data$mean_expression), .groups = "drop")

  iso_df <- dplyr::left_join(iso_df, gene_totals, by = "gene_id")
  iso_df <- dplyr::mutate(iso_df,
    proportion = ifelse(.data$gene_total > 0,
                        .data$mean_expression / .data$gene_total,
                        0))

  # For each gene, pick top isoform by mean expression
  dominant <- iso_df |>
    dplyr::filter(.data$mean_expression > 0) |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::arrange(dplyr::desc(.data$mean_expression)) |>
    dplyr::filter(dplyr::row_number() == 1L) |>
    dplyr::ungroup()

  # Apply threshold filter
  if (threshold > 0) {
    dominant <- dplyr::filter(dominant, .data$proportion >= threshold)
  }

  dplyr::select(dominant,
    "gene_id", dominant_isoform_id = "isoform_id",
    "mean_expression", "proportion")
}


#' Generate Isoform Pairs from Expression Data
#'
#' Creates pairs of isoforms for structural comparison based on expression
#' ranking within each gene.
#'
#' @param expression_matrix Numeric matrix with isoform IDs as rownames.
#' @param gene_map A data frame with columns `isoform_id` and `gene_id`.
#' @param samples Character vector of column names to use.
#' @param method Character; either `"top_two"` (pair rank 1 vs rank 2 per
#'   gene) or `"dominant_vs_rest"` (pair dominant vs all others).
#' @param dominant_threshold Numeric; proportion threshold for dominant
#'   isoform selection in `"dominant_vs_rest"` mode. Default 0.5.
#' @param min_expression Numeric; minimum mean expression for a comparator
#'   isoform. Default 0 (no filter).
#' @return A tibble with columns: gene_id, reference_isoform_id,
#'   comparator_isoform_id, source, direction, logFC, adj_p_val.
#' @examples
#' expr_file <- system.file("extdata", "example_expression.csv",
#'   package = "Isopair")
#' expr_df <- read.csv(expr_file)
#' expr_mat <- as.matrix(expr_df[, -1])
#' rownames(expr_mat) <- expr_df$isoform_id
#' data(example_structures)
#' gene_map <- unique(example_structures[, c("isoform_id", "gene_id")])
#' pairs <- generatePairsExpression(expr_mat, gene_map, colnames(expr_mat))
#' head(pairs)
#' @export
generatePairsExpression <- function(expression_matrix, gene_map, samples,
                                    method = "top_two",
                                    dominant_threshold = 0.5,
                                    min_expression = 0) {

  method <- match.arg(method, c("top_two", "dominant_vs_rest"))

  # Compute mean expression
  sub_mat <- expression_matrix[, samples, drop = FALSE]
  mean_expr <- rowMeans(sub_mat)

  iso_df <- tibble::tibble(
    isoform_id = rownames(expression_matrix),
    mean_expression = as.numeric(mean_expr)
  )
  iso_df <- dplyr::left_join(iso_df, gene_map, by = "isoform_id")
  iso_df <- dplyr::filter(iso_df, !is.na(.data$gene_id))

  if (method == "top_two") {
    # Rank by expression within each gene, pair rank 1 vs rank 2
    ranked <- iso_df |>
      dplyr::filter(.data$mean_expression > 0) |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::arrange(dplyr::desc(.data$mean_expression)) |>
      dplyr::mutate(rank = dplyr::row_number()) |>
      dplyr::ungroup()

    top1 <- dplyr::filter(ranked, .data$rank == 1L)
    top2 <- dplyr::filter(ranked, .data$rank == 2L)

    pairs <- dplyr::left_join(
      dplyr::select(top1, "gene_id",
                    reference_isoform_id = "isoform_id"),
      dplyr::select(top2, "gene_id",
                    comparator_isoform_id = "isoform_id"),
      by = "gene_id"
    )
    pairs <- dplyr::filter(pairs, !is.na(.data$comparator_isoform_id))

  } else {
    # dominant_vs_rest
    dominant <- identifyDominantIsoforms(
      expression_matrix, gene_map, samples,
      threshold = dominant_threshold
    )

    # All isoforms with expression
    all_iso <- iso_df |>
      dplyr::filter(.data$mean_expression > min_expression)

    # For each gene with a dominant, pair dominant vs all others
    pairs_list <- list()
    for (i in seq_len(nrow(dominant))) {
      g <- dominant$gene_id[i]
      dom_id <- dominant$dominant_isoform_id[i]
      others <- all_iso |>
        dplyr::filter(.data$gene_id == g, .data$isoform_id != dom_id)
      if (nrow(others) > 0L) {
        pairs_list[[i]] <- tibble::tibble(
          gene_id = g,
          reference_isoform_id = dom_id,
          comparator_isoform_id = others$isoform_id
        )
      }
    }
    pairs <- dplyr::bind_rows(pairs_list)
  }

  if (nrow(pairs) == 0L) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0), source = character(0),
      direction = character(0), logFC = numeric(0),
      adj_p_val = numeric(0)
    ))
  }

  dplyr::mutate(pairs,
    source = "expression",
    direction = NA_character_,
    logFC = NA_real_,
    adj_p_val = NA_real_
  )
}


#' Generate Isoform Pairs from Differential Expression Results
#'
#' Pairs each significantly differentially expressed isoform against the
#' dominant isoform of its gene.
#'
#' @param de_results A data frame of differential expression results.
#' @param isoform_col Character; column name for isoform IDs.
#' @param gene_col Character; column name for gene IDs.
#' @param logfc_col Character; column name for log fold change.
#' @param pval_col Character; column name for adjusted p-value.
#' @param dominant_isoforms A tibble from [identifyDominantIsoforms()]
#'   (required).
#' @param pval_threshold Numeric; adjusted p-value cutoff. Default 0.05.
#' @param logfc_threshold Numeric; absolute logFC cutoff. Default 0.
#' @return A tibble with columns: gene_id, reference_isoform_id,
#'   comparator_isoform_id, source, direction, logFC, adj_p_val.
#' @examples
#' de_file <- system.file("extdata", "example_de.csv", package = "Isopair")
#' de <- read.csv(de_file)
#' expr_file <- system.file("extdata", "example_expression.csv",
#'   package = "Isopair")
#' expr_df <- read.csv(expr_file)
#' expr_mat <- as.matrix(expr_df[, -1])
#' rownames(expr_mat) <- expr_df$isoform_id
#' data(example_structures)
#' gene_map <- unique(example_structures[, c("isoform_id", "gene_id")])
#' dom <- identifyDominantIsoforms(expr_mat, gene_map, colnames(expr_mat))
#' pairs <- generatePairsDE(de, "isoform_id", "gene_id", "logFC",
#'   "adj_p_val", dom)
#' head(pairs)
#' @export
generatePairsDE <- function(de_results, isoform_col, gene_col,
                            logfc_col, pval_col,
                            dominant_isoforms,
                            pval_threshold = 0.05,
                            logfc_threshold = 0) {

  # Validate columns exist
  required_cols <- c(isoform_col, gene_col, logfc_col, pval_col)
  missing_cols <- setdiff(required_cols, names(de_results))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Columns not found in de_results: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Extract columns by name
  de_df <- tibble::tibble(
    isoform_id = de_results[[isoform_col]],
    gene_id = de_results[[gene_col]],
    logFC = as.numeric(de_results[[logfc_col]]),
    adj_p_val = as.numeric(de_results[[pval_col]])
  )

  # Filter for significant
  sig <- dplyr::filter(de_df,
    .data$adj_p_val < pval_threshold,
    abs(.data$logFC) > logfc_threshold
  )

  if (nrow(sig) == 0L) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0), source = character(0),
      direction = character(0), logFC = numeric(0),
      adj_p_val = numeric(0)
    ))
  }

  # Join with dominant isoforms
  sig <- dplyr::left_join(
    sig,
    dplyr::select(dominant_isoforms, "gene_id",
                  reference_isoform_id = "dominant_isoform_id"),
    by = "gene_id"
  )

  # Remove genes without dominant, and self-pairs
  sig <- dplyr::filter(sig,
    !is.na(.data$reference_isoform_id),
    .data$isoform_id != .data$reference_isoform_id
  )

  if (nrow(sig) == 0L) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0), source = character(0),
      direction = character(0), logFC = numeric(0),
      adj_p_val = numeric(0)
    ))
  }

  tibble::tibble(
    gene_id = sig$gene_id,
    reference_isoform_id = sig$reference_isoform_id,
    comparator_isoform_id = sig$isoform_id,
    source = "DE",
    direction = ifelse(sig$logFC > 0, "up", "down"),
    logFC = sig$logFC,
    adj_p_val = sig$adj_p_val
  )
}


#' Generate Isoform Pairs from Differential Usage Results
#'
#' Pairs isoforms based on differential usage (e.g., from DTU analysis).
#' Two methods are available: pairing most extreme positive/negative dPSI
#' isoforms per gene, or pairing the dominant isoform against the most
#' significant DU isoform.
#'
#' @param du_results A data frame of differential usage results.
#' @param isoform_col Character; column name for isoform IDs.
#' @param gene_col Character; column name for gene IDs.
#' @param dpsi_col Character; column name for delta PSI (or similar metric).
#' @param pval_col Character; column name for adjusted p-value.
#' @param dominant_isoforms A tibble from [identifyDominantIsoforms()].
#'   Required for `method = "dominant_vs_sig"`.
#' @param pval_threshold Numeric; adjusted p-value cutoff. Default 0.05.
#' @param method Character; either `"extreme"` (pair most positive vs most
#'   negative dPSI per gene) or `"dominant_vs_sig"` (pair dominant vs most
#'   significant DU isoform per gene).
#' @return A tibble with columns: gene_id, reference_isoform_id,
#'   comparator_isoform_id, source, direction, logFC, adj_p_val.
#' @examples
#' du_file <- system.file("extdata", "example_du.csv", package = "Isopair")
#' du <- read.csv(du_file)
#' pairs <- generatePairsDU(du, "isoform_id", "gene_id", "dPSI",
#'   "adj_p_val", method = "extreme")
#' head(pairs)
#' @export
generatePairsDU <- function(du_results, isoform_col, gene_col,
                            dpsi_col, pval_col,
                            dominant_isoforms = NULL,
                            pval_threshold = 0.05,
                            method = "dominant_vs_sig") {

  method <- match.arg(method, c("extreme", "dominant_vs_sig"))

  # Validate columns exist
  required_cols <- c(isoform_col, gene_col, dpsi_col, pval_col)
  missing_cols <- setdiff(required_cols, names(du_results))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Columns not found in du_results: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  if (method == "dominant_vs_sig" && is.null(dominant_isoforms)) {
    stop("dominant_isoforms is required for method = 'dominant_vs_sig'.")
  }

  # Extract columns by name
  du_df <- tibble::tibble(
    isoform_id = du_results[[isoform_col]],
    gene_id = du_results[[gene_col]],
    dpsi = as.numeric(du_results[[dpsi_col]]),
    adj_p_val = as.numeric(du_results[[pval_col]])
  )

  # Filter for significant
  sig <- dplyr::filter(du_df, .data$adj_p_val < pval_threshold)

  if (nrow(sig) == 0L) {
    return(tibble::tibble(
      gene_id = character(0), reference_isoform_id = character(0),
      comparator_isoform_id = character(0), source = character(0),
      direction = character(0), logFC = numeric(0),
      adj_p_val = numeric(0)
    ))
  }

  empty_result <- tibble::tibble(
    gene_id = character(0), reference_isoform_id = character(0),
    comparator_isoform_id = character(0), source = character(0),
    direction = character(0), logFC = numeric(0),
    adj_p_val = numeric(0)
  )

  if (method == "extreme") {
    # Per gene: pair most-positive-dPSI vs most-negative-dPSI
    pos <- sig |>
      dplyr::filter(.data$dpsi > 0) |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::arrange(dplyr::desc(.data$dpsi)) |>
      dplyr::filter(dplyr::row_number() == 1L) |>
      dplyr::ungroup()

    neg <- sig |>
      dplyr::filter(.data$dpsi < 0) |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::arrange(.data$dpsi) |>
      dplyr::filter(dplyr::row_number() == 1L) |>
      dplyr::ungroup()

    # Only genes with both positive and negative
    shared_genes <- intersect(pos$gene_id, neg$gene_id)
    if (length(shared_genes) == 0L) return(empty_result)

    pos <- dplyr::filter(pos, .data$gene_id %in% shared_genes)
    neg <- dplyr::filter(neg, .data$gene_id %in% shared_genes)

    pairs <- dplyr::left_join(
      dplyr::select(pos, "gene_id",
                    reference_isoform_id = "isoform_id",
                    ref_dpsi = "dpsi", ref_pval = "adj_p_val"),
      dplyr::select(neg, "gene_id",
                    comparator_isoform_id = "isoform_id",
                    comp_dpsi = "dpsi", comp_pval = "adj_p_val"),
      by = "gene_id"
    )

    result <- tibble::tibble(
      gene_id = pairs$gene_id,
      reference_isoform_id = pairs$reference_isoform_id,
      comparator_isoform_id = pairs$comparator_isoform_id,
      source = "DU",
      direction = "extreme",
      logFC = pairs$ref_dpsi - pairs$comp_dpsi,
      adj_p_val = pmin(pairs$ref_pval, pairs$comp_pval)
    )

  } else {
    # dominant_vs_sig: pair dominant vs most significant DU isoform per gene
    sig <- dplyr::left_join(
      sig,
      dplyr::select(dominant_isoforms, "gene_id",
                    dominant_isoform_id = "dominant_isoform_id"),
      by = "gene_id"
    )

    # Remove genes without dominant, and self-pairs
    sig <- dplyr::filter(sig,
      !is.na(.data$dominant_isoform_id),
      .data$isoform_id != .data$dominant_isoform_id
    )

    if (nrow(sig) == 0L) return(empty_result)

    # Per gene, pick most significant (lowest p-value)
    top_sig <- sig |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::arrange(.data$adj_p_val) |>
      dplyr::filter(dplyr::row_number() == 1L) |>
      dplyr::ungroup()

    result <- tibble::tibble(
      gene_id = top_sig$gene_id,
      reference_isoform_id = top_sig$dominant_isoform_id,
      comparator_isoform_id = top_sig$isoform_id,
      source = "DU",
      direction = ifelse(top_sig$dpsi > 0, "up", "down"),
      logFC = top_sig$dpsi,
      adj_p_val = top_sig$adj_p_val
    )
  }

  result
}
