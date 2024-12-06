#' Perturb Single Cell
#'
#' This function perturbs a cell from a single-cell RNA-seq count matrix by adjusting the proportion of mitochondrial genes.
#'
#' @param count_matrix A matrix of gene expression counts with genes as rows and cells as columns.
#' @param mito_genes A character vector of mitochondrial gene names.
#' @param zero_out_constant A numeric constant to adjust the zero-out proportion (default is 0.0).
#' @return A perturbed count matrix with adjusted mitochondrial gene proportions.
#' @examples
#' # Example usage:
#' # count_matrix <- matrix(data = rpois(1000, lambda = 10), nrow = 100, ncol = 10)
#' # mito_genes <- c("MT-CO1", "MT-CO2", "MT-CO3")
#' # perturbed_matrix <- perturb_single_cell(count_matrix, mito_genes)
#' @export
perturb_single_cell <- function(count_matrix, mito_genes = NULL, zero_out_constant = 0.0) {

  # Load default mito_genes if not provided
  if (is.null(mito_genes)) {
    data("mito_genes", package = "damagesimulator")
  }

  # Convert count matrix to matrix for manipulation
  count_matrix <- as.matrix(count_matrix)

  # Identify mitochondrial and non-mitochondrial gene indices
  mito_idx <- rownames(count_matrix) %in% mito_genes
  non_mito_idx <- !mito_idx

  # Precompute mitochondrial and non-mitochondrial gene sums for each cell
  mito_counts <- colSums(count_matrix[mito_idx, , drop=FALSE])
  non_mito_counts <- colSums(count_matrix[non_mito_idx, , drop=FALSE])
  total_counts <- colSums(count_matrix)

  # Randomly assign target mitochondrial percentages for each cell
  target_mito_pct <- runif(ncol(count_matrix), min = 20, max = 98)

  for (i in seq_len(ncol(count_matrix))) {
    # Step 1: Calculate target non-mitochondrial gene count to zero
    current_non_mito_genes <- sum(count_matrix[non_mito_idx, i] > 0)

    # Define a scaling curve to control zero-out proportion
    max_zero_out_proportion <- 0.32  # testing to see if higher helps smooth it out # 32
    zero_out_proportion <- max_zero_out_proportion * (1 - exp(-0.02 * target_mito_pct[i]))  # Subtle scaling with cap

    # Add a small constant if specified
    zero_out_proportion <- min(zero_out_proportion + zero_out_constant, max_zero_out_proportion)

    # Calculate the number of genes to zero out, ensuring it doesn't exceed available genes
    genes_to_zero <- min(floor(current_non_mito_genes * zero_out_proportion), current_non_mito_genes)

    # Randomly zero out the calculated number of non-mito genes
    if (genes_to_zero > 0) {
      non_mito_gene_indices <- which(non_mito_idx & count_matrix[, i] > 0)
      genes_to_zero_indices <- sample(non_mito_gene_indices, size = genes_to_zero)
      count_matrix[genes_to_zero_indices, i] <- 0
    }

    # Step 2: Calculate new mito and non-mito counts after zeroing
    mito_sum <- sum(count_matrix[mito_idx, i])
    non_mito_sum <- sum(count_matrix[non_mito_idx, i])
    total_counts[i] <- mito_sum + non_mito_sum  # Update total counts

    # Step 3: Apply scaling to remaining non-mito genes to reach target mito percentage
    target_mito_counts <- total_counts[i] * (target_mito_pct[i] / 100)
    reduction_factor <- (total_counts[i] - target_mito_counts) / max(non_mito_sum, 1)  # No restriction on reduction factor

    # Apply the calculated reduction factor directly
    count_matrix[non_mito_idx, i] <- count_matrix[non_mito_idx, i] * reduction_factor

    # Optional: display iteration details
    cat("Processed cell", i, "- Target mito %", round(target_mito_pct[i], 2),
        "- Remaining non-mito genes:", sum(count_matrix[non_mito_idx, i] > 0), "\n")
  }

  return(count_matrix)
}
