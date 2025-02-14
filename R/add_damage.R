#' simulate_damage
#'
#' @param count_matrix Existing scRNA-seq data in the form of a matrix
#' @param organism String specifying what organism the input data originated from
#' @param damage_proportion Proportion of damage to be introduced to the sample
#' @param beta_proportion
#' @param lambda
#' @param n_samples
#'
#' @return
#' @export
#'
#' @examples
simulate_damage <- function(count_matrix,
                            organism = "Hsap",
                            damage_proportion = 0.5,
                            beta_proportion = 0.5,
                            lambda = 10,
                            n_samples = 1000
) {

  # Retrieve genes corresponding to the organism of interest
  if (organism == "Hsap"){
    data("mito_genes", package = "damagesimulator")
    ribo_pattern <- "^(RPS|RPL|MRPS|MRPL)"
  }

  if (organism == "Mmus"){
    data("mouse_mito_genes", package = "damagesimulator")
    mito_genes <- mouse_mito_genes
    ribo_pattern <- "^(Rps|Rpl|Mrps|Mrpl)"
  }

  # Randomly select proportion of cells to damage
  total_cells <- ncol(count_matrix)
  damaged_cell_number <- round(total_cells * damage_proportion)
  damaged_cell_selections <- sample(seq_len(total_cells), size = damaged_cell_number, replace = FALSE)

  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)


  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)


  # Isolate gene set indices
  # LOAD FROM PACKAGE
  mito_idx <- grep(mito_genes, rownames(count_matrix), ignore.case = FALSE)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)
  other_idx <- setdiff(seq_len(nrow(count_matrix)), union(mito_idx, ribo_idx))


  # Initialize for storing modified counts
  new_matrix <- count_matrix

  # Define the number of Monte Carlo samples
  n_samples <- n_samples # tested 1000 10000 100000 & found no helpful shift in shape but large increase in time

  # Loop over the damaged cells and apply the reduction to non-mito genes.
  for (i in seq_along(damaged_cell_selections)) {
    #i = 100  # JUST FOR TESTING

    # Index in the count matrix for the cell
    cell <- damaged_cell_selections[i]

    # Target mito proportion for this cell
    target <- damage_levels[i]

    # Compute initial gene sums in the cell
    M_i <- sum(count_matrix[mito_idx, cell])  # Mitochondrial
    R_i <- sum(count_matrix[ribo_idx, cell])  # Ribosomal
    O_i <- sum(count_matrix[other_idx, cell]) # Other (non-mito & non-ribo)
    T_i <- R_i + O_i  # Total non-mito counts (everything that must be reduced)

    # Monte Carlo sampling to estimate r
    r_samples <- runif(n_samples, 0.01, 0.7)  # Generate random r values in [0,1]


    # Compute A, B, and the absolute error for each sample
    A_values <- (r_samples * R_i) / (M_i + r_samples * T_i)
    B_values <- M_i / (M_i + r_samples * T_i)
    exp_decay_values <- exp(-lambda * A_values)
    errors <- abs(B_values - exp_decay_values)

    # Select the best r (minimizing the error)
    best_r <- r_samples[which.min(errors)]

    # Apply reduction to non-mito genes (ribo and other genes) in this cell.
    nonMitoGenes <- c(ribo_idx, other_idx)
    new_matrix[nonMitoGenes, cell] <- round(best_r * count_matrix[nonMitoGenes, cell])

  }

  # QC statistics for all cells
  qc_summary <- data.frame(
    Cell = colnames(count_matrix),
    Damaged_Status = damage_label$status[match(colnames(count_matrix), damage_label$barcode)],
    Original_Features = colSums(count_matrix != 0),
    New_Features = colSums(new_matrix != 0),
    Original_MitoProp = colSums(count_matrix[mito_idx, , drop = FALSE]) / colSums(count_matrix),
    New_MitoProp = colSums(new_matrix[mito_idx, , drop = FALSE]) / colSums(new_matrix),
    Original_RiboProp = colSums(count_matrix[ribo_idx, , drop = FALSE]) / colSums(count_matrix),
    New_RiboProp = colSums(new_matrix[ribo_idx, , drop = FALSE]) / colSums(new_matrix)
  )

  # Return a list containing the new count matrix and the QC summary.
  return(list(new_matrix = new_matrix, qc_summary = qc_summary))

}
