
# Load libraries ----
library(Matrix)       # For sparse matrices
library(ggplot2)      # For plotting
library(scales)       # For number formatting in ggplot2
library(patchwork)    # For combining multiple ggplot2 plots
library(scRNAseq)

# Damage simulation function ----
simulate_damage <- function(count_matrix, 
                                 damage_proportion = 0.2, 
                                 beta_proportion = 0.5, 
                                 lambda = 5,
                                 damage_distribution_shapes = c(a1 = 3, b1 = 20,  
                                                                a2 = 25, b2 = 5))  # Low ~ 0.05–0.15 High ~ 0.8–0.9) {
{  
  # Identify mutally exclusive gene sets
  mito_idx <- grep("^MT-", rownames(count_matrix), ignore.case = FALSE)
  ribo_idx <- grep("^(RPS|RPL)", rownames(count_matrix), ignore.case = FALSE)
  other_idx <- setdiff(seq_len(nrow(count_matrix)), union(mito_idx, ribo_idx))
  
  # Sample a number of cells to be damaged
  total_cells <- ncol(count_matrix)
  damaged_cell_number <- round(total_cells * damage_proportion)
  damaged_cell_selections <- sample(seq_len(total_cells), size = damaged_cell_number, replace = FALSE)
  
  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  udamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  udamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[udamaged_cell_number_cells], status = rep("control", length(udamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, udamaged_cell_number)
  
  # Prepare a new count matrix to store the modified counts.
  new_matrix <- count_matrix
  
  # Initialise storage for damage levels
  reduction_levels <- numeric(damaged_cell_number)  # Create a vector for damage levels
  
  #Compute damage level for each damaged cell
  for (i in seq_along(damaged_cell_selections)) {
    cell <- damaged_cell_selections[i]
    #target <- damage_levels[i]
    
    M_i <- sum(count_matrix[mito_idx, cell])  # mitochondrial counts
    R_i <- sum(count_matrix[ribo_idx, cell])  # ribosomal counts
    O_i <- sum(count_matrix[other_idx, cell]) # other counts 
    T_i <- R_i + O_i
    
    B <- R_i / (M_i + T_i)  # Ribosomal proportion
    exp_decay <- exp(-lambda * B)
    
    r_i <- (1 - exp_decay) / exp_decay  # Exponential decay model
    
    reduction_levels[i] <- r_i  # Store the reduction level for this cell
  }
  
  # Prepare matrix form for reduction
  non_mito_genes <- c(ribo_idx, other_idx)
  
  
  # Create a matrix with the damage levels repeated for each gene and each cell
  perturb_factors_matrix <- matrix(rep(reduction_levels, each = length(non_mito_genes)),
                                   nrow = length(non_mito_genes),
                                   ncol = length(damaged_cell_selections),
                                   byrow = FALSE)
  
  # Apply reduction to the new matrix
  new_matrix[non_mito_genes, damaged_cell_selections] <- round(perturb_factors_matrix * count_matrix[non_mito_genes, damaged_cell_selections])
  
  # QC statistics for all cells
  qcSummary <- data.frame(
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
  return(list(new_matrix = new_matrix, qcSummary = qcSummary))
  
}


# Testing ----

# counts <- readRDS("Projects/ReviewArticle/simulated_data/count_matrix.rds")
# Or available in R 
sce <- fetchDataset(name = "he-organs-2020", version = "2023-12-21", path = "blood")
dim(sce) # 14552 features 1407 cells 
counts <- sce@assays@data$counts
rownames(counts) <- rownames(sce)
colnames(counts) <- colnames(sce)
counts <- as(counts, "sparseMatrix") 
count_matrix <- counts

# Run the function 
results <- simulate_damage(count_matrix = counts, damage_prop = 0.9,  beta_proportion = 1)

# Plot outcome
# Plot outcome
mito_ribo_old <- ggplot(results$qcSummary, aes(x = Original_RiboProp, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  theme_classic() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_ribo_new <- ggplot(results$qcSummary, aes(x = New_RiboProp, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_old <- ggplot(results$qcSummary, aes(x = Original_Features, y = Original_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

mito_features_new <- ggplot(results$qcSummary, aes(x = New_Features, y = New_MitoProp, colour = Damaged_Status)) + 
  scale_color_manual(values = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + 
  scale_y_continuous(labels = number_format(accuracy = 0.1)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")


# umap <- DimPlot(results, group.by = "orig.ident", cols = c("damaged" = "#FE3A55", "cell" = "#C1C1C2")) + NoAxes() + theme(panel.background = element_rect(fill = NULL, color = "black"))

((mito_features_old / mito_features_new) | (mito_ribo_old / mito_ribo_new))
