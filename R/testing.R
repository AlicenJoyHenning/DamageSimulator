
# Load libraries ----
library(Matrix)       # For sparse matrices
library(ggplot2)      # For plotting
library(scales)       # For number formatting in ggplot2
library(patchwork)    # For combining multiple ggplot2 plots
library(scRNAseq)

# Damage simulation function ----

simulate_damage <- function(count_matrix, 
                            damage_proportion = 0.5, 
                            beta_proportion = 0.5, 
                            lambda = 10, 
                            damage_distribution_shapes = c(a1 = 18, b1 = 35, a2 = 40, b2 = 10)  #  Low ~ 0.05–0.15 High ~ 0.8–0.9
                            # (a + b) -> how steep, a/(a+b) where it sits      
                            ) {
  
  # Randomly select proportion of cells to damage
  total_cells <- ncol(count_matrix)
  damaged_cell_number <- round(total_cells * damage_proportion)
  damaged_cell_selections <- sample(seq_len(total_cells), size = damaged_cell_number, replace = FALSE)
  
  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)
  
  
  # Draw target damage level for the damaged cells
  high_count <- round(damaged_cell_number * beta_proportion)
  low_count  <- damaged_cell_number - high_count
  low_scaling  <- rbeta(low_count, shape1 = damage_distribution_shapes["a1"], shape2 = damage_distribution_shapes["b1"])
  high_scaling <- rbeta(high_count, shape1 = damage_distribution_shapes["a2"], shape2 = damage_distribution_shapes["b2"])
  damage_levels <- c(low_scaling, high_scaling)
  
  # View whether the shape is correct
  #test <- data.frame(damage_levels)
  #ggplot(test, aes(x = damage_levels)) + geom_histogram() + xlim(c(0, 1))
  #min(damage_levels) # 0.12482
  
  # Storage of damage levels for all cell barcodes for plotting later
  damage_label <- data.frame(barcode = colnames(count_matrix)[damaged_cell_selections], status = rep("damaged", length(damaged_cell_selections)))
  undamaged_cell_number_cells <- setdiff(seq_len(total_cells), damaged_cell_selections)
  undamaged_cell_number <- data.frame(barcode = colnames(count_matrix)[undamaged_cell_number_cells], status = rep("control", length(undamaged_cell_number_cells)))
  damage_label <- rbind(damage_label, undamaged_cell_number)
  

  # Isolate gene set indices (consistent across cells, not subsetting the matirx)
  mito_idx <- grep("^MT-", rownames(count_matrix), ignore.case = FALSE)
  ribo_idx <- grep("^(RPS|RPL)", rownames(count_matrix), ignore.case = FALSE)
  other_idx <- setdiff(seq_len(nrow(count_matrix)), union(mito_idx, ribo_idx))
  
  
  # Initialize for storing modified counts
  new_matrix <- count_matrix
  
  # Define the number of Monte Carlo samples
  n_samples <- 10000 # tested 1000 10000 100000 & found no helpful shift in shape from increasing
  
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
