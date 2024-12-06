#' Introduce Perturbation to an Entire Count Matrix
#'
#' This function perturbs an entire count matrix by adjusting the proportion of mitochondrial genes in a subset of cells.
#'
#' @param seurat A Seurat object containing the single-cell RNA-seq data.
#' @param mito_genes A character vector of mitochondrial gene names.
#' @param percent_damage A numeric value indicating the percentage of cells to be perturbed.
#' @param project_name A character string for the project name.
#' @param output_dir A character string for the output directory (default is "/home/alicen/Projects/ReviewArticle/damage_perturbation/scDesign2").
#' @return A Seurat object with perturbed and unperturbed cells combined.
#' @examples
#' # Example usage:
#' # seurat <- CreateSeuratObject(counts = matrix(data = rpois(1000, lambda = 10), nrow = 100, ncol = 10))
#' # mito_genes <- c("MT-CO1", "MT-CO2", "MT-CO3")
#' # perturbed_seurat <- perturb_matrix(seurat, mito_genes, 20, "example_project")
#' @export
perturb_matrix <- function(seurat,
                             mito_genes,
                             percent_damage,
                             project_name,
                             output_dir) {

  # Extract percentage of cells at random, these will be altered and then reintroduced
  damaged_cell_number <- round(((percent_damage / 100) * (dim(seurat)[2])), 0)
  damaged_cell_selections <- sample(colnames(seurat), damaged_cell_number)
  to_be_perturbed <- subset(seurat, cells = damaged_cell_selections)
  to_be_perturbed_matrix <- to_be_perturbed@assays$RNA$counts

  # Remove these cells from the seurat object
  all_cells <- colnames(seurat)
  remaining_cells <- setdiff(all_cells, damaged_cell_selections)
  unperturbed <- subset(seurat, cells = remaining_cells)
  unperturbed_matrix <- unperturbed@assays$RNA$counts

  # Perturb the selected cells
  perturbed_matrix <- apply(to_be_perturbed_matrix, 2, perturb_single_cell, mito_genes = mito_genes)
  rownames(perturbed_matrix) <- rownames(to_be_perturbed_matrix)

  # Reorder the rows of perturbed_matrix to match the order of unperturbed_matrix
  perturbed_matrix <- perturbed_matrix[match(rownames(unperturbed_matrix), rownames(perturbed_matrix)), ]

  # Combine perturbed and unperturbed cells
  combined_matrix <- cbind(unperturbed_matrix, perturbed_matrix)
  combined_matrix <- as(combined_matrix, "dgCMatrix")
  combined_seurat <- CreateSeuratObject(counts = combined_matrix, assay = "RNA", project = paste0(seurat@project.name, "_perturbed"))

  # Edit the Seurat object to contain meaningful information
  combined_seurat$orig.ident <- ifelse(colnames(combined_seurat) %in% damaged_cell_selections, "damaged", "cell")

  # Calculate standard quality control metrics -----

  # Load default mito_genes if not provided
  data("rb_genes", package = "damagesimulator")


  # Calculate the feature percentages and normalized expressions
  combined_seurat$mt.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = intersect(mito_genes, rownames(combined_seurat@assays$RNA)),
    assay    = "RNA"
  )

  combined_seurat$rb.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = intersect(rb_genes, rownames(combined_seurat@assays$RNA)),
    assay    = "RNA"
  )

  combined_seurat$malat1.percent <- PercentageFeatureSet(
    object   = combined_seurat,
    features = "MALAT1",
    assay    = "RNA"
  )

  # Calculate pseudo-nuclear fraction score for DropletQC
  combined_seurat$nf_malat1 <- FetchData(combined_seurat, vars = "MALAT1", layer = "counts")
  combined_seurat$nf_malat1 <- (combined_seurat$nf_malat1 - min(combined_seurat$nf_malat1)) / (max(combined_seurat$nf_malat1) - min(combined_seurat$nf_malat1))

  # Transfer cell annotations
  combined_seurat$celltype <- seurat$celltype[match(rownames(combined_seurat@meta.data), rownames(seurat@meta.data))]
  combined_seurat$celltype <- ifelse(combined_seurat$orig.ident == "damaged", "damaged", combined_seurat$celltype)

  # Save outputs
  saveRDS(combined_seurat, file.path(output_dir, paste0(project_name, ".rds")))
  write.csv(combined_matrix, file.path(output_dir, paste0(project_name, "_matrix.csv")))

  # Visualize the output
  test <- NormalizeData(combined_seurat) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)

  # Note that these markers are applicable for immune/PBMCs
  colours <- c("T" = "#A6BEAE",
               "Monocyte" = "#88A1CD",
               "DC" = "#A799C9",
               "B" = "#BDC5EE",
               "NK" = "#9DBD73",
               "damaged" = "red"
  )

  plot <- DimPlot(test, group.by = "celltype")  +
    NoAxes() + ggtitle(project_name) + scale_color_manual(values = colours) +
    theme(panel.border = element_rect(colour = "black"))

  ggsave(filename = file.path(output_dir, paste0(project_name, "_perturbed.png")),
         plot,
         width = 5,
         height = 4,
         units = "in",
         dpi = 300)

  return(combined_seurat)
}
