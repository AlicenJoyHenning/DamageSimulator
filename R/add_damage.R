#' add_damage
#'
#' @param data ScRNA-seq count matrix in the form of a matrix, `Seurat` object, or `SingleCellExperiment` object.
#' @param damage_fraction Numeric of the proportion of damage to be introduced to the sample, default is 0.20.
#' @param beta_proportion Numeric of the proportion of damaged cells belonging to the beta distribution, distribution with a high manifestation of damage, default is 0.5.
#' @param damage_distribution_shapes Vector of four shape parameters for two beta distributions, default is c(a1 = 3, b1 = 12, a2 = 12, b2 = 3).
#' @param damage_scale_factor Numeric value multiplied by the damage level assigned to a cell to generate the damage perturbation factor, which is then applied to all cytoplasmic genes, default is 0.15.
#' @param organism A character string for the organism of origin, either "human" or "mouse".
#' @param output_format A character string specifying what data class the output count matrix be give, either "Seurat", "SingleCellExperiment" or "matrix", default matches the input class.
#' @param generate_plots A boolean value specifying whether to plot QC metrics and clusters of the resulting damage simulated data, default is FALSE.
#' @param view_PBMC_markers A boolean value specifying whether input data is of PBMC origin and marker genes can be viewed
#' @param project_name A character string for the project name.
#' @param save_path A character string specifying the output directory.
#'
#' @return Single cell count matrix with a proportion of damage perturbed cells housed in either a "Seurat", "SingleCellExperiment" or "matrix" object.
#' @export
#'
#' @examples
#' # Load example data
#' library(scRNAseq)
#' pbmc <- fetchDataset("he-organs-2020", "2023-12-21", path = "blood")
#'
#' # Run the
#' damaged_PBMC <- add_damage(
#'   data = pbmc,
#'   damage_fraction = 0.1,
#'   organism = "human",
#'.  output_format = "SingleCellExperiment",
#'   project_name = "test
#'  )
#'


data <- readRDS("/Users/alicen/Projects/ReviewArticle/automating-workflow/extracted_data/datasets/non-groundtruth/hPBMC_processed.rds")


damaged_PBMC_0.5 <- add_damage(
  data = data,
  damage_fraction = 0.5,
  organism = "human",
  output_format = "Seurat",
  generate_plots = TRUE,
  view_PBMC_markers = TRUE,
  project_name = "test_data_0.5",
  save_path = "~/Projects/ReviewArticle/"
)

mt <- as.data.frame(damaged_PBMC_0.1$mt.percent)
mt_0.8 <- as.data.frame(damaged_PBMC_0.8$mt.percent)
plot_0.1 <- ggplot(mt, aes(y = mt$`damaged_PBMC_0.1$mt.percent`, x = "")) + geom_boxplot()
plot_0.8 <- ggplot(mt, aes(y = mt_0.8$`damaged_PBMC_0.8$mt.percent`, x = "")) + geom_boxplot()

plot_0.1 | plot_0.8

add_damage <- function(
    data,
    damage_fraction = 0.2,
    beta_proportion = 0.5,
    damage_distribution_shapes = c(a1 = 3, b1 = 12, a2 = 12, b2 = 3),
    damage_scale_factor = 0.15,
    organism = "human",
    output_format = "Seurat",
    generate_plots = FALSE,
    view_PBMC_markers = FALSE,
    project_name = "DamageSimulation",
    save_path = "./"
) {

  # Reproducibility
  set.seed(42)

  # Identify input data type and prepare accordingly
  if (class(data) == "Seurat"){

    # Step 1: Select cells to perturb
    total_cell_number <- ncol(data@assays$RNA$counts)
    damaged_cell_number <- round(damage_fraction * total_cell_number)
    damaged_cell_selections <- sample(colnames(data@assays$RNA$counts), damaged_cell_number)

    # Isolate into two matrices
    to_be_perturbed <- subset(data, cells = damaged_cell_selections)
    to_be_perturbed <- as.matrix(to_be_perturbed@assays$RNA$counts)
    unperturbed <- subset(data, cells = setdiff(colnames(data), damaged_cell_selections))
    unperturbed <- as.matrix(unperturbed@assays$RNA$counts)

  }

  if (class(data) == "matrix"){

    # Step 1: Select cells to perturb
    total_cell_number <- ncol(data)
    damaged_cell_number <- round(damage_fraction * total_cell_number)
    damaged_cell_selections <- sample(colnames(data), damaged_cell_number)

    # Isolate into two matrices
    to_be_perturbed <- data[, damaged_cell_selections]
    unperturbed <- data[, setdiff(colnames(data), damaged_cell_selections)]

  }

  if (class(data) == "SingleCellExperiment"){

    # Step 1: Select cells to perturb
    total_cell_number <- ncol(data)
    damaged_cell_number <- round(damage_fraction * total_cell_number)
    damaged_cell_selections <- sample(colnames(data), damaged_cell_number)

    # Isolate into two matrices
    counts_matrix <- data@assays@data$counts
    colnames(counts_matrix) <- colnames(data)
    rownames(counts_matrix) <- rownames(data)
    to_be_perturbed <- counts_matrix[, damaged_cell_selections]
    unperturbed_cells <- setdiff(colnames(data), damaged_cell_selections)
    unperturbed <- counts_matrix[, unperturbed_cells]

    # CHECKS :
    # dim(to_be_perturbed)[2] == length(damaged_cell_selections)
    # dim(unperturbed)[2] == length(unperturbed_cells)

  }


  # Step 2: Generate damage scaling values
  # Using two beta distributions (two distinct populations, outcome of exponential relationship between. damage level (x) and physical manifestation of damage (y)) i.e. once a threshold is reached damage becomes 'visible'
  high_count <- round(damaged_cell_number * beta_proportion)
  low_count <- damaged_cell_number - high_count # Ensure proportion maintained


  # Draw damage_levels from the distribution given it satisfies the constraints
  if (length(damage_distribution_shapes) == 4) {

    # Check the shape parameters satisfy the target distribution
    unimodal_constraint <- (damage_distribution_shapes[["a1"]] > 1) &
      (damage_distribution_shapes[["b1"]] > 1) &
      (damage_distribution_shapes[["a2"]] > 1) &
      (damage_distribution_shapes[["b2"]] > 1)

    bimodal_constraint <- (damage_distribution_shapes[["a1"]] < damage_distribution_shapes[["b1"]]) &
      (damage_distribution_shapes[["a2"]] > damage_distribution_shapes[["b2"]])

    if (unimodal_constraint & bimodal_constraint){
      # Generate the values for the distribution
      low_scaling <- rbeta(low_count, shape1 = damage_distribution_shapes[["a1"]], shape2 = damage_distribution_shapes[["b1"]])
      high_scaling <- rbeta(high_count, shape1 = damage_distribution_shapes[["a2"]], shape2 = damage_distribution_shapes[["b2"]])
      damage_levels <- c(low_scaling, high_scaling)
      #damaged_cell_number == length(damage_levels) # TRUE

    } else {
      stop("Please ensure damage_distribution_shapes satisfy constraints.")
    }

  } else {
    stop("Please ensure damage_distribution_shapes has 4 entries.")
  }

  # Step 3: Perform perturbation

  # Retrieve mitochondrial data
  if (organism == "human"){

    mito_genes <- get("human_mito_genes", envir = asNamespace("damagesimulator"))

  } else if (organism == "mouse") {

    mito_genes <- get("mouse_mito_genes", envir = asNamespace("damagesimulator"))

  }

  # Identify non-mitochondrial genes (these expression values will be altered)
  mito_idx <- rownames(to_be_perturbed) %in% mito_genes
  non_mito_idx <- which(!mito_idx)

  # Extract mitochondrial and non-mitochondrial expression separately
  to_be_perturbed_non_mito <- to_be_perturbed[non_mito_idx, , drop = FALSE]
  #to_be_perturbed_mito <- to_be_perturbed[mito_idx, , drop = FALSE]

  # Define a lambda constant for scaling the damage levels (adjustable value)
  lambda <- damage_scale_factor  # Smaller lambda leads to less reduction

  # Generate a single perturbation factor per cell based on damage level, adjusted with lambda
  perturb_factors <- exp(-lambda * damage_levels)

  # Expand perturb_factors to match non-mito gene matrix dimensions
  perturb_factors_matrix <- matrix(rep(perturb_factors, each = nrow(to_be_perturbed_non_mito)),
                                   nrow = nrow(to_be_perturbed_non_mito),
                                   ncol = damaged_cell_number,
                                   byrow = FALSE)

  perturb_factors_matrix <- matrix(perturb_factors,
                                   nrow = nrow(to_be_perturbed_non_mito),
                                   ncol = length(perturb_factors),
                                   byrow = TRUE)


  # Apply the same reduction factor across all non-mito genes for each cell
  perturbed_expression <- to_be_perturbed_non_mito * perturb_factors_matrix


  # No change to mitochondrial counts
  perturbed_combined <- as.matrix(to_be_perturbed) # identical(colnames(perturbed_combined), damaged_cell_selections) # TRUE
  perturbed_combined[non_mito_idx, ] <- perturbed_expression


  # Step 4: Recombine with unperturbed data
  combined_matrix <- cbind(perturbed_combined, as.matrix(unperturbed))
  combined_matrix <- combined_matrix[rownames(unperturbed), ]
  combined_matrix <- as(combined_matrix, "dgCMatrix") # Essential for speed in view_output!
  combined_seurat <- CreateSeuratObject(counts = combined_matrix, assay = "RNA")

  # Adding metadata
  combined_seurat$orig.ident <- ifelse(rownames(combined_seurat@meta.data) %in% damaged_cell_selections, "damaged", "cell")

  if (organism == "human"){

    # Mitochondrial percent
    m_genes <- rownames(combined_seurat)[grep("^MT-", rownames(combined_seurat))]
    combined_seurat$mt.percent <- PercentageFeatureSet(
      object = combined_seurat,
      features = intersect(m_genes, rownames(combined_seurat@assays$RNA)),
      assay = "RNA"
    )

    # Ribosomal percent
    rb_genes <- rownames(combined_seurat)[grep("^RP[SL]", rownames(combined_seurat))]
    combined_seurat$rb.percent <- PercentageFeatureSet(
      object = combined_seurat,
      features = intersect(rb_genes, rownames(combined_seurat@assays$RNA)),
      assay = "RNA"
    )

  }

  if (organism == "mouse"){

    # Mitochondrial percent
    true_mito_genes <- mito_genes[grepl("^mt-", mito_genes)] # isolate mito genes from nuclear localised genes

    combined_seurat$mt.percent <- PercentageFeatureSet(
      object = combined_seurat,
      features = intersect(true_mito_genes, rownames(combined_seurat@assays$RNA)),
      assay = "RNA"
    )

    # Ribosomal percent
    rb_genes <- rownames(combined_seurat)[grep("^RP[SL]",  rownames(combined_seurat), ignore.case = TRUE)]
    combined_seurat$rb.percent <- PercentageFeatureSet(
      object = combined_seurat,
      features = intersect(rb_genes, rownames(combined_seurat@assays$RNA)),
      assay = "RNA"
    )

  }


  # Call helper function to generate diagnostic plots
  perturbed_by_reduction_plot <- view_output(combined_seurat, "Perturbed by reduction", "orig.ident", PBMC = view_PBMC_markers)
  metrics_plot <- perturbed_by_reduction_plot$metrics # 1000 x 480
  cluster_plot <- perturbed_by_reduction_plot$clusters # 800 x 380

  if (generate_plots == TRUE & is.na(save_path)) {
    stop("Please give save_path.")
  }

  if (generate_plots){

  ggsave(
    filename = file.path(save_path, paste0(project_name, "_metrics.png")),
    plot = metrics_plot,
    width = 3000, height = 1480, units = "px"
  )

  ggsave(
    filename = file.path(save_path, paste0(project_name, "_clusters.png")),
    plot = cluster_plot,
    width = 2400, height = 1100, units = "px"
  )
 }


  # Save the data ----
  if (!is.na(output_format)) {
    if (!output_format %in% c("Seurat", "SingleCellExperiment", "matrix")) {
      stop("Please ensure output_format is one of 'Seurat', 'SingleCellExperiment' or 'matrix'.")
    }
  } else {
    # If no final data type is specified, use the same that was used as input
    output_format <- class(data)[[1]]
  }

  # How to treat each type of output data
  if (output_format == "Seurat") {

    output <- combined_seurat

  }

  if (output_format == "matrix") {

    # Extract & return matrix from Seurat
    output <- as.matrix(combined_seurat@assays$RNA$counts)

  }

  if (output_format == "SingleCellExperiment") {

    # Convert to sce
    counts <- combined_seurat@assays$RNA
    sce <- SingleCellExperiment(list(counts = counts))
    rowData(sce)$featureType <- rownames(combined_seurat@assays$RNA)
    colData(sce)$subsets_mito_percent <- combined_seurat$mt.percent
    mainExpName(sce) <- 'gene'
    colData(sce)$detected <- combined_seurat$nFeature_RNA
    output <- sce

  }

  return(output)

}
