#' view_output
#'
#' Helper function for plotting
#'
#' @param data Seurat object resulting from add_damage damaged cells.
#' @param organism A character string for the organism of origin, either "human" or "mouse".
#' @param project_name A character string for the project name.
#' @param PBMC A boolean specifying whether the input data is of PBMC origin, default is FALSE.
#' @param group_by A character string specifying which annotated column should be used to colour the UMAP plot.
#' @param PC_num Numeric specifying how many principal components (PCs) should be included in the PCA plot, default = 50.
#' @param top_genes Numeric specifying how many genes should be considered for the PCA plot, default = 50.
#'
#' @return QC and UMAP plot
#' @export
#'
#' @examples
view_output <- function(
    data,
    organism = "human",
    project_name,
    PBMC = FALSE,
    group_by = "orig.ident",
    PC_num = 50,
    top_genes = 50
){
  # Preparations for plotting ----
  # Define colours of clusters
  colours <- c("0" = "#A6BEAE", "1" = "#88A1CD", "2" = "#A799C9", "3" = "#BDC5EE", "4" = "#9DBD73",
               "5" = "#A6BEAE", "6" = "#88A1CD", "7" = "#A799C9", "8" = "#BDC5EE", "9" = "#9DBD73",
               "10" = "#A6BEAE", "11" = "#88A1CD", "12" = "#A799C9", "13" = "#BDC5EE", "14" = "#9DBD73",
               "15" = "#A6BEAE", "16" = "#88A1CD", "17" = "#A799C9", "18" = "#BDC5EE", "19" = "#9DBD73",
               "20" = "#A6BEAE", "damaged" = "red", "cell" = "#E1E1E1")

  if (organism == "human"){ data("mito_genes", package = "damagesimulator")}

  if (organism == "mouse"){
    # Call mouse mito genes
    data("mouse_mito_genes", package = "damagesimulator")
    mito_genes <- mouse_mito_genes

  }

  # Define non-mitochondrial genes using mito genes
  non_mito_genes <- rownames(data@assays$RNA)[!(rownames(data@assays$RNA) %in% mito_genes)]


  # Plotting main UMAP clusters ----

  # Typical dimensionality reduction Seurat workflow
  data <- NormalizeData(data) %>%
    FindVariableFeatures() %>% # maybe this is where my interest should be, not in PCs?
    ScaleData() %>% # scale data according to non-mito. genes?
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.1) %>% # very low because I don't want many divisions within clusters, broad overview
    RunUMAP(dims = 1:30)

  # Generate plots (default Seurat plots)
  clusters <- DimPlot(data, group.by = group_by)  +
    NoAxes() + scale_color_manual(values = colours) + ggtitle(project_name) +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1))

  if (PBMC){

    # Define cell types
    celltypes <- FeaturePlot(data,
                             features = c("MS4A1", "CD3E", "NKG7", "CD14"), # Only for humans
                             cols = c("#E1E1E1", "#64A493"))

    for (plot in 1:length(celltypes)) {
      celltypes[[plot]] <- celltypes[[plot]] + NoLegend() + NoAxes() +
        theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
              plot.title = element_text(size = 12))
    }

    # Combine with general UMAP labeling damaged and non-damaged cells
    cluster_plot <- clusters | celltypes

  } else {

    # Else just output the general UMAP
    cluster_plot <- clusters

  }


  # QC Metrics plots ----

  # Pre-calculated metrics
  mean_counts <- colMeans(data@assays$RNA$counts)
  var_counts <- colVars(data@assays$RNA$counts)
  mito_percent <- data$mt.percent
  ribo_percent <- data$rb.percent
  cells <- rownames(data@meta.data)

  # Manual metrics
  mito_counts <- colSums(data@assays$RNA$counts[rownames(data@assays$RNA$counts) %in% mito_genes, ]) # Calculate mitochondrial counts
  zero_counts <- colSums(data@assays$RNA$counts == 0) # Calculate zero counts (genes with 0 expression per cell)
  feature_counts <- rep(dim(data)[1], length(zero_counts))   # Calculate total features per cell (total number of genes in the dataset)
  dropout <- zero_counts / feature_counts    # Calculate dropout as a proportion of zero counts over total feature counts (between 0 and 1)

  # Create a data frame with all the metrics
  df <- data.frame(
    Cell = cells,
    Total_counts = data$nCount_RNA,
    Mean_counts = mean_counts,
    Features = data$nFeature_RNA,
    Variance = var_counts,
    Mito_percent = mito_percent,
    Ribo_percent = ribo_percent,
    Mito_counts = mito_counts,
    Zero_counts = zero_counts,
    Dropout = dropout
  )


  df <- df[df$Variance <= 150, ]

  # Existing metrics plot
  count_variance <- ggplot(df, aes(x = Mean_counts, y = Variance)) +
    #geom_point() +
    #geom_smooth(method=lm) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Cell Count") +
    geom_bin2d(bins = 30) + # granularity
    theme_classic() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  dropout_count <- ggplot(df, aes(x = Total_counts, y = Dropout)) +
    #geom_point() +
    #geom_smooth(method=lm) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Cell Count") +
    geom_bin2d(bins = 30) + # granularity
    theme_classic() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  mt_feature <- ggplot(df, aes(x = Features, y = Mito_percent)) +
    #geom_point() +
    #geom_smooth(method=lm) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Cell Count") +
    geom_bin2d(bins = 30) + # granularity
    theme_classic() +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


  mt_rb <- ggplot(df, aes(x = Ribo_percent, y = Mito_percent)) +
    #geom_point() +
    #geom_smooth(method=lm) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Cell Count") +
    geom_bin2d(bins = 30) + # granularity
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

  # PC plot (how much mito. genes contribute to PCA)
  pc_loadings <- data@reductions$pca@feature.loadings
  pc_stdev <- data@reductions$pca@stdev
  mito_gene_indices <- which(rownames(pc_loadings) %in% mito_genes)
  mito_percentages <- numeric(PC_num)

  for (pc_idx in seq_len(PC_num)) {
    top_genes_indices <- order(abs(pc_loadings[, pc_idx]), decreasing = TRUE)[1:top_genes]
    mito_percentages[pc_idx] <- sum(top_genes_indices %in% mito_gene_indices) / top_genes * 100
  }

  pc_data <- data.frame(
    PC = seq_len(PC_num),
    Percent.mito = mito_percentages,
    Standard_Deviation = pc_stdev[seq_len(PC_num)]
  )

  # Calculate scaled Standard Deviation
  y_max <- max(25, max(pc_data$Percent.mito, na.rm = TRUE))
  pc_data$Scaled_Standard_Deviation <- pc_data$Standard_Deviation * y_max / max(pc_data$Standard_Deviation, na.rm = TRUE)


  # Create the plot
  mito_pc_plot <- ggplot(pc_data, aes(x = PC)) +
    geom_bar(aes(y = Percent.mito, fill = "Mitochondrial %"), stat = "identity", color = "black") + # Mitochondrial % as bars
    geom_point(aes(y = Scaled_Standard_Deviation, color = "Standard Deviation"), size = 2) + # Standard Deviation as points
    scale_y_continuous(
      limits = c(0, y_max),
      name = "Mitochondrial Contribution (%)",
      sec.axis = sec_axis(~ . * max(pc_data$Standard_Deviation, na.rm = TRUE) / y_max,
                          name = "Standard Deviation")
    ) +
    scale_fill_manual(values = c("Mitochondrial %" = "lightgrey")) +
    scale_color_manual(values = c("Standard Deviation" = "blue")) +
    labs(
      x = "Principal Component",
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  mito_pc_plot


  # Combine QC plots
  plots <-  (count_variance | dropout_count | mt_feature | mt_rb ) / mito_pc_plot

  return(list(clusters = cluster_plot,
              metrics = plots))

}
