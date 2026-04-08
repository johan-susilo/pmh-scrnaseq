suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(Seurat)
})

#' @title Safe Save Plot
#' @description Wraps plot saving in a tryCatch to prevent the pipeline from 
#' crashing if a graphical device error occurs. Saves BOTH PDF and PNG.
#' @param plot_obj The ggplot object.
#' @param base_filepath String. Path to save the file (without extension).
#' @param w Numeric. Width in inches.
#' @param h Numeric. Height in inches.
#' @export
safe_save_plot <- function(plot_obj, base_filepath, w = 15, h = 15) {
  pdf_path <- paste0(base_filepath, ".pdf")
  tryCatch({
    pdf(pdf_path, width = w, height = h)
    print(plot_obj)
    dev.off()
  }, error = function(e) {
    if (length(dev.list()) > 0) dev.off()
  })

  png_path <- paste0(base_filepath, ".png")
  tryCatch({
    png(png_path, width = w, height = h, units = "in", res = 300)
    print(plot_obj)
    dev.off()
  }, error = function(e) {
    if (length(dev.list())> 0) dev.off()
  })

}

#' @title Create QC Violin Plot
#' @description Generates a robust violin plot for QC metrics.
#' @param seurat_obj A Seurat object.
#' @param features Character vector of features to plot.
#' @param title String for the main plot title.
#' @return A ggplot object.
#' @export
create_qc_violin_plot <- function(seurat_obj, features, title) {
  
  # figure out if the cells have been grouped yet.
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cluster_id <- "seurat_clusters"
  } else if (length(levels(Idents(seurat_obj))) > 1){
    seurat_obj@meta.data$temp_ident <- Idents(seurat_obj)
    cluster_id <- "temp_ident"
  } else {
    cluster_id <- NULL
  }

  # this runs if we DO have clusters
  if (!is.null(cluster_id)) {

    # ggplot needs data in a "long" format, so that
    # all the values are stacked in one column called 'value'.
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(c(features, cluster_id))) %>%
      mutate(cell = rownames(seurat_obj@meta.data)) %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')

    # we standardize the column name to "Identity" so ggplot knows
    # exactly what to look for, and convert it to a factor so clusters plot in numerical order.
    colnames(qc_data)[colnames(qc_data) == cluster_id] <- "Identity"
    qc_data$Identity <- as.factor(qc_data$Identity)

    # this dynamically generates enough unique colors based on the number of clusters
    n_clusters <- length(unique(qc_data$Identity))
    colors <- if(n_clusters <= 12) colorRampPalette(brewer.pal(min(n_clusters, 12), "Paired"))(n_clusters) else colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)

    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = Identity)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.1, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      scale_fill_manual(values = colors) +
      theme_bw(base_size = 12) +
      theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
      labs(title = title, x = 'Identity', y = 'Value')

  } else {

    # hardcode a new column called "Identity" and set every single cell to "All"
    qc_data <- seurat_obj@meta.data %>%
      select(all_of(features)) %>%
      mutate(cell = rownames(seurat_obj@meta.data), Identity = "All") %>%
      pivot_longer(cols = all_of(features), names_to = 'metric', values_to = 'value')


    # here, we map 'fill' to 'metric' (not Identity) so the different
    # violins (nCount, nFeature) are colored differently, since there is only one "All" group.
    p <- ggplot(qc_data, aes(x = Identity, y = value, fill = metric)) +
      geom_violin(trim = FALSE, scale = "width") +
      geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
      facet_wrap(~metric, scales = 'free', ncol = length(features)) +
      theme_bw(base_size = 12) +
      theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
      labs(title = title, x = '', y = 'Value')
  }

  return(p)
}

#' @title Plot Cluster Proportions
#' @description Creates a stacked barplot of cluster proportions.
#' @param seurat_obj A Seurat object.
#' @param group_by The metadata column to group by (e.g., "orig.ident1").
#' @param custom_colors A vector of colors to use.
#' @return A ggplot object.
plot_cluster_proportions <- function(seurat_obj, group_by, custom_colors) {

  # calculate proportions dynamically using the provided group_by variable
  counts <- table(Idents(seurat_obj), seurat_obj[[group_by]][,1])
  props <- prop.table(counts, margin = 2)

  prop_df <- as.data.frame(props)
  colnames(prop_df) <- c("Cluster", "Group", "Proportion")
  prop_df <- prop_df %>% filter(Proportion > 0)

  p <- ggplot(prop_df, aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.2) + 
    geom_text(aes(label = as.character(Cluster), size = Proportion), 
              position = position_stack(vjust = 0.5), color = "black") +
    scale_size_continuous(range = c(0.5, 4), guide = "none") + 
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    ) + 
    labs(x = "Group", y = "Percentage of Total Cells") +
    scale_y_continuous(labels = scales::percent)


}