# Universal Heatmap Function with p-values
create_heatmap <- function(metabo_mat, pvalue, labels_group, output_path = "heatmap.pdf", g) {
  library(ComplexHeatmap)
  library(circlize)
  
  # Validate input
  if (missing(metabo_mat)) stop("Please provide the metabolomics matrix 'metabo_mat'.")
  if (missing(pvalue)) stop("Please provide a vector of p-values.")
  if (missing(labels_group)) stop("Please provide the group labels.")
  
  # Convert p-values to significance levels
  pch <- sapply(pvalue, function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p <= 0.05) return("*")
    else return(".")
  })
  
  # Define color gradient for p-values
  pvalue_col_fun <- colorRamp2(c(0.05, 0.01, 0.001), c("#f2e5f6", "#d9a8e3", "#9b6d9d"))
  
  # Define heatmap annotation
  ha <- HeatmapAnnotation(
    pvalue = anno_simple(pvalue, col = pvalue_col_fun, pch = pch, pt_gp = gpar(col = "black"), pt_size = unit(7, "mm")),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # Define heatmap colors
  base_colors <- c("#003333", "#f2e5f6", "#3c1518", "#d9a8e3", "#9b6d9d", "#ff8c61")
  color_palette <- colorRampPalette(base_colors)
  
  # Create the heatmap
  ht <- Heatmap(as.matrix(metabo_mat), name = "Intensity", col = color_palette(100),
                row_split = g, column_split = factor(colnames(metabo_mat)),
                row_gap = unit(5, "mm"), cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                border = TRUE, column_names_rot = 45, column_dend_side = "bottom",
                row_names_gp = gpar(fontsize = 15, fontface = "bold"),
                column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
                column_title = "Proteins", row_title = "Group", 
                top_annotation = ha,
                left_annotation = rowAnnotation(foo = anno_block(
                  gp = gpar(fill = c("#d77a61", "#d8b4a0", "#B56576")),
                  labels = labels_group,
                  labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold")
                ))
  )
  
  # Save the heatmap to a file
  pdf(output_path, width = 10, height = 8)
  draw(ht)
  dev.off()
}

# Example usage:
# create_heatmap(metabolomics_matrix, p_values, group_labels, "output_heatmap.pdf", factor(grouping_variable))






# ++++++++++++++++++++++++++++
# Generate a Heatmap for Metabolomics Data with p-values
# ++++++++++++++++++++++++++++
create_heatmap <- function(metabo_mat, pvalue, labels_group, output_path = "results/heatmap.pdf", g) {
  library(ComplexHeatmap)
  library(circlize)
  
  # Convert p-values to significance levels (***, **, *)
  pch <- sapply(pvalue, map_p_value_to_pch)
  
  # Define color gradient for p-values
  pvalue_col_fun <- colorRamp2(c(0.05, 0.01, 0.001), c("#f2e5f6", "#d9a8e3", "#9b6d9d"))
  
  # Define the heatmap annotation for p-values
  ha <- HeatmapAnnotation(
    pvalue = anno_simple(pvalue, col = pvalue_col_fun, pch = pch, pt_gp = gpar(col = "black"), pt_size = unit(7, "mm")),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # Define the heatmap color palette
  base_colors <- c("#003333", "#f2e5f6", "#3c1518", "#d9a8e3", "#9b6d9d", "#ff8c61")
  color_palette <- colorRampPalette(base_colors)
  
  # Create the heatmap
  ht <- Heatmap(as.matrix(metabo_mat), name = "Intensity", col = color_palette(100),
                row_split = g, column_split = factor(colnames(metabo_mat)),
                row_gap = unit(5, "mm"), cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                border = TRUE, column_names_rot = 45, column_dend_side = "bottom",
                row_names_gp = gpar(fontsize = 15, fontface = "bold"),
                column_names_gp = gpar(fontsize = 16, fontface = "bold"), 
                column_title = "Proteins", row_title = "Group", 
                top_annotation = ha,
                left_annotation = rowAnnotation(foo = anno_block(
                  gp = gpar(fill = c("#d77a61", "#d8b4a0", "#B56576")),
                  labels = labels_group,
                  labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold")
                ))
  )
  
  # Save the heatmap to a file
  pdf(output_path, width = 10, height = 8)
  draw(ht)
  dev.off()
}

# Example usage:
# create_heatmap(metabo_matrix, pval_vector, c("Group1", "Group2"), "heatmap_output.pdf", g = factor(1:10))

# ++++++++++++++++++++++++++++
# Generate Boxplots with Annotations and p-values
# ++++++++++++++++++++++++++++
generate_boxplot <- function(df_expression, df_pvalues = NULL, y_expand_value = 0.1, output_path = "results/boxplot", 
                             plot_title = "Boxplot", groups, plot_annotations = TRUE) {
  
  library(ggplot2)
  library(reshape2)
  
  # Convert data to long format
  long_df <- melt(df_expression, id = 'Group')
  colnames(long_df)[2] <- "Variable"
  colnames(long_df)[3] <- "Intensity"
  
  # Generate summary statistics
  df_summary <- long_df %>%
    group_by(Group, Variable) %>%
    summarise(len = mean(Intensity, na.rm = TRUE), se = sd(Intensity, na.rm = TRUE) / sqrt(n()), .groups = 'drop')
  
  # Merge p-values if available
  if (!is.null(df_pvalues)) {
    df_summary <- df_summary %>%
      left_join(df_pvalues %>% select(Gene_Name, adj.pvalues), by = c("Variable" = "Gene_Name"))
  }
  
  # Annotations for p-values
  if (plot_annotations && !is.null(df_pvalues)) {
    annotations <- df_summary %>%
      filter(Group %in% groups) %>%
      group_by(Variable) %>%
      summarise(
        y_max = max(len) + y_expand_value,
        pval = case_when(
          adj.pvalues < 0.001 ~ "***",
          adj.pvalues < 0.01 ~ "**",
          adj.pvalues < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
  }
  
  # Generate the plot
  p <- ggplot(df_summary, aes(x = Group, y = len, group = Variable)) +
    geom_line(aes(colour = Group, group = interaction(Variable, Group)), size = 1.3, linetype = "dotted") +
    geom_errorbar(aes(ymin = len - se, ymax = len + se), width = 0.3, size = 1.3, colour = "gray50") +
    geom_point(aes(colour = Group), size = 6, shape = 16) +
    scale_colour_manual(values = c("Group1" = "turquoise3", "Group2" = "violetred2")) +
    facet_wrap(~Variable) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "gray", size = 0.5, linetype = 2),
      axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 12)
    ) +
    ggtitle(plot_title)
  
  # Add annotations for p-values
  if (plot_annotations && !is.null(df_pvalues)) {
    p <- p + geom_text(data = annotations, aes(x = 1.5, y = y_max, label = pval), size = 5, vjust = -0.5)
  }
  
  # Save the plot
  ggsave(filename = paste0(output_path, ".png"), plot = p, width = 10, height = 6, dpi = 300)
  ggsave(filename = paste0(output_path, ".pdf"), plot = p, width = 10, height = 6)
  
  return(p)
}

# Example usage:
# generate_boxplot(df_expression, df_pvalues, 0.1, "results/boxplot_example", "Protein Expression", c("Group1", "Group2"))

# ++++++++++++++++++++++++++++
# Generate Combined Boxplots for Multiple Protein Categories
# ++++++++++++++++++++++++++++
generate_combined_boxplots <- function(data_aarhus, data_aalborg, protein_lists, protein_category_names, output_path = "results/combined_boxplot", labels) {
  
  library(ggplot2)
  library(dplyr)
  
  for (i in seq_along(protein_lists)) {
    selected_proteins <- protein_lists[[i]]
    category_name <- protein_category_names[[i]]
    
    for (protein in selected_proteins) {
      if (!(protein %in% colnames(data_aarhus)) || !(protein %in% colnames(data_aalborg))) {
        next  # Skip if the protein doesn't exist in both datasets
      }
      
      # Prepare data for Aarhus
      df_aarhus <- data_aarhus %>%
        select(Group, all_of(protein)) %>%
        mutate(Cohort = "Aarhus")
      
      # Prepare data for Aalborg
      df_aalborg <- data_aalborg %>%
        select(Group, all_of(protein)) %>%
        mutate(Cohort = "Aalborg")
      
      # Combine data
      df_combined <- bind_rows(df_aarhus, df_aalborg) %>%
        rename(Intensity = all_of(protein))
      
      # Generate the boxplot
      p <- ggplot(df_combined, aes(x = Group, y = Intensity, fill = Cohort)) +
        geom_boxplot() +
        scale_fill_manual(values = c("Aarhus" = "red", "Aalborg" = "blue")) +
        theme_minimal() +
        labs(
          title = paste0("Protein: ", protein, " (", category_name, ")"),
          x = "Group",
          y = "Intensity"
        )
      
      # Save the plot
      ggsave(paste0(output_path, "_", category_name, "_", protein, ".png"), plot = p, width = 10, height = 6, dpi = 300)
    }
  }
}

# Example usage:
# generate_combined_boxplots(df_aarhus, df_aalborg, protein_lists, protein_category_names, "results/combined_boxplots")

