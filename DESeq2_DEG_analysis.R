# DESeq2 for DEGs analysis and making plots in R
rm(list = ls())

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tibble)
library(dplyr)
library(ggpubr)
library(ggrepel)

work_dir <- "RNAseq/DESeq2_results3"
if (!file.exists(work_dir)) {dir.create(work_dir, recursive = TRUE)}
setwd(work_dir)

directory <- "RNAseq/htseq_results"

# in the MA plot, label out the NK2 AAEL003794, KAT8 AAEL007672, PP6 AAEL026630, RMT1 AAEL001037

# Loop for groups 1-6 (PE72), 7-12 (PBM24), 13-18 (PBM60), 19-24 (PBM36)
for(i in seq(1, 24, by = 6)) {

  if (i == 1) {
    plot_title <- "PE72 WT vs KO"
  } else if (i == 7) {
    plot_title <- "PBM24 WT vs KO"
  } else if (i == 13) {
    plot_title <- "PBM60 WT vs KO"
  } else if (i == 19) {
    plot_title <- "PBM36 WT vs KO"
  } else {
    plot_title <- "error: not correct number ID"
  }
  
  group_files <- paste0(seq(i, i+5), "_counts.txt") %>% as.factor()
  group_conditions <- c(rep("WT", 3), rep("KO", 3)) %>% as.factor()
  
  # Create a data frame for one group
  group_sample_table <- data.frame(sampleName = group_files,
                                   fileName = group_files,
                                   condition = group_conditions)
  
  group_sample_table$condition <- factor(group_sample_table$condition, levels = c("WT","KO"))
  
  # DESeq2 analysis
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = group_sample_table,
                                         directory = directory,
                                         design= ~ condition)
  
  dds <- DESeq(ddsHTSeq)
  
  #Note on factor levels
  dds$condition <- relevel(dds$condition, ref = "WT")
  
  #Differential expression analysis
  res <- results(dds)
  
  #Log fold change shrinkage for visualization and ranking
  resLFC <- lfcShrink(dds, coef="condition_KO_vs_WT", type="apeglm")
  # Exclude NAs and filter by padj
  resLFC <- resLFC[!is.na(resLFC$padj),]
  
  # PCA Plot Preparation vst: Variance Stabilizing Transformation
  vsd_group <- vst(dds, blind = FALSE)
  pca_plot <- plotPCA(vsd_group, intgroup = "condition", returnData = F) +
    scale_y_continuous(limits = c(-5, 5))

  # Set the filename for the PCA plot
  pca_plot_filename <- paste0("PCA_plot_sample_", plot_title, ".png")
  
  # Open a new png device to save the plot
  png(file = pca_plot_filename, width = 1600, height = 1200, res = 300)
  plot(pca_plot)
  dev.off()

  print("Plot generation completed.")
  print(paste0("Saved PCA plot to", pca_plot_filename))
  
  topVarGenes <- head(order(rowVars(assay(vsd_group)), decreasing = TRUE), 20)
  mat  <- assay(vsd_group)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd_group)[, c("condition", "sizeFactor")]) # The sizeFactor is 
  
  heatmap_plot <- pheatmap(mat, annotation_col = anno, show_colnames = F)
  
  # Set the filename for the heatmap plot
  heatmap_plot_filename <- paste0("Heatmap_plot_sample_", plot_title, ".png")
  
  # Open a new png device to save the plot
  png(file = heatmap_plot_filename, width = 1600, height = 1600, res = 300)
  print(heatmap_plot)
  dev.off()
  
  res <- as.data.frame(res)
  res_filtered <- res %>% filter(!is.na(padj))
  
  # Define a named vector with custom gene labels
  myGenes <- c(AAEL003794 = "NK2", 
               AAEL007672 = "KAT8", 
               AAEL026630 = "PP6", 
               AAEL001037 = "RMT1")
  # Check if rownames are set properly; if not, set them (this step assumes gene IDs are in a column)
  if (is.null(rownames(res_filtered))) {
    res_filtered$gene_id <- as.character(res_filtered$gene_id)  # Ensure gene IDs are characters, not factors
    rownames(res_filtered) <- res_filtered$gene_id}
  
  # Update gene names in the dataset for the plot
  # This step maps gene IDs to custom names only for those specified in 'myGenes'
  res_filtered$display_name <- ifelse(rownames(res_filtered) %in% names(myGenes), myGenes[rownames(res_filtered)], NA)
  
  # ggmaplot with updated labels with padj < 0.05
  ma_plot <- ggmaplot(res_filtered, main = plot_title,
           fdr = 0.05, fc = 2, size = 1,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           #genenames = as.vector(res_filtered$display_name),
           legend = "right",
           font.legend = c("bold", 14),
           font.main = c("bold", 18),
           top = 0,
           ggtheme = ggplot2::theme_minimal() +
             ggplot2::theme(axis.title.x = element_text(size = 16, face = "bold"),
                            axis.title.y = element_text(size = 16, face = "bold"),
                            axis.text.x = element_text(size = 14, face = "bold"),
                            axis.text.y = element_text(size = 14, face = "bold")))
  
  # Adding geom_text_repel from ggrepel package for labeling points
  ma_plot <- ma_plot +
    geom_label_repel(size = 6, # adjust the label font size
                     fontface="bold",
                     data = res_filtered, # Applying label filtering here
                     aes(label = as.vector(res_filtered$display_name),
                         x = log2(baseMean), y = log2FoldChange),
                     min.segment.length = 0,
                     segment.color = 'red')
  # Customize the y-axis
  #ma_plot <- ma_plot + scale_y_continuous(limits = c(-5, 5), 
  #                                       breaks = seq(-5, 5, by = 1))
  
  # Set the filename for the MA plot
  ma_plot_filename <- paste0("MA_plot_sample_", plot_title, ".png")
  
  # Open a new png device to save the MA plot
  png(file = ma_plot_filename, width = 2700, height = 1800, res = 300)
  
  print(ma_plot)
  dev.off()
  
  # Add protein description to res_filtered table
  GeneAddProtein_df <- 
    read.delim("/Users/iEcR_RNAseq/GeneAddProtein.txt", sep = "\t")
  GeneAddProtein_df <- GeneAddProtein_df[,c(1,2,3)]
  
  res_filtered$GeneID <- rownames(res_filtered)
  res_filtered_addpro_df <- merge(GeneAddProtein_df, res_filtered, by = "GeneID", all.x = TRUE)
  res_filtered_addpro_df <- res_filtered_addpro_df %>% filter(!is.na((padj)), padj < 0.05)
  res_filtered_addpro_dfname <- paste0("DESeq2_res_", plot_title, ".csv")
  
  write.csv(res_filtered_addpro_df, file = res_filtered_addpro_dfname, row.names = F)
}
