rm(list=ls())

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(readxl)
library(ComplexHeatmap)
library(gridExtra)
library(forcats)
library(grid)
library(reshape2)
library(umap)

# Specify the directory path
dir_path <- "RNAseq/genes_clutering"

# Check if the directory exists, and create it if it doesn't
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
setwd(dir_path)
directory <- "RNAseq/htseq_results"

# Define group files and factors
group_files <- paste0(c(1,2,3,4,5,6,
                        7,8,9,10,11,12,
                        19,20,21,22,23,24,
                        13,14,15,16,17,18), "_counts.txt") %>% as.factor()

strain <- factor(c(rep("wt", 3), rep("FoxA_KO", 3),
                   rep("wt", 3), rep("FoxA_KO", 3),
                   rep("wt", 3), rep("FoxA_KO", 3),
                   rep("wt", 3), rep("FoxA_KO", 3)))

hour <- factor(c(rep("PE72", 3), rep("PE72", 3),
                 rep("PBM24", 3), rep("PBM24", 3),
                 rep("PBM36", 3), rep("PBM36", 3),
                 rep("PBM60", 3), rep("PBM60", 3)))

replicate <- factor(c(1,2,3,1,2,3,
                      1,2,3,1,2,3,
                      1,2,3,1,2,3,
                      1,2,3,1,2,3))

# Create a data frame for one group
group_sample_table <- data.frame(sampleName = group_files,
                                 fileName = group_files,
                                 strain = strain,
                                 hour = hour,
                                 replicate = replicate)

# DESeq2 analysis
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = group_sample_table,
                                       directory = directory,
                                       design= ~ strain + hour)

ddsHTSeq$hour <- relevel(ddsHTSeq$hour, ref = "PE72")
ddsHTSeq <- DESeq(ddsHTSeq, test = "LRT", reduced = ~ hour)
res <- results(ddsHTSeq)

# Normalize counts
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)

# Average the groups as previously done
group_sizes <- c(3, 3, 3, 3, 3, 3, 3, 3)
averages <- list()
current_col <- 1

for(size in group_sizes) {
  group_average <- rowMeans(normalized_counts[, current_col:(current_col + size - 1)], na.rm = TRUE)
  averages[[paste("Group", current_col, "-", current_col + size - 1)]] <- group_average
  current_col <- current_col + size
}

average_matrix <- do.call(cbind, averages)

# Scaling the rows
cal_z_score <- function(x){
  s <- sd(x, na.rm = TRUE)
  if(s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}

mt_norm <- t(apply(average_matrix, 1, cal_z_score))

# Use the differentially expressed (adjusted P < 0.05) across the 4 conditions by ANOVA
# intersect the 4 groups DEGs, WT VS KO for the four time points PE72, PBM24, PBM36, PBM60.
PBM24_DEGs_list <- read.csv("RNAseq/Deseq2_results2/DESeq2_res_PBM24 WT vs KO.csv")
PBM36_DEGs_list <- read.csv("RNAseq/Deseq2_results2/DESeq2_res_PBM36 WT vs KO.csv")
PBM60_DEGs_list <- read.csv("RNAseq/Deseq2_results2/DESeq2_res_PBM60 WT vs KO.csv")
PE72_DEGs_list <- read.csv("RNAseq/Deseq2_results2/DESeq2_res_PE72 WT vs KO.csv")

# Extract the first column
PBM24_genes <- unique(PBM24_DEGs_list[[1]])
PBM36_genes <- unique(PBM36_DEGs_list[[1]])
PBM60_genes <- unique(PBM60_DEGs_list[[1]])
PE72_genes <- unique(PE72_DEGs_list[[1]])

# Combine all gene lists into a single list
all_genes <- list(PBM24_genes, PBM36_genes, PBM60_genes, PE72_genes)

# Create a table that counts the occurrence of each gene across all lists
gene_count <- table(unlist(all_genes))

# Select genes that appear in at least three of the four lists
common_genes <- names(gene_count[gene_count >= 3])

# Print the number of common genes and the list of common genes
cat("Number of common genes:", length(common_genes), "\n")
print(common_genes)

mt_heatmap <- mt_norm[rownames(mt_norm) %in% common_genes, ]
colnames(mt_heatmap) <- c("PE72_wt", "PE72_FoxA_KO",
                          "PBM24_wt", "PBM24_FoxA_KO",
                          "PBM36_wt", "PBM36_FoxA_KO",
                          "PBM60_wt", "PBM60_FoxA_KO")

mt_heatmap_wt <- mt_heatmap[, c("PE72_wt", "PBM24_wt", "PBM36_wt", "PBM60_wt")]
mt_heatmap_ko <- mt_heatmap[, c("PE72_FoxA_KO", "PBM24_FoxA_KO", "PBM36_FoxA_KO", "PBM60_FoxA_KO")]

interested_genes <- c("AAEL003794", "AAEL026630")

# get the row index
# Get the indices of rows matching the interested genes
row_indices <- which(rownames(mt_heatmap) %in% interested_genes)
row_indices
# Create a row annotation object
row_anno <- rowAnnotation(foo = anno_mark(at = row_indices, 
                                          labels = interested_genes),
                          annotation_name_gp= gpar(fontsize = 20))

# Plot the heatmap with the filtered and normalized data
heatmap_plot <- Heatmap(mt_heatmap_wt,
                        name = "Global Heatmap",
                        cluster_rows = TRUE,
                        split = 10,
                        cluster_columns = FALSE,
                        column_labels = colnames(mt_heatmap_wt),
                        show_row_names = FALSE,
                        column_names_rot = 45,
                        right_annotation = row_anno)

plot_filename <- paste0("wt_heatmap.png")
png(file = plot_filename, width = 2000, height = 4000, res = 300)
print(heatmap_plot)
dev.off()

# PCA
# Read data
data <- mt_heatmap

# Perform PCA
pca_results <- prcomp(data, scale. = TRUE)

# Get the variance explained by each principal component
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2)
percentage_explained <- variance_explained * 100

# Extract the scores
scores <- as.data.frame(pca_results$x)

# Create a dataframe for plotting
plot_data <- data.frame(PC1 = scores[,1], PC2 = scores[,2])

# Generate the PCA plot
pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.4, color = "blue") +  # Adjust color and transparency as needed
  theme_minimal() +
  labs(title = sprintf("PCA of Gene Expression Data (PC1: %.2f%%, PC2: %.2f%%)",
                       percentage_explained[1], percentage_explained[2]),
       x = sprintf("Principal Component 1 (%.2f%%)", percentage_explained[1]),
       y = sprintf("Principal Component 2 (%.2f%%)", percentage_explained[2]))

# Print the plot
print(pca_plot)

# Save the plot as a PNG file
ggsave("pca_plot.png", pca_plot, width = 10, height = 8, dpi = 300)

# UMAP 
# Assuming 'mt_heatmap' is already loaded in your R session
umap_results <- umap::umap(mt_heatmap)

# Extract the UMAP coordinates
umap_data <- as.data.frame(umap_results$layout)

# Add row names as a column for labeling or identification
umap_data$Gene = row.names(umap_data)

# Create a UMAP plot
umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, label = Gene)) +
  geom_point(alpha = 0.6, color = "blue") +  # Adjust color and transparency as needed
  theme_minimal() +
  labs(title = "UMAP Projection of Gene Expression Data",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2")

# Print the plot
print(umap_plot)

# Optionally save the plot
ggsave("umap_plot.png", umap_plot, width = 10, height = 8, dpi = 300)


# Plot the heatmap with the filtered and normalized data
heatmap_plot <- Heatmap(mt_heatmap_ko,
                        name = "Global Heatmap",
                        cluster_rows = TRUE,
                        split = 10,
                        cluster_columns = FALSE,
                        column_labels = colnames(mt_heatmap_ko),
                        show_row_names = FALSE,
                        column_names_rot = 45,
                        right_annotation = row_anno)

plot_filename <- paste0("ko_heatmap.png")
png(file = plot_filename, width = 2000, height = 4000, res = 300)
print(heatmap_plot)
dev.off()

row_anno <- rowAnnotation(foo = anno_mark(at = row_indices, 
                                          labels = interested_genes,
                                          labels_gp = gpar(fontsize = 20)))
# Plot the heatmap with the filtered and normalized data
heatmap_plot <- Heatmap(mt_heatmap,
                        name = "Z-score",
                        cluster_rows = TRUE,
                        split = 10,
                        cluster_columns = FALSE,
                        column_labels = colnames(mt_heatmap),
                        show_row_names = FALSE,
                        column_names_rot = 45,
                        right_annotation = row_anno,
                        row_names_gp = gpar(fontsize = 20),  # Increase row names text size
                        column_names_gp = gpar(fontsize = 20),  # Increase column names text size
                        heatmap_legend_param = list(legend_direction = "horizontal", 
                                                    labels_gp = gpar(fontsize = 15)  # Increase legend labels text size
                        ))

plot_filename <- paste0("all_heatmap.pdf")
pdf(file = plot_filename, width = 8, height = 10)
draw(heatmap_plot,
     heatmap_legend_side = "bottom")
dev.off()


library(ComplexHeatmap)
library(grid)

# Lists of transcription activation and repression genes with their names
interested_genes <- c("AAEL007672", "AAEL001037", "AAEL026630", "AAEL003794")
interested_names <- c("KAT8", "RMT1", "PP6", "NK2")

# Get indices of these genes in the heatmap matrix
all_row_indices <- which(rownames(mt_heatmap) %in% interested_genes)

# Filter the genes and names to include only those present in the heatmap
present_genes <- rownames(mt_heatmap)[all_row_indices]
present_labels <- interested_names[match(present_genes, interested_genes)]

# Concatenate gene IDs with names
present_labels <- paste(present_genes, present_labels)

# Create row annotation object with geneID and gene name
row_anno <- rowAnnotation(foo = anno_mark(at = all_row_indices, 
                                          labels = present_labels,
                                          labels_gp = gpar(fontsize = 20)),
                          annotation_name_gp = gpar(fontsize = 20))

# Plot the heatmap with the filtered and normalized data
heatmap_plot <- Heatmap(mt_heatmap,
                        name = "Z-score",
                        cluster_rows = TRUE,
                        split = 7,
                        cluster_columns = FALSE,
                        column_labels = colnames(mt_heatmap),
                        show_row_names = FALSE,
                        column_names_rot = 45,
                        right_annotation = row_anno,
                        row_names_gp = gpar(fontsize = 20),  # Increase row names text size
                        column_names_gp = gpar(fontsize = 20),  # Increase column names text size
                        heatmap_legend_param = list(legend_direction = "horizontal", 
                                                    labels_gp = gpar(fontsize = 15)  # Increase legend labels text size
                        ))

plot_filename <- paste0("all_heatmap.pdf")
pdf(file = plot_filename, width = 8, height = 10)
draw(heatmap_plot,
     heatmap_legend_side = "bottom")
dev.off()


# Check which repression genes are present in the rownames of mt_heatmap
present_repression_genes <- repression_genes[repression_genes %in% rownames(mt_heatmap)]

# Print the results
print(present_repression_genes)

cluster_ids <- row_order(heatmap_plot)

for (i in seq_along(cluster_ids)) {
  # Extract the gene IDs for the current cluster
  gene_ids <- rownames(mt_heatmap)[cluster_ids[[i]]]
  
  # Create a filename for the cluster
  filename <- paste0("cluster_", i, "_genes.csv")
  
  # Save the gene IDs to a CSV file without the first row as "x"
  write.table(gene_ids, file = filename, row.names = FALSE, col.names = FALSE,
              quote = FALSE, sep = ",")
}



# TO DO
# 1. separate the wt and KO to two clustering plots
# 2. use all DEGs to make the clustering plots.
# 3. for each clustering set, do the GO enrichemnt, and the boxplot
# of average genes expression level.
# 4. label the important genes in the corresponding clustering set.
# 5. Use clusterProfiler to add GO enrichment results to the heatmap.

# GO enrichment
# Install and load the necessary packages if you haven't already
library(wordcloud)
library(RColorBrewer)

# Assuming go_list is already loaded with 7 GO enrichment results

# Define the cluster sizes as per the given ratios
cluster_sizes <- c(530, 692, 805, 511, 417, 585, 732)

# Loop through each GO enrichment result, create a word cloud, and save as a PDF
for (i in seq_along(go_list)) {
  # Extract GO terms and calculate their weights
  go_terms <- go_list[[i]]$Name
  go_weights <- -log10(go_list[[i]]$P.value)
  
  # Define the filename for the PDF
  png_filename <- paste0("word_cloud_cluster_", i, ".png")
  par(mar = c(0, 0, 0, 0))
  # Open the PDF device
  png(file = png_filename, width = 20, height = 15, units = "in", res = 300)
  
  # Generate the word cloud
  set.seed(123)  # For reproducibility
  wordcloud(
    words = go_terms,          # The GO term names
    freq = go_weights,         # The weights (importance) of the terms
    min.freq = 3,              # Minimum frequency for a term to be included
    max.words = 6,           # Maximum number of words in the cloud
    random.order = TRUE,      # Place the most frequent words at the center
    rot.per = 0,            # Percentage of words with 90-degree rotation
    colors = brewer.pal(8, "Dark2")  # Color palette for the words
  )
  
  # Close the PDF device
  dev.off()
}

# Methods
# Upload gene list such as Down_18PBM.txt to the website VectorBase 
# Then perform the GO enrichment with online tools.
# Download the BP, CC, MF three results and save them as, for example,
# BP_18PBM_UP.tsv, CC_18PBM_UP.tsv, MF_18PBM_UP.tsv, 
# BP_18PBM_DOWN.tsv, CC_18PBM_DOWN.tsv, MF_18PBM_DOWN.tsv
# For 18PBM MF UP first
# red uses #ECA8A9, blue uses #74AED4
# Make the 96h plots first

bp_data <- read.delim("cluster6_bp.tsv", sep = "\t") %>% 
  filter(P.value < 0.05, Benjamini < 0.05, Bonferroni < 0.05) %>% mutate(class = "BP")
cc_data <- read.delim("cluster6_cc.tsv", sep = "\t") %>% 
  filter(P.value < 0.05, Benjamini < 0.05, Bonferroni < 0.05) %>% mutate(class = "CC")
mf_data <- read.delim("cluster6_mf.tsv", sep = "\t") %>% 
  filter(P.value < 0.05, Benjamini < 0.05, Bonferroni < 0.05) %>% mutate(class = "MF")
    
combined_data <- bind_rows(bp_data, cc_data, mf_data)
    
# Order the data by class and then by Result.count within each class
combined_data <- combined_data %>%
      group_by(class) %>%
      top_n(7, Result.count) %>%
      ungroup() %>%
      arrange(class, desc(Result.count)) %>%
      mutate(Name = factor(Name, levels = unique(Name)))
    
# Set the space between the bar number label and the border automatically
# Automatically adjust space based on the maximum value
max_value <- max(combined_data$Result.count)
space_needed <- max_value * 0.1
    
    # Now plot using the adjusted combined_data
GO_plot <- ggplot(combined_data, aes(x = Name, y = Result.count, fill = class)) +
    geom_bar(stat = "identity") +
    coord_flip() + # Make it horizontal
    geom_text(aes(label = Result.count), hjust = -0.1, fontface = "bold") + # Adjust text position
    scale_fill_manual(values = c("BP" = "#F8766D", "CC" = "#00BA38", "MF" = "#619CFF"),
                      name = "GO class") +
    labs(x = "GO Term", y = "Gene Count", title = paste("Cluster 6 GO Enrichment", sep = " ")) +
    theme(
      panel.background = element_blank(), # Remove panel background
      panel.grid.major = element_line(colour = "lightgrey"), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      axis.title.x = element_text(size = 16, color = "black", face = "bold"),
      axis.title.y = element_text(size = 16, color = "black",face = "bold"),
      axis.text.x = element_text(size = 16, color = "black",face = "bold"),
      axis.text.y = element_text(size = 16, color = "black",face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Make plot title bold and centered
      legend.title = element_text(face = "bold"), # Make legend title bold
      legend.text = element_text(face = "bold"),
      ) +
      scale_y_continuous(expand = expansion(add = c(0, space_needed)))
    
    # Set the filename for the plot
    GO_plot_filename <- c("GO_plot_cluster6.png")
    png(file = GO_plot_filename, width = 3300, height = 1600, res = 300)
    print(GO_plot)
    dev.off()
    
    # Make dot plot
    dot_plot <- ggplot(combined_data, aes(x = Pct.of.bgd, y = reorder(Name, Pct.of.bgd), 
                                          color = P.value, size = Result.count)) +
      geom_point() + # Use geom_point for dot plots, now with size aesthetic
      scale_color_gradient(low = "blue", high = "red") + # Color gradient for P.values
      labs(title = paste(plot_title, "Genes GO dotplot", sep = " "), x = "Percentage of Background", y = "", color = "P value", size = "Counts") +
      theme_minimal() +
      theme(legend.position = "right",
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
            axis.title.x = element_text(size = 16, color = "black", face = "bold"),
            axis.title.y = element_text(size = 16, color = "black",face = "bold"),
            axis.text.x = element_text(size = 16, color = "black",face = "bold"),
            axis.text.y = element_text(size = 16, color = "black",face = "bold"),
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Make plot title bold and centered
            legend.title = element_text(face = "bold"), # Make legend title bold
            legend.text = element_text(face = "bold"))
    
    dot_plot_filename <- paste0("Dot_plot_", plot_title, ".png")
    png(file = dot_plot_filename, width = 2700, height = 1600, res = 300)
    print(dot_plot)
    dev.off()
