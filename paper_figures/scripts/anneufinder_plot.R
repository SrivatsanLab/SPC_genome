project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)


library(dplyr)
library(purrr)
library(tibble)
library(pbapply)  # for progress bar
library(ggplot2)
library(pheatmap)
setwd(file.path(project_root, "paper_figures/data/Anneufinder/"))

mat = read.csv("result.csv")

mat[1:5,1:5]
dist_mat <- dist(mat[, 5:ncol(mat)])
hc <- hclust(dist_mat, method = "complete")

mat_data <- as.matrix(mat[, 5:ncol(mat)])
mat_data_capped <- pmin(mat_data, 10)
colors <- c(
  "#f7fcf0",
  "#e0f3db",
  "#ccebc5",
  "#a8ddb5",
  "#7bccc4",
  "#4eb3d3",
  "#2b8cbe",
  "#0868ac",
  "#084081",
  "#041C60"
)

# Create heatmap with clustering
ph = pheatmap(mat_data_capped,
         cluster_rows = FALSE,
         color = colors,
         breaks = seq(0, 10, length.out = 11),  
         cluster_cols = TRUE,show_colnames = F,show_rownames = F,
         scale = "none",  # change to "row" or "column" if you want z-score normalization
         fontsize_row = 6,
         fontsize_col = 6)

col_clusters <- cutree(ph$tree_col, k = 10)
table(col_clusters)

mean(apply(mat_data[,col_clusters == 1],MARGIN = 2,FUN = mean))
mean(apply(mat_data[,col_clusters == 2],MARGIN = 2,FUN = mean))
mean(apply(mat_data[,col_clusters == 3],MARGIN = 2,FUN = mean))
mean(apply(mat_data[,col_clusters == 4],MARGIN = 2,FUN = mean))


annotation_col <- data.frame(Cluster = factor(col_clusters))

rownames(annotation_col) <- colnames(mat_data_capped)

pheatmap(mat_data_capped,
         cluster_rows = FALSE,
         color = colors,
         breaks = seq(0, 10, length.out = 11),  
         cluster_cols = TRUE,show_colnames = F,show_rownames = F,
         scale = "none",  # change to "row" or "column" if you want z-score normalization
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_col = annotation_col)
