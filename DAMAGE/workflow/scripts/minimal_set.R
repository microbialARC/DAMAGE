# Find the smallest set of diagnostic genes that could identify every cluster.
library(dplyr)
library(purrr)
# Import from Snakemake
minimal_set_dir <- snakemake@params[["minimal_set_dir"]]
primer_groups_path <- snakemake@input[["primer3_qc"]]
primers_list_path <- snakemake@input[["primers_list"]]
unique_alleles_groups_path <- snakemake@input[["unique_alleles_groups"]]
minimal_primers_list_path <- snakemake@output[["minimal_primers_list"]]
minimal_primer_sanke_path <- snakemake@output[["minimal_primer_sanke"]]
minimal_primer_heatmap_path <- snakemake@output[["minimal_primer_heatmap"]]
# Input 
system(paste0("mkdir -p ",minimal_set_dir))
primer_groups <- readRDS(primer_groups_path)
unique_alleles_groups <- readRDS(unique_alleles_groups_path)
primers_list <- read.csv(primers_list_path)
# Get matrix for lpsolve ----
get_alleles_cluster_matrix <- function(primer_groups,
                                       unique_alleles_groups){
  
  primer_clusters <- lapply(primer_groups,
                            function(group){
                              
                              group_entry <- which(sapply(unique_alleles_groups,function(allele) allele$group_id == as.integer(gsub("Allele_Group","",group$group))))
                              group_clusters <- unique_alleles_groups[[group_entry]]$clusters
                              return(list(
                                primer_group = group$group,
                                primer_clusters = group_clusters
                              ))
                            })
  
  covered_clusters <- unique(unlist(sapply(primer_clusters,function(primer) primer$primer_clusters)))
  
  do.call(rbind,
          sapply(primer_groups,
                 function(group){
                   group_entry <- which(sapply(unique_alleles_groups,function(allele) allele$group_id == as.integer(gsub("Allele_Group","",group$group))))
                   group_clusters <- unique_alleles_groups[[group_entry]]$clusters
                   
                   alleles_cluster_matrix_tmp <- as.data.frame(matrix(nrow = 1,
                                                                      ncol = length(covered_clusters)))
                   
                   colnames(alleles_cluster_matrix_tmp) <- covered_clusters
                   #fill 1 to the clusters covered by this allele of the allele
                   alleles_cluster_matrix_tmp[,colnames(alleles_cluster_matrix_tmp) %in% group_clusters] <- 1
                   #fill 0 to the rest of the clusters
                   alleles_cluster_matrix_tmp[is.na(alleles_cluster_matrix_tmp[1,])] <- 0
                   #rowname 
                   rownames(alleles_cluster_matrix_tmp) <- group$group
                   return(alleles_cluster_matrix_tmp)
                 },
                 simplify = FALSE))
  
  
}

alleles_cluster_matrix <- get_alleles_cluster_matrix(primer_groups,
                                                     unique_alleles_groups)

# Find the smallest set of alleles identifying all clusters ----
library(lpSolve)

get_minimal_set_primers <- function(alleles_cluster_matrix){
  
  f.obj <- rep(1, nrow(alleles_cluster_matrix))
  f.con <- t(alleles_cluster_matrix)
  f.dir <- rep(">=", ncol(alleles_cluster_matrix))
  f.rhs <- rep(1, ncol(alleles_cluster_matrix))
  lp_result <- lp("min", f.obj, f.con, f.dir, f.rhs, all.bin = TRUE)
  # Extract the selected alleles
  selected_alleles <- which(lp_result$solution == 1)
  selected_alleles <- rownames(alleles_cluster_matrix[selected_alleles,])
  
  selected_alleles_clusters_count <- sapply(selected_alleles,
                                            function(allele){
                                              sum(alleles_cluster_matrix[rownames(alleles_cluster_matrix) == allele,] == 1)
                                            })
  # Reorder
  selected_alleles_clusters_count <- selected_alleles_clusters_count[order(-selected_alleles_clusters_count)]
  return(selected_alleles_clusters_count)
}

selected_alleles <- get_minimal_set_primers(alleles_cluster_matrix)

# Export the minimal_set
minimal_set_primers <- do.call(rbind,
                               lapply(seq_along(selected_alleles),function(index){
                                 cbind(step = index,
                                       primers_list[primers_list$primer_group == names(selected_alleles[index]),])
                               }))

write.csv(minimal_set_primers,
          minimal_primers_list_path,
          quote = FALSE,
          row.names = FALSE)
# Make sanke plot ----
## Make the data frame asinput for sanke plot

make_sanke_plot_df <- function(unique_alleles_groups,
                               selected_alleles,
                               alleles_cluster_matrix){
  
  # a in-function variable recording what's left in the clusters
  clusters_left <- colnames(alleles_cluster_matrix)
  # data frame as input for sanke plot
  plot_df <- do.call(rbind,
                     sapply(seq_along(selected_alleles),
                            function(allele_pos){
                              
                              alllele_group_id <- as.integer(gsub("Allele_Group","",names(selected_alleles[allele_pos])))
                              
                              allele_entry <- which(unlist(sapply(unique_alleles_groups,function(group) group$group_id == alllele_group_id)))
                              
                              clusters_identifed_by_allele <- unique_alleles_groups[[allele_entry]]$clusters
                              
                              plot_df_tmp <- data.frame(step = allele_pos,
                                                        cluster = clusters_left,
                                                        result = ifelse(clusters_left %in% clusters_identifed_by_allele,
                                                                        "Identified by this allele",
                                                                        "Unidentified"))
                              clusters_left <<- clusters_left[!(clusters_left %in% clusters_identifed_by_allele)]
                              return(plot_df_tmp)
                            },
                            simplify = FALSE))
  
  plot_df <- rbind(data.frame(step = "Start",
                              cluster = colnames(alleles_cluster_matrix),
                              result = "Unidentified"),
                   plot_df)
  
  plot_df$cluster <- factor(plot_df$cluster,
                            levels = paste0("Cluster",as.integer(gsub("Cluster","",unique(plot_df$cluster)))[order(as.integer(gsub("Cluster","",unique(plot_df$cluster))))]))
  
  plot_df$step <- factor(plot_df$step,
                         levels = c("Start",
                                    seq_along(selected_alleles)))
  
  return(plot_df)
}

plot_df <- make_sanke_plot_df(unique_alleles_groups,
                              selected_alleles,
                              alleles_cluster_matrix)

## Sanke plot 
library(ggplot2)
library(ggalluvial)
library(Polychrome)

make_sanke_plot <- function(plot_df){
  # Allele result color:
  result_colors <- c("#B3C16D",
                     "#E3E0DC")
  names(result_colors) <- c("Identified by this allele",
                            "Unidentified")
  #Cluster colors: 
  cluster_colors <- Polychrome::createPalette(length(levels(plot_df$cluster)),
                                              c("#ff0000",
                                                "#00ff00",
                                                "#0000ff"))
  
  names(cluster_colors) <- levels(plot_df$cluster)
  
  plot <- ggplot(plot_df,
                 aes(x = step,
                     stratum = result,
                     alluvium = cluster)) +
    geom_stratum(aes(fill = result),
                 color = "#00000000") +
    geom_flow(stat = "alluvium",
              lode.guidance = "frontback",
              aes(color = "transparent",
                  fill = cluster)) +
    scale_fill_manual(values = c(cluster_colors, result_colors)) +
    scale_color_manual(values = result_colors) +
    guides(fill=guide_legend(title="Cluster")) + 
    scale_size_identity() +
    labs(x = "Gene",
         y = "Remaining Clusters") +
    scale_y_continuous(expand = c(0, 0),
                       limits = function(x) c(0, max(x) * 1.05)) +
    theme(axis.text.x = element_text(face = "bold",
                                     angle = 0,
                                     size = 27.5),
          axis.text.y = element_text(face = "bold",
                                     size = 25),
          axis.title.x = element_text(face = "bold", size = 30),
          axis.title.y = element_text(face = "bold", size = 30),
          plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 12.5),
          legend.title = element_text(size = 15),
          axis.line.x = element_line(linewidth = 1.5),
          axis.line.y = element_line(linewidth = 1.5),
          axis.ticks.length=unit(0.3,"cm"))
  
  pdf(file=minimal_primer_sanke_path,
      width=27.5,
      height=12.5)
  print(plot)
  dev.off()
}

make_sanke_plot(plot_df)

# Make heatmap to show the coverage ----
make_heatmap <- function(selected_alleles,
                         alleles_cluster_matrix){
  
  # Make the data frame as input for heatmap
  heatmap_df <- do.call(rbind,
                        sapply(seq_along(selected_alleles),
                               function(allele_pos){
                                 allele <- names(selected_alleles)[allele_pos]
                                 allele_matrix <- alleles_cluster_matrix[rownames(alleles_cluster_matrix) == allele,]
                                 true_clusters <- colnames(allele_matrix[which(allele_matrix == 1)])
                                 return(data.frame(
                                   allele = allele,
                                   gene_id = allele_pos,
                                   clusters = gsub("Cluster","",colnames(alleles_cluster_matrix)),
                                   result = ifelse(colnames(alleles_cluster_matrix) %in% true_clusters,
                                                   "Identified",
                                                   "Not identified")
                                 ))
                               },
                               simplify = FALSE))
  
  heatmap_df$gene_id <- factor(heatmap_df$gene_id,
                               levels = rev(sort(unique(heatmap_df$gene_id))))
  
  heatmap_df$clusters <- factor(heatmap_df$clusters,
                                levels = sort(unique(as.integer(heatmap_df$clusters))))
  
  heatmap_df$result <- factor(heatmap_df$result,
                              levels = c("Identified",
                                         "Not identified"))
  # Make the heatmap
  result_color <- c("#B3C16D",
                    "#E3E0DC")
  names(result_color) <- c("Identified",
                           "Not identified")
  
  heatmap_plot <- ggplot(heatmap_df,
                         aes(x = clusters,
                             y = gene_id,
                             fill = result)) +
    geom_tile(color = "black",
              size = 0.5) +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = result_color,
                      guide = "none") + 
    guides(fill=guide_legend(title="Result")) + 
    labs(x = "Clusters",
         y = "Gene") + 
    theme(legend.position = "top",
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(face = "bold", size = 15),
          axis.title.y = element_text(face = "bold", size = 15))
  
  pdf(file=minimal_primer_heatmap_path,
      width=12.5,
      height=12.5)
  print(heatmap_plot)
  dev.off()
}

make_heatmap(selected_alleles,
             alleles_cluster_matrix)
