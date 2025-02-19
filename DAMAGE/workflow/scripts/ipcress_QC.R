#Use ipcress to test the primers
library(dplyr)
library(purrr)
library(ggplot2)
# Import from snakemake
primer3_primers_path <- snakemake@input[["primer3_output"]]
ipcress_output_path <- snakemake@input[["ipcress_output"]]
unique_alleles_groups_path <- snakemake@input[["unique_alleles_groups"]]
primer3_primers_qc_path <- snakemake@output[["primer3_qc"]]
primer_list_path <- snakemake@output[["primers_list"]]
primers_visual <- snakemake@output[["primers_visual"]]
# Input ----
primer3_primers <- readRDS(primer3_primers_path)
unique_alleles_groups <- readRDS(unique_alleles_groups_path)
ipcress_output_list <- readLines(ipcress_output_path)
# Primers QC ----
primer3_primers <- primer3_primers[which(sapply(primer3_primers,function(primer) length(primer$primers)>0))]

primer_qc <- function(output_path){
  
  output_file <- readLines(output_path)
  
  output_file_experiment_unique <- sort(unique(gsub(" Experiment: ",
                                                    "",
                                                    output_file %>% 
                                                      purrr::keep(function(line) grepl("Experiment:",line)))))
  #Count the product within the genome using the primers
  disqualified_group_primer <- lapply(output_file_experiment_unique,
                                function(group_primer){
                                  product_count = length(output_file %>% 
                                                           purrr::keep(function(line) grepl(paste0(" Experiment: ",group_primer),line)))
                                  # If the product_count > 1 or = 0, the primer is disqualified
                                  # return to the output
                                  if(product_count > 1){
                                    return(group_primer)
                                  }
                                })
  
  return(disqualified_group_primer)
}

disqualified_primers <- unlist(purrr::compact(sapply(ipcress_output_list, primer_qc)))

# Remove disqualified_primers_output from sorted_primers
primer3_primers_qc <- primer3_primers[sapply(primer3_primers, function(primer_group) !(primer_group$group %in% disqualified_primers))]
primer3_primers_qc_path <- saveRDS(primer3_primers_qc,
                                   primer3_primers_qc_path)
# Export primers as spread sheet
primers_list_df <- do.call(rbind,
                      lapply(primer3_primers_qc,
                             function(primer_group){
                               if(primer_group$category == "sanger"){
                                 
                                 primer_entry <- which(sapply(unique_alleles_groups,
                                                              function(allele) 
                                                                allele$group_id == as.integer(gsub("Allele_Group",
                                                                                                   "",
                                                                                                   primer_group$group))))
                                 return(data.frame(
                                   primer_id = paste0(primer_group$group,"_",seq_along(primer_group$primers)),
                                   left_primer = sapply(primer_group$primers,function(primer) primer$left_primer),
                                   left_primer_gc = sapply(primer_group$primers,function(primer) primer$left_primer_gc),
                                   right_primer = sapply(primer_group$primers,function(primer) primer$right_primer),
                                   right_primer_gc = sapply(primer_group$primers,function(primer) primer$right_primer_gc),
                                   product_size = sapply(primer_group$primers,function(primer) primer$product_size),
                                   primer_gene = unique_alleles_groups[primer_entry][[1]]$gene,
                                   primer_gene_description = unique_alleles_groups[primer_entry][[1]]$description,
                                   clusters = paste(unique_alleles_groups[primer_entry][[1]]$clusters,
                                                    collapse = ", ")
                                 ))
                               }else{
                                 return(NULL)
                               }
                             }))

primers_list_df <- do.call(rbind,
                      lapply(primer3_primers_qc,
                             function(primer_group){
                               if(primer_group$category == "sanger"){
                                 
                                 primer_entry <- which(sapply(unique_alleles_groups,
                                                              function(allele) 
                                                                allele$group_id == as.integer(gsub("Allele_Group",
                                                                                                   "",
                                                                                                   primer_group$group))))
                                 return(data.frame(
                                   primer_id = paste0(primer_group$group,"_",seq_along(primer_group$primers)),
                                   primer_group = primer_group$group,
                                   left_primer = sapply(primer_group$primers,function(primer) primer$left_primer),
                                   left_primer_gc = sapply(primer_group$primers,function(primer) primer$left_primer_gc),
                                   right_primer = sapply(primer_group$primers,function(primer) primer$right_primer),
                                   right_primer_gc = sapply(primer_group$primers,function(primer) primer$right_primer_gc),
                                   product_size = sapply(primer_group$primers,function(primer) primer$product_size),
                                   primer_gene = unique_alleles_groups[primer_entry][[1]]$gene,
                                   primer_gene_description = unique_alleles_groups[primer_entry][[1]]$description,
                                   clusters = paste(unique_alleles_groups[primer_entry][[1]]$clusters,
                                                    collapse = ", ")
                                 ))
                               }else{
                                 return(NULL)
                               }
                             }))

# Export the primers passing QC
write.csv(primers_list_df,
          primer_list_path,
          quote = FALSE,
          row.names = FALSE)


# Visualize the number of primers for each transmission cluster ----
cluster_primer_visual <- function(primer3_primers_qc,unique_alleles_groups){
  # Total clusters in alleles_sum
  total_clusters <- sort(unique(unlist(sapply(unique_alleles_groups, function(group) group$clusters))))
  # Summarize the coverage of primers
  primer_clusters <- lapply(primer3_primers_qc,
                            function(primer_group){
                              if(primer_group$category == "sanger"){
                                
                                primer_entry <- which(sapply(unique_alleles_groups,
                                                             function(allele) 
                                                               allele$group_id == as.integer(gsub("Allele_Group",
                                                                                                  "",
                                                                                                  primer_group$group))))
                                
                                primer_clusters <- unique_alleles_groups[primer_entry][[1]]$clusters
                                
                                return(list(
                                  primer_group = primer_group$group,
                                  primer_clusters = primer_clusters
                                ))
                              }else{
                                return(
                                  NULL
                                )
                              }  
                            })
  # Find how many primers detect each of the cluster in total_clusters
  clusters_coverage <- lapply(total_clusters,
                              function(cluster){
                                cluster_primer_entry <- which(sapply(primer_clusters,
                                                                     function(primer) 
                                                                       cluster %in% primer$primer_clusters))
                                return(list(
                                  cluster = cluster,
                                  num_primers = length(cluster_primer_entry),
                                  primers = primer_clusters[cluster_primer_entry]
                                ))
                              })
  
  covered_clusters <- unique(unlist(sapply(primer_clusters,function(primer) primer$primer_clusters)))
  
  missed_clsuters <- setdiff(total_clusters,
                             covered_clusters)
  
  #generate the plots
  plot_df <- do.call(rbind,
                     sapply(clusters_coverage,
                            function(cluster){
                              
                              plot_df_tmp <- data.frame(
                                cluster = as.integer(gsub("Cluster","",cluster$cluster)),
                                num_primers = cluster$num_primers
                              )
                              
                              return(plot_df_tmp)
                            },
                            simplify = FALSE))
  
  plot_df$cluster <- factor(plot_df$cluster,
                            levels = as.integer(plot_df$cluster)[order(as.integer(plot_df$cluster))])
  
  plot <- ggplot(plot_df,
                 aes(x=cluster,
                     y=num_primers)) + 
    scale_y_continuous(expand = c(0, 0),
                       limits = function(x) c(0, max(x) * 1.05)) +
    geom_bar(fill = "#86C2E8",
             stat = "identity") + 
    geom_text(aes(label = num_primers), 
              color = "#305583",
              vjust = -0.5,
              size = 5,
              fontface = "bold") + 
    labs(x = "Cluster",
         y = "Genes") +
    theme(axis.text.x = element_text(face = "bold",
                                     hjust = 1,
                                     size = 15),
          axis.text.y = element_text(face = "bold",
                                     size = 25),
          legend.position = "top",
          axis.title.x = element_text(face = "bold", size = 30),
          axis.title.y = element_text(face = "bold", size = 30),
          plot.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.key.size = unit(1,"cm"),
          axis.line.x = element_line(linewidth = 1.5),
          axis.line.y = element_line(linewidth = 1.5),
          axis.ticks.length=unit(0.3,"cm"))
  
  pdf(file=primers_visual,
      width=22.5,
      height=8)
  print(plot)
  dev.off()
  
}

# Execute
cluster_primer_visual(primer3_primers_qc,
                      unique_alleles_groups)