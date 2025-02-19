#Use ipcress to test the primers
library(dplyr)
library(purrr)
# Import fron snakemake
thresher_dir <- snakemake@params[["thresher_dir"]]
primers_path <- snakemake@input[["primers_RDS"]]
unique_alleles_groups_path <- snakemake@input[["unique_alleles_groups"]]
clusters_summary_path <- snakemake@input[["clusters_summary"]]
mismatch <- snakemake@params[["mismatch"]]
max_product_length <- snakemake@params[["max_product_length"]]
min_product_length <- snakemake@params[["min_product_length"]]
input_dir <- snakemake@params[["input_dir"]]
scripts_dir <- snakemake@params[["scripts_dir"]]
output_dir <- snakemake@params[["output_dir"]]
ipcress_cmd_path <- snakemake@output[["ipcress_cmd_path"]]
# Input ----
system(paste0("mkdir -p ",input_dir))
system(paste0("mkdir -p ",scripts_dir))
system(paste0("mkdir -p ",output_dir))

primers <- readRDS(primers_path)
unique_alleles_groups <- readRDS(unique_alleles_groups_path)
clusters_summary <- readRDS(clusters_summary_path)
# Only keep the groups with primers present from primer3
primers <- primers %>% 
  purrr::keep(function(primer) !is.null(primer$primers))

get_ipcress_cmd <- function(primer,
                            thresher_dir,
                            input_dir,
                            scripts_dir,
                            unique_alleles_groups,
                            clusters_summary,
                            min_product_length,
                            max_product_length,
                            mismatch){
  
  # The group ID of the primers
  group_id <- as.integer(gsub("Allele_Group",
                              "",
                              primer$group))
  
  # Input for ipcress
  input <- unlist(sapply(seq_along(primer$primers),
                         function(i){
                           primer_set <- primer$primers[[i]]
                           paste(paste0(primer$group,
                                        "_",
                                        i),
                                 primer_set$left_primer,
                                 primer_set$right_primer,
                                 min_product_length,
                                 max_product_length,
                                 collapse = " ")
                         }))
  
  writeLines(input,
             file.path(input_dir,
                       paste0(primer$group,
                              "_ipcress_input")))
  
  # The position of primers in allele_groups
  group_entry <- which(sapply(unique_alleles_groups,function(group) group$group_id == group_id))
  
  # IDs of the clusters detected by the primers
  cluster_id <- as.numeric(gsub("Cluster",
                                "",
                                unique_alleles_groups[[group_entry]]$clusters))
  
  # The isolates in the clusters
  clusters_entry <- which(sapply(clusters_summary,function(cluster) cluster$cluster == cluster_id))
  clusters_isolates <- unlist(sapply(clusters_entry,function(entry) clusters_summary[[entry]]$genomes))
  
  # Path to the genomes of the isolates
  clusters_isolates_path <- sapply(clusters_isolates,function(isolate) file.path(thresher_dir,"bakta_annotation",isolate,paste0(isolate,".fna")))
  
  # Command for ipcress
  
  group_ipcress_cmd <- unlist(sapply(clusters_isolates_path,
                                     function(isolate_path){
                                       paste("ipcress --mismatch",
                                             mismatch,
                                             "-i",
                                             file.path(input_dir,
                                                       paste0(primer$group,
                                                              "_ipcress_input")),
                                             "-s",
                                             isolate_path,
                                             ">",
                                             file.path(output_dir,
                                                       paste0(primer$group,
                                                              "_",
                                                              basename(isolate_path),
                                                              "_ipcress_output")),
                                             collapse = " ")
                                     }))
  
  writeLines(group_ipcress_cmd,
             file.path(scripts_dir,paste0("ipcress_scripts_",primer$group,".sh")))
  
  return(file.path(scripts_dir,paste0("ipcress_scripts_",primer$group,".sh")))
}

ipcress_cmd_all <- unlist(sapply(primers,
                             function(primer) get_ipcress_cmd(primer,
                                                              thresher_dir,
                                                              input_dir,
                                                              scripts_dir,
                                                              unique_alleles_groups,
                                                              clusters_summary,
                                                              min_product_length,
                                                              max_product_length,
                                                              mismatch)))

writeLines(ipcress_cmd_all,
           ipcress_cmd_path)
