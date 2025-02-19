library(parallel)

generate_fasta <- function(group,
                           gene_data,
                           output_dir){
  
  # Extract nt sequence from gene_data
  ## Cluster alleles
  cluster_alleles_fasta <- unique(Biostrings::DNAStringSet(sapply(seq_len(length(unlist(group$cluster_annotation_ids))),
                                                           function(entry){
                                                             # Allele annotation ID: group$cluster_annotation_ids[[entry]][1]
                                                             # nt sequence of the allele
                                                             return(gene_data$dna_sequence[gene_data$annotation_id == unlist(group$cluster_annotation_ids)[entry]])
                                                           })))
  
  names(cluster_alleles_fasta) <- rep(group$cluster_alleles,length(cluster_alleles_fasta))
  
  ## Non-cluster alleles
  if(length(unlist(group$non_cluster_annotation_ids)) > 0){
    
    non_cluster_alleles_fasta <- unique(DNAStringSet(sapply(seq_len(length(group$non_cluster_annotation_ids)),
                                                            function(entry){
                                                              # Allele annotation ID: group$non_cluster_annotation_ids[entry]
                                                              # nt sequence of the allele
                                                              annotation_id <- paste(head(unlist(strsplit(group$non_cluster_annotation_ids[entry], "_")), 2), collapse = "_")
                                                              return(Biostrings::DNAString(gene_data$dna_sequence[gene_data$annotation_id == annotation_id]))
                                                            })))
    
    
    names(non_cluster_alleles_fasta) <- paste0("Ortholog_Allele_Non_Cluster_",
                                               seq_len(length(non_cluster_alleles_fasta)))
    
    combined_fasta <- c(cluster_alleles_fasta,
                        non_cluster_alleles_fasta)
  }else{
    combined_fasta <- cluster_alleles_fasta
  }
  
  Biostrings::writeXStringSet(combined_fasta,
                              filepath = file.path(output_dir,
                                                   paste0("Allele_Group",
                                                group$group_id,
                                                ".fasta")),
                              format = "fasta")
  
  return(combined_fasta)
}



# Import from snakemake

unique_allele_groups <- readRDS(snakemake@input[["unique_alleles_groups"]])
gene_data <- read.csv(snakemake@input[["panaroo_gene_data"]])
output_dir <- snakemake@params[["fasta_dir"]]
system(paste0("mkdir -p ",output_dir))

# Execute the function
cl <- makeCluster(detectCores())

clusterExport(cl,
              c("unique_allele_groups",
                "output_dir",
                "gene_data",
                "generate_fasta"))

clusterEvalQ(cl, {
  library(Biostrings)
})


combined_fasta <- parLapplyLB(cl,
                              unique_allele_groups,
                              function(group){
                                generate_fasta(group,
                                               gene_data,
                                               output_dir)
                              })

names(combined_fasta) <- sapply(unique_allele_groups, function(group) paste0("Allele_Group",group$group_id))

saveRDS(combined_fasta,
        snakemake@output[["alleles_fasta_RDS"]])

stopCluster(cl)
rm(cl)

