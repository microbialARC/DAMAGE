library(parallel)
library(compiler)
# Functions to summarize alleles
# This function will find the cluster-specific alleles(100% sensitivity and specificity) but the ortholog group shared by many cluster/strains ---- 
summarize_ortholog_alleles <- function(i) {
  
  cluster_id <- gsub("cluster|_alleles_scoring.csv",
                     "",
                     basename(scoring_list[i]))
  
  output <- read.csv(scoring_list[i])
  
  summary_allele_df_tmp <- list(Allele_ID = NA,
                                Cluster = cluster_id,
                                GNUscore = output$GNU.score,
                                Sensitivity = output$sensitivity,
                                Specificity = output$specificity,
                                Prot_sequence = output$allele,
                                DNA_sequence = list(),
                                Gene = list(),
                                Annotation = list(),
                                AnnotationID = list())
  
  #the id of the allele
  summary_allele_df_tmp$Allele_ID <- paste0("Ortholog_Allele_",
                                            "Cluster",
                                            cluster_id,
                                            "_",
                                            seq_len(nrow(output)))
  
  #iterate through rows in output 
  for(n in seq_len(nrow(output))){
    
    allele_list <- strsplit(output$Hits[n],
                            split = " ")[[1]]
    
    allele_list <- allele_list[allele_list %in% panaroo_gene_data$annotation_id]
    
    #DNA sequences of this allele
    summary_allele_df_tmp$DNA_sequence[[n]] <- unique(sapply(X = seq_len(length(allele_list)),
                                                             function(X){
                                                               panaroo_gene_data$dna_sequence[panaroo_gene_data$annotation_id == allele_list[X]]
                                                             }))
    
    #the gene name of the alleles
    summary_allele_df_tmp$Gene[[n]] <- unique(sapply(X = seq_len(length(allele_list)),
                                                     function(X){
                                                       panaroo_gene_data$gene_name[panaroo_gene_data$annotation_id == allele_list[X]]
                                                     }))
    
    #the gene annotation of the alleles
    summary_allele_df_tmp$Annotation[[n]] <- unique(sapply(X = seq_len(length(allele_list)),
                                                           function(X){
                                                             panaroo_gene_data$description[panaroo_gene_data$annotation_id == allele_list[X]]
                                                           }))
    
    #the gene annotation ID of the alleles in genomes
    summary_allele_df_tmp$AnnotationID[[n]] <- lapply(summary_allele_df_tmp$DNA_sequence[[n]], function(seq) {
      panaroo_gene_data$annotation_id[panaroo_gene_data$dna_sequence == seq]
    })
  }
  #return the list 
  return(summary_allele_df_tmp)
}


# load the input from snakemake ----
scoring_list <- list.files(path = snakemake@params[["scoring_dir"]],
                          pattern = "*alleles_scoring.csv",
                          all.files = TRUE,
                          full.names = TRUE,
                          recursive = FALSE)

panaroo_gene_data <- read.csv(snakemake@input[["panaroo_gene_data"]])
panaroo_gene_data <- read.csv(snakemake@input[["panaroo_gene_data"]])

#export the data needed for the parallel process
cl <- makeCluster(detectCores())
clusterExport(cl,
              c("panaroo_gene_data",
                "scoring_list"),
              envir = environment())

#Execute the parallel process of identifying the strains of different groups
cluster_ortholog_alleles_summary <- parLapplyLB(cl,
                                                seq_len(length(scoring_list)),
                                                summarize_ortholog_alleles)
cluster_alleles_summary <- list(
  ortholog = cluster_ortholog_alleles_summary
)

saveRDS(cluster_alleles_summary,
        snakemake@output[["alleles_summary"]])


stopCluster(cl)
rm(cl)
