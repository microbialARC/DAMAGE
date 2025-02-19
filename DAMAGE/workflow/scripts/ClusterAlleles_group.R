# Libraries ----
library(dplyr)
library(parallel)
library(compiler)
# Import from snakemake
alleles_summary_path <- snakemake@input[["alleles_summary"]]
gene_presence_absence_path <- snakemake@input[["panaroo_gene_presence_absence"]]
clusters_summary_path <- snakemake@input[["clusters_summary_RDS"]]
# Load 
clusters_summary <- readRDS(alleles_summary_path)
gene_presence_absence <- read.csv(gene_presence_absence_path)
alleles_summary <- readRDS(alleles_summary_path)
# For now only ortholog will be used
alleles_summary <- alleles_summary[["ortholog"]]
all_alleles <- as.character(unlist(lapply(alleles_summary,
                                          function(cluster_alleles){
                                            cluster_alleles$Allele_ID
                                          })))


  
find_allele_groups <- function(all_alleles, alleles_summary, gene_presence_absence, clusters_summary) {
  redundant_allele_groups <- lapply(all_alleles,
                                    function(allele){
                                      # Find the annotation IDs for this allele
                                      ## First find the clusters
                                      cluster_entry <- which(sapply(alleles_summary,
                                                                    function(cluster_alleles) any(grepl(allele,cluster_alleles$Allele_ID))))
                                      ## Then find the position in the cluster 
                                      annotation_entry <- which(sapply(alleles_summary[[cluster_entry]]$Allele_ID,
                                                                       function(allele_position) allele == allele_position))
                                      ## Get the gene names
                                      annotation_gene <- unlist(alleles_summary[[cluster_entry]]$Gene[[annotation_entry]])
                                      ## Get the annotation description 
                                      annotation_description <- unlist(alleles_summary[[cluster_entry]]$Annotation[[annotation_entry]])
                                      ## Find the annotation IDs
                                      annotation_ids <- unlist(alleles_summary[[cluster_entry]]$AnnotationID[[annotation_entry]])
                                      
                                      # Use the annotation IDs to find the group determined by panaroo
                                      ## Find the row in gene_presence_absence
                                      row_entry <- which(apply(gene_presence_absence, 1, function(x) any(grep(annotation_ids[1], x, value = FALSE))))
                                      # These are all annotation IDs within the same group
                                      group_annotation_ids <- as.character(unlist(gene_presence_absence[row_entry,4:ncol(gene_presence_absence)]))
                                      #No empty annotation IDs
                                      group_annotation_ids <- group_annotation_ids[group_annotation_ids != ""]
                                      #No refound or psudo genes
                                      group_annotation_ids <- group_annotation_ids[!(grepl("refound",group_annotation_ids) |
                                                                                       grepl("pseudo",group_annotation_ids))]
                                      # Split the ; if any
                                      if(any(grepl("\\;",group_annotation_ids))){
                                        group_annotation_ids <- c(group_annotation_ids[-grep("\\;",group_annotation_ids)],
                                                                  unlist(strsplit(group_annotation_ids[grep("\\;",group_annotation_ids)],split = "\\;")))
                                      }
                                      # Exclude annotation ID of this allele from group_annotation_ids
                                      group_annotation_ids <- group_annotation_ids[!(group_annotation_ids %in% annotation_ids)]
                                      # Use the strain_df to tell what group_annotation_ids are from isolates in the clusters
                                      # Only those in the clusters will be analyzed
                                      cluster_annotation_ids <- group_annotation_ids[unlist(lapply(group_annotation_ids,function(id) gsub("_\\d{5}", "", id))) %in% unlist(sapply(clusters_summary,function(cluster) cluster$genomes))]
                                      non_cluster_annotation_ids <- group_annotation_ids[!(group_annotation_ids %in% cluster_annotation_ids)]
                                      
                                      #only when this isolate is in a cluster can the alleles possibly be in alleles_summary
                                      if(length(cluster_annotation_ids) > 0){
                                        # Use those annotation IDs to find the allele IDs in alleles_summary with Single-linkage clustering method
                                        homologous_alleles <- unique(unlist(sapply(cluster_annotation_ids,
                                                                                   function(annotation_id){
                                                                                     
                                                                                     isolate <-  gsub("_\\d{5}", "", annotation_id)
                                                                                     isolate_cluster_id <-  unlist(sapply(clusters_summary, \(x) if(isolate %in% x$genomes) x$cluster))
                                                                                     isolate_cluster_entry <- which(sapply(alleles_summary, function(cluster) isolate_cluster_id == cluster$Cluster))
                                                                                     # Note: 
                                                                                     # We are using the annotation_id to find the alleles in alleles_summary
                                                                                     # The reason why annotation_id can be absent in alleles_summary is that:
                                                                                     # alleles_summary only records the alleles unique to the clusters
                                                                                     # while annotation_id in the same panaroo group with other unique alleles might not be unique to their own cluster
                                                                                     # For example:
                                                                                     # We look at Allele_Cluster1_1, this allele is in alleles_summary because Allele_Cluster1_1 is 100% specificity and sensitivity to Cluster1
                                                                                     # and we have annotation_id Example_00000 in the same panaroo group with Allele_Cluster1_1, but Example_00000 is in Cluster 7.
                                                                                     # When we look into alleles from Cluster 7, we don't find it in alleles_summary, 
                                                                                     # This is because Example_00000 is not 100% specificity and sensitivity to Cluster 7.
                                                                                     
                                                                                     if(annotation_id %in% unlist(alleles_summary[[isolate_cluster_entry]]$AnnotationID)){
                                                                                       
                                                                                       homologous_allele_entry <- which(unlist(lapply(alleles_summary[[isolate_cluster_entry]]$AnnotationID, function(annotation_position) annotation_id %in% unlist(annotation_position))))
                                                                                       
                                                                                       alleles_summary[[isolate_cluster_entry]]$Allele_ID[homologous_allele_entry]
                                                                                       
                                                                                     }
                                                                                   })))
                                        
                                        
                                      }else{
                                        homologous_alleles <- NULL
                                      }  
                                      
                                      
                                      # This step might be redundant and can be integrated into previous step in the future
                                      # Aim is that after establishing cluster_alleles
                                      # We put their corresponding annotation ids as lists
                                      
                                      cluster_final_alleles <- sort(c(homologous_alleles,allele))
                                      
                                      cluster_final_alleles_annotation_ids <- lapply(cluster_final_alleles,
                                                                                     function(final_allele){
                                                                                       
                                                                                       final_cluster_entry <- which(sapply(alleles_summary,
                                                                                                                           function(cluster_alleles) any(grepl(final_allele,cluster_alleles$Allele_ID))))
                                                                                       
                                                                                       final_annotation_entry <- which(sapply(alleles_summary[[final_cluster_entry]]$Allele_ID,
                                                                                                                              function(allele_position) final_allele == allele_position))
                                                                                       
                                                                                       unlist(alleles_summary[[final_cluster_entry]]$AnnotationID[[final_annotation_entry]])
                                                                                       
                                                                                     })
                                      
                                      
                                      return(list(
                                        cluster_alleles = cluster_final_alleles,
                                        gene = annotation_gene,
                                        description = annotation_description,
                                        clusters = sort(unique(as.character(sapply(cluster_final_alleles, function(allele) strsplit(allele,split = "\\_")[[1]][3])))),
                                        cluster_annotation_ids = cluster_final_alleles_annotation_ids,
                                        non_cluster_annotation_ids = non_cluster_annotation_ids
                                      ))
                                    })
  return(redundant_allele_groups)
}

# Setup parallel processing
cl <- makeCluster(detectCores())

clusterExport(cl,
              c("alleles_summary",
                "gene_presence_absence",
                "clusters_summary",
                "find_allele_groups"))

clusterEvalQ(cl, {
  library(dplyr)
})


redundant_allele_groups <- parLapplyLB(cl, all_alleles, function(allele) {
  find_allele_groups(allele, 
                     alleles_summary, 
                     gene_presence_absence, 
                     clusters_summary)[[1]]
})

stopCluster(cl)
rm(cl)



sort_redundant_groups <- function(redundant_allele_groups){
  
  unique_entry <- !duplicated(lapply(redundant_allele_groups, function(group) paste(sort(group$cluster_alleles), collapse = ",")))
  
  unique_groups <- redundant_allele_groups[unique_entry]
  
  analyzed_alleles <- unlist(lapply(unique_groups,function(group) group$cluster_alleles))
  
  check_frequency <- as.data.frame(table(analyzed_alleles))
  
  if(all(check_frequency$Freq == 1)){
    # Add an ID to each group
    unique_groups <- lapply(seq_along(unique_groups), function(i) {
      unique_groups[[i]]$group_id <- i
      return(unique_groups[[i]])
    })
  }else{
    stop("Sort redundant groups failed")
  }
}

unique_allele_groups <- sort_redundant_groups(redundant_allele_groups)

saveRDS(unique_allele_groups,
        snakemake@output[["unique_alleles_groups"]])
