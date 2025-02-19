# Libraries
library(Biostrings)
library(purrr)
library(parallel)
library(compiler)
# Import from snakemake
msa_list_path <- snakemake@input[["msa_list_path"]]
min_length <- snakemake@params[["min_length"]]
msa_target_path <- snakemake@output[["msa_target"]]
# Function to find the PCR target sequences ----
msa_list <- readLines(msa_list_path)

get_msa_target <- function(msa,
                           min_length){
  
  min_length <- as.numeric(min_length)
  msa_seq <- as.character(readDNAMultipleAlignment(msa))
  msa_length <- as.integer(nchar(msa_seq[1]))
  
  # Only when consensus sequences of genes with length >= 100 + 100 + min_length(bp) will be used for Sanger Sequencing
  if(msa_length >= (min_length + 200)){
    # Find the position and compositions of variants in msa_seq
    variant_pos <- lapply(seq_len(msa_length),
                          function(pos){
                            nt_pos <- msa_seq %>% purrr::map_chr(function(seq) substr(seq, pos, pos))
                            if(length(unique(nt_pos)) > 1){
                              return(list(
                                pos = pos,
                                nt_comp = nt_pos
                              ))
                            }else{
                              return(NULL)
                            }
                          }) %>% purrr::compact()
    
    # Find a range(min_length), which contains most SNVs,
    # and at least one SNV should be present in the range to differentiate all clusters involved
    # Such range is called a valid range
    
    # The sliding window starts from 100th nt position
    valid_range <- lapply(100:(msa_length-min_length-100), function(pos){
      range <- pos:(pos+min_length)
      variant_range <- variant_pos %>% purrr::keep(function(range_pos) range_pos$pos %in% range)
      if(length(variant_range) > 0){
        #To tell if the variants in this range can differentiate all clusters involved,
        #and if all cluster sequences different from non-cluster strain sequences
        msa_range <- msa_seq %>% purrr::map_chr(function(seq) substr(seq,pos,pos+min_length))
        cluster_msa_range <- msa_range[grepl("Ortholog_Allele_Cluster", names(msa_range))]
        non_clustere_msa_range <- msa_range[!grepl("Ortholog_Allele_Cluster", names(msa_range))]
        #length(unique(cluster_msa_range)) == length(cluster_msa_range)
        #all sequences are different, then proceed 
        if(length(unique(cluster_msa_range)) == length(cluster_msa_range)){
          if(!any(sapply(cluster_msa_range,function(seq) seq %in% non_clustere_msa_range))){
            return(
              list(start = pos,
                   end = pos + min_length,
                   variants = variant_range,
                   variants_num = length(variant_range))
            )
          }
        }
      }
    }) %>% purrr::compact()
    
    if(length(valid_range) > 0){
      
      # if we have multiple valid ranges, pick the one with the most number of variants
      max_variant_range_entry <- which.max(unlist(lapply(valid_range,function(range) range$variants_num)))
      max_variant_range <- valid_range[[max_variant_range_entry]]
      
      return(list(group = gsub("_msa.fasta","",basename(msa)),
                  msa = msa_seq,
                  msa_length = msa_length,
                  cluster_msa_consensus = consensusString(msa_seq[grepl("Ortholog_Allele_Cluster", names(msa_seq))],
                                                          ambiguityMap="N",
                                                          threshold=1),
                  target_seq = substr(consensusString(msa_seq[grepl("Ortholog_Allele_Cluster", names(msa_seq))],
                                                      ambiguityMap="N",
                                                      threshold=1),
                                      max_variant_range$start,
                                      max_variant_range$end),
                  target_seq_start =  max_variant_range$start,
                  target_seq_end = max_variant_range$end,
                  target_seq_variant = max_variant_range$variants,
                  target_seq_variant_num = max_variant_range$variants_num))
      
    }else{
      return(list(group = gsub("_msa.fasta","",basename(msa)),
                  msa = "N/A",
                  msa_length = msa_length,
                  cluster_msa_consensus = "N/A",
                  target_seq = "N/A",
                  target_seq_start = "N/A",
                  target_seq_end = "N/A",
                  target_seq_variant = "N/A",
                  target_seq_variant_num = "N/A"))
    }
    
  }else{
    return(list(group = gsub("_msa.fasta","",basename(msa)),
                msa = "N/A",
                msa_length = msa_length,
                cluster_msa_consensus = "N/A",
                target_seq = "N/A",
                target_seq_start = "N/A",
                target_seq_end = "N/A",
                target_seq_variant = "N/A",
                target_seq_variant_num = "N/A"))
  }
}

# Setup parallel processing ----
cl <- makeCluster(detectCores())

clusterExport(cl,
              c("msa_list",
                "min_length",
                "get_msa_target"))

clusterEvalQ(cl, {
  library(Biostrings)
  library(purrr)
})

msa_target <- parLapplyLB(cl,
                          msa_list,
                          function(msa) {
                            get_msa_target(msa,
                                           min_length)
                          })


saveRDS(msa_target,
        msa_target_path)

stopCluster(cl)
rm(cl)


