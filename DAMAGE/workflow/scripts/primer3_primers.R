#Analysis of the primer3 output and get the primers
library(purrr)
library(dplyr)
# Import from snakemake
output_dir <- snakemake@params[["output_dir"]]
output_path <- snakemake@input[["output_path"]]
primer3_primers <- snakemake@output[["primers_RDS"]]

output_list <- readLines(output_path)

sort_groups_primers <- function(output){
  
  output_file <- read.table(file.path(output))
  
  pairs_of_primer <- as.integer(gsub("PRIMER_PAIR_NUM_RETURNED=",
                                     "",
                                     output_file[grep("PRIMER_PAIR_NUM_RETURNED",output_file$V1),]))
  
  if(pairs_of_primer > 0){
    primer_list <- lapply((seq_len(pairs_of_primer) - 1),
                          function(primer_num){
                            list(
                              left_primer = gsub(paste0("PRIMER_LEFT_",
                                                           primer_num,
                                                           "_SEQUENCE="),
                                                    "",
                                                    output_file[grepl(paste0("PRIMER_LEFT_",
                                                                             primer_num,
                                                                             "_SEQUENCE="),
                                                                      output_file$V1),]),
                              left_primer_tm = gsub(paste0("PRIMER_LEFT_",
                                                              primer_num,
                                                              "_TM="),
                                                       "",
                                                       output_file[grepl(paste0("PRIMER_LEFT_",
                                                                                primer_num,
                                                                                "_TM="),
                                                                         output_file$V1),]),
                              left_primer_gc = gsub(paste0("PRIMER_LEFT_",
                                                              primer_num,
                                                              "_GC_PERCENT="),
                                                       "",
                                                       output_file[grepl(paste0("PRIMER_LEFT_",
                                                                                primer_num,
                                                                                "_GC_PERCENT="),
                                                                         output_file$V1),]),
                              right_primer = gsub(paste0("PRIMER_RIGHT_",
                                                           primer_num,
                                                           "_SEQUENCE="),
                                                    "",
                                                    output_file[grepl(paste0("PRIMER_RIGHT_",
                                                                             primer_num,
                                                                             "_SEQUENCE="),
                                                                      output_file$V1),]),
                              right_primer_tm = gsub(paste0("PRIMER_RIGHT_",
                                                              primer_num,
                                                              "_TM="),
                                                       "",
                                                       output_file[grepl(paste0("PRIMER_RIGHT_",
                                                                                primer_num,
                                                                                "_TM="),
                                                                         output_file$V1),]),
                              right_primer_gc = gsub(paste0("PRIMER_RIGHT_",
                                                              primer_num,
                                                              "_GC_PERCENT="),
                                                       "",
                                                       output_file[grepl(paste0("PRIMER_RIGHT_",
                                                                                primer_num,
                                                                                "_GC_PERCENT="),
                                                                         output_file$V1),]),
                              
                              product_size = gsub(paste0("PRIMER_PAIR_",
                                                         primer_num,
                                                         "_PRODUCT_SIZE="),
                                                  "",
                                                  output_file[grepl(paste0("PRIMER_PAIR_",
                                                                           primer_num,
                                                                           "_PRODUCT_SIZE="),
                                                                    output_file$V1),]),
                              product_tm = gsub(paste0("PRIMER_PAIR_",
                                                       primer_num,
                                                       "_PRODUCT_TM="),
                                                "",
                                                output_file[grepl(paste0("PRIMER_PAIR_",
                                                                         primer_num,
                                                                         "_PRODUCT_TM="),
                                                                  output_file$V1),])
                                
                            )
                          })
  }else{
    primer_list <- NULL
  }
  
  
  group_name <- if(grepl("pre_abs_output_",output)){
    gsub("pre_abs_output_","",basename(output)) 
  }else if(grepl("sanger_output_",output)){
    gsub("sanger_output_","",basename(output))
  }
  
  category <- if(grepl("pre_abs_output_",output)){
    "presence_absence"
  }else if(grepl("sanger_output_",output)){
    "sanger"
  }
  
  

  return(list(
    group = group_name,
    category = category,
    primers = primer_list
  ))
}

sorted_primers <- lapply(output_list, function(output) sort_groups_primers(output))

names(sorted_primers) <- sapply(output_list,
                                function(output){
                                  if(grepl("pre_abs_output_",output)){
                                    gsub("pre_abs_output_","",basename(output)) 
                                  }else if(grepl("sanger_output_",output)){
                                    gsub("sanger_output_","",basename(output))
                                  }
                                }
                                  )
  
saveRDS(sorted_primers,
        primer3_primers)
