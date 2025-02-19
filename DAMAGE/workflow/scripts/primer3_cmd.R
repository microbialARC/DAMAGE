# Libraries
library(dplyr)
library(purrr)
library(Biostrings)
# Import from snakemake ----
msa_target_path <- snakemake@input[["msa_target"]]
max_length <- snakemake@params[["max_length"]]
min_length <- snakemake@params[["min_length"]]
unique_alleles_fasta_path <- snakemake@input[["unique_alleles_groups"]]
msa_target_path <- snakemake@input[["msa_target"]]
scripts_dir <- snakemake@params[["scripts_dir"]]
input_dir <- snakemake@params[["input_dir"]]
output_dir <- snakemake@params[["output_dir"]]
primer3_input_path <- snakemake@output[["primer3_input"]]
cmd_path <- snakemake@output[["cmd_path"]]

# Input ----
system(paste0("mkdir -p ",input_dir))
system(paste0("mkdir -p ",output_dir))
system(paste0("mkdir -p ",scripts_dir))
unique_alleles_fasta <- readRDS(unique_alleles_fasta_path)
msa_target <- readRDS(msa_target_path)

# Function for generating the input and commands for primer3 using for sanger sequencing PCR ----
generate_sanger_input <- function(group,
                                  min_length,
                                  max_length){
 if(group$target_seq != "N/A"){
      input <- c("PRIMER_TASK=generic",
                 "PRIMER_PICK_LEFT_PRIMER=1",
                 "PRIMER_PICK_RIGHT_PRIMER=1",
                 paste0("PRIMER_PRODUCT_SIZE_RANGE=",min_length,"-",max_length),
                 # sequence ID
                 paste0("SEQUENCE_ID=",
                        group$group),
                 # template sequence
                 paste0("SEQUENCE_TEMPLATE=",
                        gsub("\\-","N",group$cluster_msa_consensus)),
                 # target sequence
                 paste0("SEQUENCE_TARGET=",
                        group$target_seq_start,
                        ",",min_length),
                 "=")
      
      writeLines(input,
                 file.path(input_dir,
                           paste0("sanger_input_",
                                  group$group)),
                 sep = "\n")
      return(input)
    }
}
sanger_groups_input <- lapply(msa_target,
                               function(group)
                                 generate_sanger_input(group,
                                                       min_length,
                                                       max_length)) %>%
  compact()

# Function for generating the input and commands for primer3 using for gene presence/absence PCR ----
pre_abs_groups <- unique_alleles_fasta[which(sapply(unique_alleles_fasta,function(group) length(group) == 1))]
generate_pre_abs_input <- function(group,
                                   input_dir){
    input <- c("PRIMER_TASK=generic",
               "PRIMER_PICK_LEFT_PRIMER=1",
               "PRIMER_PICK_RIGHT_PRIMER=1",
               "PRIMER_PRODUCT_SIZE_RANGE=100-1000",
               # sequence ID
               paste0("SEQUENCE_ID=",
                      names(group)),
               # template sequence
               paste0("SEQUENCE_TEMPLATE=",
                      as.character(group)),
               "=")
    writeLines(input,
               file.path(input_dir,
                         paste0("pre_abs_input_",
                                names(group))),
               sep = "\n")
    return(input)
}

pre_abs_groups_input <- lapply(pre_abs_groups,
                               function(group)
                                 generate_pre_abs_input(group,
                                                        input_dir))


all_input <- list(
  sanger = sanger_groups_input,
  pre_abs = pre_abs_groups_input
)

saveRDS(all_input,
        primer3_input_path)

# Sanger commands ----
sanger_groups_cmd <- function(group,
                              input_dir,
                              output_dir,
                              scripts_dir){
  
  cmd <- paste0("primer3_core < ",
                file.path(input_dir,paste0("sanger_input_",group$group)),
                " --output=",file.path(output_dir,paste0("sanger_output_",group$group)))
  writeLines(cmd,
             file.path(scripts_dir,paste0("sanger_cmd_",group$group,".sh")))
}

invisible(lapply(msa_target,
                 function(group) 
                   if(group$target_seq != "N/A"){
                   sanger_groups_cmd(group,
                                     input_dir,
                                     output_dir,
                                     scripts_dir)}
                 ))


# Presence absence commands ----


pre_abs_groups_cmd <- function(group,
                               input_dir,
                               output_dir,
                               scripts_dir){
  
  cmd <- paste0("primer3_core < ",
                file.path(input_dir,paste0("pre_abs_input_",names(group))),
                " --output=",file.path(output_dir,paste0("pre_abs_output_",names(group))))
  writeLines(cmd,
             file.path(scripts_dir,paste0("pre_abs_cmd_",names(group),".sh")))
}

invisible(lapply(pre_abs_groups,
       function(group) pre_abs_groups_cmd(group,
                                          input_dir,
                                          output_dir,
                                          scripts_dir)))


# Export command path ----
all_cmd <- list.files(path = scripts_dir,
                      pattern = "*.sh",
                      all.files = TRUE,
                      full.names = TRUE,
                      recursive = FALSE)
writeLines(all_cmd,
           cmd_path)
