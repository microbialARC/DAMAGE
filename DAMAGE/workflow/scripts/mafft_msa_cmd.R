# Generate the mafft commands for alleles
# Libraries
library(Biostrings)
# Import from snakemake
scripts_dir <- snakemake@params[["scripts_dir"]]
msa_dir <- snakemake@params[["msa_dir"]]
fasta_dir <- snakemake@params[["fasta_dir"]]
cmd_path <- snakemake@output[["cmd_path"]]
unique_alleles_fasta_path <- snakemake@input[["unique_alleles_fasta_path"]]
# Get the mafft commands 
unique_alleles_fasta <- readRDS(unique_alleles_fasta_path)

get_mafft_cmd <- function(unique_alleles_fasta,
                          scripts_dir,
                          msa_dir,
                          fasta_dir){
  
  system(paste0("mkdir -p ",scripts_dir))
  system(paste0("mkdir -p ",msa_dir))
  
  for(i in seq_along(unique_alleles_fasta)){
    if(length(unique_alleles_fasta[[i]]) > 1){
      
      cmd <- paste0("mafft --auto --thread 1 ",
                    file.path(fasta_dir,paste0(names(unique_alleles_fasta[i]),".fasta")),
                    " > ",
                    file.path(msa_dir,paste0(names(unique_alleles_fasta[i]),"_msa.fasta")))
      
      writeLines(cmd,
                 file.path(scripts_dir,paste0(names(unique_alleles_fasta[i]),"_mafft.sh")))
    }
  }
}


get_mafft_cmd(unique_alleles_fasta,
              scripts_dir,
              msa_dir,
              fasta_dir)


all_cmd <- list.files(path = scripts_dir,
                      pattern = "*_mafft.sh",
                      all.files = TRUE,
                      full.names = TRUE,
                      recursive = FALSE)

writeLines(all_cmd,
           cmd_path)
