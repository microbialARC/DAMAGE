rule allele_fasta:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS"),
        panaroo_gene_data = os.path.join(config["input"], "panaroo", "gene_data.csv")
    params:
        fasta_dir = os.path.join(config["output"], "cluster_alleles", "fasta")
    output:
        alleles_fasta_RDS = os.path.join(config["output"], "cluster_alleles", "fasta" ,"unique_alleles_fasta.RDS")
    script:
        os.path.join(BASE_PATH,"scripts","ClusterAlleles_fasta.R")