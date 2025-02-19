rule mafft_msa_cmd:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        unique_alleles_fasta_path = os.path.join(config["output"], "cluster_alleles", "fasta", "unique_alleles_fasta.RDS"),
        unique_alleles_groups_path = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS")
    output:
        cmd_path = os.path.join(config["output"], "mafft", "scripts", "mafft_cmd_path.txt")
    params:
        fasta_dir = os.path.join(config["output"], "cluster_alleles", "fasta"),
        scripts_dir = os.path.join(config["output"], "mafft", "scripts"),
        msa_dir = os.path.join(config["output"], "mafft", "msa")
    script:
        os.path.join(BASE_PATH, "scripts", "mafft_msa_cmd.R")