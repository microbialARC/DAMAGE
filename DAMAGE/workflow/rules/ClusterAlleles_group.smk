rule allele_group:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        alleles_summary = os.path.join(config["output"], "cluster_alleles", "cluster_alleles_summary.RDS"),
        panaroo_gene_presence_absence = os.path.join(config["input"], "panaroo", "gene_presence_absence.csv"),
        clusters_summary_RDS = os.path.join(config["input"], "thresher", "output", "clusters_summary.RDS")
    params:
        scoring_dir = os.path.join(config["output"], "cluster_alleles", "scoring")
    output:
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS")
    script:
        os.path.join(BASE_PATH,"scripts","ClusterAlleles_group.R")