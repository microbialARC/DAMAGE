rule allele_summary:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        cluster_scoring = expand(os.path.join(config["output"], "cluster_alleles", "scoring", "cluster{cluster}_alleles_scoring.csv"), cluster=clusters_list),
        panaroo_gene_data = os.path.join(config["input"], "panaroo", "gene_data.csv"),
        panaroo_gene_presence_absence = os.path.join(config["input"], "panaroo", "gene_presence_absence.csv"),
        clusters_summary = os.path.join(config["input"], "thresher", "output", "clusters_summary.csv")
    params:
        scoring_dir = os.path.join(config["output"], "cluster_alleles", "scoring")
    output:
        alleles_summary = os.path.join(config["output"], "cluster_alleles", "cluster_alleles_summary.RDS")
    script:
        os.path.join(BASE_PATH,"scripts","ClusterAlleles_summary.R")
