rule allele_visual:
    conda:
        os.path.join(BASE_PATH,"envs/R_env.yaml")
    input:
        cluster_alleles_summary = os.path.join(config["output"], "cluster_alleles", "cluster_alleles_summary.RDS")
    output:
        ortholog_alleles_visual = os.path.join(config["output"], "cluster_alleles", "Cluster_ortholog_alleles.pdf")
    params:
        base_dir = os.path.join(config["output"], "cluster_alleles")
    script:
        os.path.join(BASE_PATH,"scripts","ClusterAlleles_visual.R")