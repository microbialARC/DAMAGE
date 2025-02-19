rule ipcress_QC:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        ipcress_output = os.path.join(config["output"], "ipcress", "output", "ipcress_output_path.txt"),
        primer3_output = os.path.join(config["output"], "primer3", "output", "primer3_primers.RDS"),
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS")
    output:
        primer3_qc = os.path.join(config["output"], "primer3", "output", "primer3_qc.RDS"),
        primers_list = os.path.join(config["output"], "primer3", "primers_list.csv"),
        primers_visual = os.path.join(config["output"], "primer3", "primers_genes_clusters.pdf")
    params:
        output_dir = os.path.join(config["output"], "primer3", "output")
    script:
        os.path.join(BASE_PATH, "scripts", "ipcress_QC.R")