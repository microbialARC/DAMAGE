rule icpress_cmd:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        primers_RDS = os.path.join(config["output"], "primer3", "output", "primer3_primers.RDS"),
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS"),
        clusters_summary = os.path.join(config["input"], "thresher", "output", "clusters_summary.RDS")
    output:
        ipcress_cmd_path = os.path.join(config["output"], "ipcress", "scripts", "ipcress_cmd_path.txt")
    params:
        thresher_dir = config["input"],
        mismatch = config["mismatch"],
        max_product_length = config["max_product_length"],
        min_product_length = config["min_product_length"],
        input_dir = os.path.join(config["output"], "ipcress", "input"),
        scripts_dir = os.path.join(config["output"], "ipcress", "scripts"),
        output_dir = os.path.join(config["output"], "ipcress", "output")
    script:
        os.path.join(BASE_PATH, "scripts", "ipcress_cmd.R")