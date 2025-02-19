rule primer3_cmd:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        msa_target = os.path.join(config["output"], "mafft", "msa","msa_target.RDS"),
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS")
    output:
        cmd_path = os.path.join(config["output"], "primer3", "scripts", "primer3_cmd_path.txt"),
        primer3_input = os.path.join(config["output"], "primer3", "input", "primer3_input.RDS")
    params:
        input_dir = os.path.join(config["output"], "primer3", "input"),
        output_dir = os.path.join(config["output"], "primer3", "output"),
        scripts_dir = os.path.join(config["output"], "primer3", "scripts"),
        max_length = config["max_product_length"],
        min_length = config["min_product_length"]
    script:
        os.path.join(BASE_PATH, "scripts", "primer3_cmd.R")