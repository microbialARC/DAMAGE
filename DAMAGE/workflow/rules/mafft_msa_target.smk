rule mafft_msa_target:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        msa_list_path = os.path.join(config["output"], "mafft", "msa","msa_path.txt")
    output:
        msa_target = os.path.join(config["output"], "mafft", "msa","msa_target.RDS")
    params:
        min_length = config["min_product_length"]
    script:
        os.path.join(BASE_PATH, "scripts", "mafft_msa_target.R")