rule primer3_primers:
    conda:
        os.path.join(BASE_PATH,"envs","R_env.yaml")
    input:
        output_path = os.path.join(config["output"], "primer3", "output","primer3_output_path.txt")
    output:
        primers_RDS = os.path.join(config["output"], "primer3", "output", "primer3_primers.RDS")
    params:
        output_dir = os.path.join(config["output"], "primer3", "output")
    script:
        os.path.join(BASE_PATH, "scripts", "primer3_primers.R")