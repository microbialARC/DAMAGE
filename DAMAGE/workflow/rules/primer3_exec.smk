rule primer3_exec:
    conda:
        os.path.join(BASE_PATH,"envs","primer3.yaml")
    input:
        cmd_path = os.path.join(config["output"], "primer3", "scripts", "primer3_cmd_path.txt")
    output:
        output_path = os.path.join(config["output"], "primer3", "output","primer3_output_path.txt")
    threads:
        config["threads"]
    params:
        output_dir = os.path.join(config["output"], "primer3", "output")
    shell:
        """
        sed -i 's#//#/#g' {input.cmd_path}
        module load parallel
        parallel --jobs {threads} bash :::: {input.cmd_path}
        ls {params.output_dir}/*output* > {output.output_path}
        """