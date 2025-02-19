rule ipcress_exec:
    conda:
        os.path.join(BASE_PATH,"envs","exonerate.yaml")
    input:
        ipcress_cmd_path = os.path.join(config["output"], "ipcress", "scripts", "ipcress_cmd_path.txt")
    output:
        ipcress_output_path = os.path.join(config["output"], "ipcress", "output","ipcress_output_path.txt")
    threads:
        config["threads"]
    params:
        ipcress_output_dir = os.path.join(config["output"], "ipcress", "output")
    shell:
        """
        mkdir -p {params.ipcress_output_dir}
        sed -i 's#//#/#g' {input.ipcress_cmd_path}
        module load parallel
        parallel --jobs {threads} bash :::: {input.ipcress_cmd_path}
        ls {params.ipcress_output_dir}/*_ipcress_output > {params.ipcress_output_dir}/ipcress_output_path.txt
        """