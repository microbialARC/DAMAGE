rule mafft_msa_exec:
    conda:
        os.path.join(BASE_PATH,"envs","mafft.yaml")
    input:
        cmd_path = os.path.join(config["output"], "mafft", "scripts", "mafft_cmd_path.txt")
    output:
        msa_path = os.path.join(config["output"], "mafft", "msa","msa_path.txt")
    threads:
        config["threads"]
    params:
        fasta_dir = os.path.join(config["output"], "cluster_alleles", "fasta"),
        msa_dir = os.path.join(config["output"], "mafft", "msa")
    shell:
        """
        sed -i 's#//#/#g' {input.cmd_path}
        module load parallel
        parallel --jobs {threads} bash :::: {input.cmd_path}

        ls {params.msa_dir}/*_msa.fasta > {params.msa_dir}/msa_output.txt
        msa_output=$(cat {params.msa_dir}/msa_output.txt | xargs -n 1 basename | sed 's#_msa.fasta##g')
        cmd_output=$(cat {input.cmd_path} | xargs -n 1 basename | sed 's#_mafft.sh##g')

        if [ "$(echo ${{msa_output}} | tr ' ' '\n' | sort)" == "$(echo ${{cmd_output}} | tr ' ' '\n' | sort)" ]; then
            ls {params.msa_dir}/*_msa.fasta > {output.msa_path}
            rm {params.msa_dir}/msa_output.txt
        fi
        """