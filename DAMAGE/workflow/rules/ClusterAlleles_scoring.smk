rule allele_scoring:
    conda:
        os.path.join(BASE_PATH,"envs/allele_scoring.yaml")
    input:
        cluster_traits = expand(os.path.join(config["output"], "cluster_alleles", "traits", "whatsgnu_traits_cluster{cluster}.csv"), cluster=clusters_list)
    output:
        cluster_scoring = expand(os.path.join(config["output"], "cluster_alleles", "scoring", "cluster{cluster}_alleles_scoring.csv"), cluster=clusters_list)
    params:
        base_dir = os.path.join(config["output"], "cluster_alleles"),
        scoring_dir = os.path.join(config["output"], "cluster_alleles", "scoring"),
        report_dir = os.path.join(config["output"], "cluster_alleles", "reports"),
        whatsgnu_dir = os.path.join(config["input"], "whatsgnu"),
        scoring_script = os.path.join(BASE_PATH,"scripts","WhatsGNU_allele_scoring.py")
    shell:
        """
        # Import from snakemake
        scoring_script={params.scoring_script}
        whatsgnu_dir={params.whatsgnu_dir}
        report_dir={params.report_dir}
        scoring_dir={params.scoring_dir}
        base_dir={params.base_dir}
        cluster_traits=({input.cluster_traits})
        
        # Create report dir
        mkdir -p ${{report_dir}}

        # Copy reports from whatsgnu in THRESHER to DAMAGE pipeline
        cp -r ${{whatsgnu_dir}}/*/*_WhatsGNU_report.txt ${{report_dir}}

        # Create and go to scoring dir
        mkdir -p ${{scoring_dir}}
       
        # Run WhatsGNU allele scoring script through cluster traits files
        for trait_file in "${{cluster_traits[@]}}"; do
            cluster_id=$(basename "${{trait_file}}" .csv | sed 's#whatsgnu_traits_cluster##g')
            python ${{scoring_script}} ${{trait_file}} ${{report_dir}}
            awk -F',' 'NR==1 || ($10 == 100 && $11 == 100)' "${{base_dir}}/traits_allele_scoring.csv" > "${{base_dir}}/cluster${{cluster_id}}tmp"
            mv "${{base_dir}}/cluster${{cluster_id}}tmp" "${{scoring_dir}}/cluster${{cluster_id}}_alleles_scoring.csv"
            rm "${{base_dir}}/traits_allele_scoring.csv"
        done
        """