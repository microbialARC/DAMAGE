rule minimal_set:
    input:
        primer3_qc = os.path.join(config["output"], "primer3", "output", "primer3_qc.RDS"),
        unique_alleles_groups = os.path.join(config["output"], "cluster_alleles", "unique_alleles_groups.RDS"),
        primers_list = os.path.join(config["output"], "primer3", "primers_list.csv")
    output:
        minimal_primers_list = os.path.join(config["output"], "minimal_set", "minimal_primers.csv"),
        minimal_primer_sanke = os.path.join(config["output"], "minimal_set", "minimal_sanke.pdf"),
        minimal_primer_heatmap = os.path.join(config["output"], "minimal_set", "minimal_heatmap.pdf")
    params:
        minimal_set_dir = os.path.join(BASE_PATH, "minimal_set")
    script:
        os.path.join(BASE_PATH, "scripts", "minimal_set.R")