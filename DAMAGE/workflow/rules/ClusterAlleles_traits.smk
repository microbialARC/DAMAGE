rule cluster_traits:
    input: 
        clusters_summary = os.path.join(config["input"], "thresher", "output", "clusters_summary.csv"),
        strains = os.path.join(config["input"], "thresher", "output", "plateau_strains.csv")
    output:
        cluster_traits = expand(os.path.join(config["output"], "cluster_alleles", "traits", "whatsgnu_traits_cluster{cluster}.csv"), cluster=clusters_list),
        traits_path = os.path.join(config["output"], "cluster_alleles", "traits", "cluster_traits_path.txt")
    params:
        traits_dir = os.path.join(config["output"], "cluster_alleles", "traits")
    script:
        os.path.join(BASE_PATH,"scripts/ClusterAlleles_traits.py")