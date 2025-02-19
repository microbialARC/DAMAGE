import glob
import pandas as pd
import os

# Import from snakemake
clusters_summary = pd.read_csv(snakemake.input["clusters_summary"])
strains = pd.read_csv(snakemake.input["strains"])
traits_dir = snakemake.params["traits_dir"]
os.makedirs(traits_dir, exist_ok=True)

# Loop through clusters and create traits files
clusters_list = clusters_summary["cluster"].tolist()
for cluster in clusters_list:
    cluster_genomes =  clusters_summary[clusters_summary["cluster"] == cluster]["genomes"].values[0].split("|")
    non_cluster_genomes = strains[~strains["genome"].isin(cluster_genomes)]["genome"].tolist()
    # Create traits dataframe
    df = pd.DataFrame({"genome": cluster_genomes + non_cluster_genomes,
                       "traits": [1] * len(cluster_genomes) + [0] * len(non_cluster_genomes)})
    # Save traits dataframe
    df.to_csv(os.path.join(traits_dir,
                           f"whatsgnu_traits_cluster{cluster}.csv"),
                           sep=",",
                           index=False)

# Check if all traits files are created, and list the path of the traits files
all_traits_files = glob.glob(os.path.join(traits_dir, "whatsgnu_traits_cluster*.csv"))
traits_files_count = list()
for traits_file in all_traits_files:
    traits_files_count.append(os.path.basename(traits_file).replace("whatsgnu_traits_cluster", "").replace(".csv", ""))
    traits_files_count = list(map(int, traits_files_count))

if all(cluster in traits_files_count for cluster in clusters_list):
    with open(os.path.join(traits_dir, "cluster_traits_path.txt"), "w") as f:
        for cluster in clusters_list:
            f.write(os.path.join(traits_dir, f"whatsgnu_traits_cluster{cluster}.csv\n"))
