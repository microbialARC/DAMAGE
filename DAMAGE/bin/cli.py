#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )
    # Add the arguments
    # Argument for the input directory(output directory of Thresher)
    parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="""Path to the Thresher output directory."""
    )
    # Argument for DAMAGE output directory
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="""Path to the DAMAGE output directory."""
    )

    # Argument for the max mismatch allowed in the PCR primers
    parser.add_argument(
        "--mismatch",
        type=int,
        default=0,
        required=False,
        help="""The maximum number of mismatch allowed in the PCR primers.
Default is 0"""
    )

    # The minimum length for PCR product
    parser.add_argument(
        "-pmin",
        "--min_product_length",
        type=int,
        default=500,
        help="""The minimum length for PCR product (bp).
Default is 500"""
    )
    # The maximum length for PCR product
    parser.add_argument(
        "-pmax",
        "--max_product_length",
        type=int,
        default=750,
        help="""The maximum length for PCR product (bp).
Default is 750"""
    )


    # Argument for the number of threads
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="""Number of threads to use.
Default is 8"""
    )
    
    # Parse the arguments and check if there is any clusters identified by Thresher
    args = parser.parse_args()    
    if args.input:
        if not os.path.exists(args.input):
            raise SystemExit(
                "The directory does not exist."
            )
        clusters_summary = pd.read_csv(os.path.join(args.input, "thresher", "output", "clusters_summary.csv"))
        if clusters_summary.iloc[0, 0] == "No Clusters Found":
            raise SystemExit("No Clusters Found by Thresher, DAMAGE will NOT proceed")
    elif not args.input:
        raise SystemExit(
            "Please provide the Thresher output directory as the input"
        )
    return args

def create_config(args):

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    # Convert input and output path to absolute path
    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","workflow","config"), exist_ok=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","workflow","config", "config.yaml"), "w") as f:
        f.write("input: {}\n".format(input_path))
        f.write("output: {}\n".format(output_path))
        f.write("mismatch: {}\n".format(args.mismatch))
        f.write("threads: {}\n".format(args.threads))
        f.write("min_product_length: {}\n".format(args.min_product_length))
        f.write("max_product_length: {}\n".format(args.max_product_length))

# Main function
def main():
    args = parse_args()
    create_config(args)
    snakefile_abs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","workflow", "Snakefile")
    # Execute the workflow using the Snakefile
    # This needs to be fixed in the future
    # Because I haven't figured out how to use the Snakemake API for the newer version of Snakemake
    # https://snakemake-api.readthedocs.io/en/stable/api_reference/snakemake_api.html
    # So for now I will just use the os.system to execute the workflow
    os.system(f"snakemake --snakefile {snakefile_abs_path} --use-conda -f")

if __name__ == "__main__":
    sys.exit(main())