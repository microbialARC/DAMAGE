# DAMAGE
Diagnostic Allele Method for Actionable Genomic Epidemiology

---
## **Development Status**
Please be advised that this pipeline is still in its early development stage. It is subject to significant changes in terms of options, outputs, and other functionalities. Users should be prepared for potential modifications and updates in future releases.

## **Workflow**  
*Coming soon.*

---

## **Installation**  
DAMAGE utilizes a Snakemake pipeline complemented by a Python script for input validation, configuration file generation, and execution of the workflow.

### **Bioconda Installation**  
*Coming soon.*

### **Manual Installation via Git**  
1. Clone the repository:  
   ```bash
   git clone https://github.com/microbialARC/DAMAGE
   cd DAMAGE
2. Install dependencies:
    ```bash
    conda env create -f damage.yml
3. Activate DAMAGE environment and install:
    ```bash
    conda activate damage
    bash install.sh
## **Usage**
### Command line options
```
usage: DAMAGE [-h] -i INPUT -o OUTPUT [--mismatch MISMATCH] [-pmin MIN_PRODUCT_LENGTH] [-pmax MAX_PRODUCT_LENGTH] [-t THREADS]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the Thresher output directory.
  -o OUTPUT, --output OUTPUT
                        Path to the DAMAGE output directory.
  --mismatch MISMATCH   The maximum number of mismatch allowed in the PCR primers.
                        Default is 0
  -pmin MIN_PRODUCT_LENGTH, --min_product_length MIN_PRODUCT_LENGTH
                        The minimum length for PCR product (bp).
                        Default is 500
  -pmax MAX_PRODUCT_LENGTH, --max_product_length MAX_PRODUCT_LENGTH
                        The maximum length for PCR product (bp).
                        Default is 750
  -t THREADS, --threads THREADS
                        Number of threads to use.
                        Default is 8

```

**Note: DAMAGE will only be run when there is cluster identified in THRESHER**

## Outputs
### Diagnostic Alleles
* Visualization of alleles with 100% sensitivity and 100% specificity to each transmission cluster: `cluster_alleles/Cluster_ortholog_alleles.pdf`

### Primer Design (Primer3)
* List of primers for identified transmission clusters: `primer3/primers_list.csv`
    * **Columns:**
        * primer_id - Unique identifier for each primer pair
        * primer_group - Grouping of primers
        * left_primer - Forward primer sequence(5' to 3')
        * left_primer_gc - GC content of the forward primer
        * right_primer - Reverse primer sequence(5' to 3')
        * right_primer_gc - GC content of the reverse primer
        * product_size - Expected amplicon size in base pairs
        * primer_gene - Target gene for the primer pair
        * primer_gene_description - Functional description of the target gene
        * clusters - Target transmission clusters for the primer pair
* Visualization of primers across transmission clusters: `primer3/primers_genes_clusters.pdf` 

### Minimal Primer Set
* Minimal set of primers to detect all clusters if applicable: `minimal_set/minimal_primers.csv`
    * **Columns:**
        * step - Step in the Minimal Primer Set
        * primer_id - Unique identifier for each primer pair
        * primer_group - Grouping of primers
        * left_primer - Forward primer sequence(5' to 3')
        * left_primer_gc - GC content of the forward primer
        * right_primer - Reverse primer sequence(5' to 3')
        * right_primer_gc - GC content of the reverse primer
        * product_size - Expected amplicon size in base pairs
        * primer_gene - Target gene for the primer pair
        * primer_gene_description - Functional description of the target gene
        * clusters - Target transmission clusters for the primer pair
* Diagram of minimal primer set: `minimal_set/minimal_sanke.pdf`
* Heatmap visualization of minimal primer set coverage: `minimal_set/minimal_heatmap.pdf`


## **Author**
Qianxuan(Sean) She

[![PennMedicine](data/PennMedicine.png)](https://www.pennmedicine.org/)  
[![CHOP_Research](data/CHOP_Research.png)](https://www.research.chop.edu/)  
[![PennCHOP](data/PennCHOP.png)](https://www.research.chop.edu/pennchop-microbiome-program)
