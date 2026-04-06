# Conservation-and-Purging-of-Wall-lizards
---

Scripts and data analysis used for the MSc project: "Shaped by Isolation: Genetic Load and Purging in Non-Native Populations of Wall Lizards (*Podarcis muralis*)."

This repository contains the scripts developed to analyze population structure, genetic diversity, inbreeding, and the purging of deleterious mutations using whole-genome data from non-native and native populations of wall lizards.

---

## List of Contents

### Bash Scripts

This directory includes the pipeline, starting from raw reads to the analysis of genomic data used in the study, specifically:

--- Data Preparation ---
* **00. Index the reference genome** for downstream analyses.
* **01. Trim adapters** and prepare raw reads for mapping.
* **02. Map reads** to the *Podarcis muralis* reference genome and generate mapping statistics.
* **03. Perform variant calling** on the BAM files after mapping.
* **04. Genotype and hard-filter variants.**
* **05. Organize and perform further filtering** on the VCFs for downstream use. The final step merge ingroups and outgroups in a single VCF used in polarisation.
* **06. Polarisation (uSFS)**  using est-sfs. The final results is a VCF with high-quality polarised SNPs. 

--- Analyses ---
* **Genetic structure:**
    * **07. Principal Component Analysis.**
    * **08. Genetic clustering analysis (ADMIXTURE).**
    * **09. Maximum Likelihood-tree.**
* **Genetic diversity and inbreeding:**
    * **10. Genetic diversity** (autosomal observed heterozygosity).
    * **11. Runs of homozygosity** (BCFtools and PLINK were used).
* **Purging of deleterious mutations:**
    * **12. Annotation of variants** based on predicted impact.
    * **13. Relative frequencies calculation** for derived allele ratios and block jackknifing.
* **Recombination sites analysis** 
    * **14. Recombination hotspots finding** 5kb windows with 5x recombinarion in three different scenarions (40Kb, 1Mb and chromosme backgrounds)
    * **15. Enrichement-depletion analysis** using a permutation test to evaluate deleterious allles in hotspots vs random coding regions in the genome. 


### R Scripts

This directory includes:

* **`Compiled_Results.R`**: This script compiles and visualizes all results from the previous analyses, including plots for genetic structure, diversity, inbreeding, and purging metrics.
* **`Sample_Locations.R`**: This script maps the geographic locations of the samples used in the project across their non-native (England) and native ranges (France and Italy), providing spatial context for the genetic data.

### Data

This directory contains the results from the bioinformatic analyses and the metadata used in the R scripts for visualization.

---