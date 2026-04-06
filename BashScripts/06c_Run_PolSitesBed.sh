#!/bin/bash -l

# Objective: Run the python script to get a bed with the polarised sites on a HCP system (Dardel)

# The resulting bed will contain only Ancestral-Reference SNPs for downstream analyses.

# modules 
module load python/3.12.3
module load PDCOLD/23.12
module load pysam/0.22.1-cpeGNU-23.12
 
python3 /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Conservation-and-Purging-of-Wall-lizards/BashScripts/14b_PolarisedSitesBED.py # Adjust path