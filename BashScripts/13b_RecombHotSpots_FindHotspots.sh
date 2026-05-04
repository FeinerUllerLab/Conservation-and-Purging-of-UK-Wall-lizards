#!/bin/bash

############################################## SET UP ##############################################

# Path to the input directory:
indir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step1

# Path to the output directory:
outdir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step2

# Suffix of the input files (assuming the input file names follow the structure
# <chromosome>_<suffix>.txt):
insuffix=w2kb_s1kb_f40kb

# Path to the genome file (a tab-delimited file with the structure <chromName><TAB><chromSize>):
genome_file=/home/feiner/Projects/UKwallies/scripts/Chr_Length_GenBank.txt

# Minimum fold difference between the mean recombination rate of the focal window and that of its
# flanking regions:
min_fold_diff_r=5

# Minimum number of SNPs in the focal window (set to 0 for no filtering):
min_snps_window=0

# Minimum number of SNPs in the flanking regions (set to 0 for no filtering):
min_snps_flanks=0

####################################################################################################

echo "Step 2: running (`date`)."

# Create the output directory if it does not already exist:
mkdir -p ${outdir}

# Create a suffix for the output files:
outsuffix=$(awk -v insuffix=${insuffix} -v min_fold_diff_r=${min_fold_diff_r} -v min_snps_window=${min_snps_window} -v min_snps_flanks=${min_snps_flanks} 'BEGIN {OFS = ""; print insuffix, "_fd", min_fold_diff_r, "_sw", min_snps_window, "_sf", min_snps_flanks}')

# Loop over the chromosomes:
while read -r chromosome chromosome_length; do

    echo "Processing chromosome ${chromosome} (${chromosome_length} bp)..."

    awk -v min_fold_diff_r=${min_fold_diff_r} -v min_snps_window=${min_snps_window} -v min_snps_flanks=${min_snps_flanks} '
    BEGIN {

        # Print the header:
        getline
        print $0

    }
    {
        # Get information about the current window:
        #chromosome = $1
        #window_start = $2
        #window_end = $3
        mean_r = $4
        n_snps = $5
        #upstream_flank_start = $6
        #upstream_flank_end = $7
        #downstream_flank_start = $8
        #downstream_flank_end = $9
        flanks_mean_r = $10
        flanks_n_snps = $11

        # If the current window qualifies as a hotspot, print the window:
        if (mean_r / flanks_mean_r >= min_fold_diff_r && n_snps >= min_snps_window && flanks_n_snps >= min_snps_flanks) {
            print $0
        }
    }
    ' ${indir}/${chromosome}_${insuffix}.txt > ${outdir}/${chromosome}_${outsuffix}.txt

done < ${genome_file}

echo "Step 2: done (`date`)."
