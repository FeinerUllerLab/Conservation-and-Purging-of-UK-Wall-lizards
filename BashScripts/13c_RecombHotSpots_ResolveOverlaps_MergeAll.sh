#!/bin/bash

############################################## SET UP ##############################################

# Path to the input directory:
indir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step2

# Path to the output directory:
outdir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step3b

# Suffix of the input files (assuming the input file names follow the structure
# <chromosome>_<suffix>.txt):
insuffix=w2kb_s1kb_f40kb_fd5_sw0_sf0

# Path to the genome file (a tab-delimited file with the structure <chromName><TAB><chromSize>):
genome_file=/home/feiner/Projects/UKwallies/scripts/Chr_Length_GenBank.txt

####################################################################################################

echo "Step 3b: running (`date`)."

# Create the output directory if it does not already exist:
mkdir -p ${outdir}

# Create a suffix for the output files:
outsuffix=$(awk -v insuffix=${insuffix} 'BEGIN {OFS = ""; print insuffix, "_overlaps-mergeall"}')

# Loop over the chromosomes:
while read -r chromosome chromosome_length; do

    echo "Processing chromosome ${chromosome} (${chromosome_length} bp)..."

    awk -v chromosome=${chromosome} '
    BEGIN {

        # Set the output field separator to tab:
        OFS = "\t"
        
        # Skip the header of the input file:
        getline

        # Print a header for the output file:
        print "chr", "start", "end"

        # Set initial values for the start and end positions of the current block of windows:
        block_start = 0
        block_end = 0

    }
    {
        # Get information about the current window:
        #chromosome = $1
        window_start = $2
        window_end = $3
        #mean_r = $4
        #n_snps = $5
        #upstream_flank_start = $6
        #upstream_flank_end = $7
        #downstream_flank_start = $8
        #downstream_flank_end = $9
        #flanks_mean_r = $10
        #flanks_n_snps = $11

        # If the current window overlaps with or is adjacent to the current block of windows:
        if (window_start <= block_end) {

            # Update the end position of the current block of windows:
            block_end = window_end

        }
        # If the current window does not overlap with or is not adjacent to the current block of
        # windows:
        else {

            # Print the current block of windows:
            if (block_end != block_start) {
                print chromosome, block_start, block_end
            }

            # Reset the current block of windows:
            block_start = window_start
            block_end = window_end

        }
    }
    END {

        # Print the last block of windows:
        if (block_end != block_start) {
            print chromosome, block_start, block_end
        }

    }
    ' ${indir}/${chromosome}_${insuffix}.txt > ${outdir}/${chromosome}_${outsuffix}.txt

done < ${genome_file}

echo "Step 3b: done (`date`)."
