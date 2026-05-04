#!/bin/bash

############################################## SET UP ##############################################

# Path to the input directory:
indir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step3b

# Path to the output directory:
outdir=/home/feiner/Projects/UKwallies/RecombHotSpots/Step4

# Suffix of the input files (assuming the input file names follow the structure
# <chromosome>_<suffix>.txt):
insuffix=w2kb_s1kb_f40kb_fd5_sw0_sf0_overlaps-mergeall

# Path to the genome file (a tab-delimited file with the structure <chromName><TAB><chromSize>):
genome_file=/home/feiner/Projects/UKwallies/scripts/Chr_Length_GenBank.txt

####################################################################################################

echo "Step 4: running (`date`)."

# Create the output directory if it does not already exist:
mkdir -p ${outdir}

# Create a suffix for the output file:
outsuffix=${insuffix}

# Remove the output file if it already exists (otherwise the output of the script will be appended
# to the current content of the output file, rather than overwriting it):
if [ -f ${outdir}/allchrom_${outsuffix}.txt ]; then
    rm ${outdir}/allchrom_${outsuffix}.txt
    echo "The output file already existed and was deleted before rewriting."
fi

# Concatenate the input files (skipping the header in all files except the first):
awk '
{
    # If this is the first line of the file:
    if (FNR == 1) {
        # If this is the first file:
        if (NR == 1) {
            # Print the line:
            print $0
        }
        # If this is not the first file:
        else {
            # Skip the line.
        }
    }
    # If this is not the first line of the file:
    else {
        # Print the line:
        print $0
    }
}
' $(awk -v indir=${indir} -v insuffix=${insuffix} 'BEGIN {OFS=""; ORS=" "} {print indir, "/", $1, "_", insuffix, ".txt"}' ${genome_file}) >> ${outdir}/allchrom_${outsuffix}.txt

echo "Step 4: done (`date`)."
