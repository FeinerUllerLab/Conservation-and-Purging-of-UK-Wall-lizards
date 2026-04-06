#!/bin/bash

INDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step1_WindowSummary/"
OUTDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step2_Filtered/"
mkdir -p $OUTDIR

echo "Filtering windows with sites == 0 and reporting counts..."

for FILE in ${INDIR}/*_window_summary.txt
do
    BASENAME=$(basename "$FILE")
    OUTFILE="${OUTDIR}/${BASENAME/.txt/_filtered.txt}"

    echo "Processing $BASENAME"

    # Initialize counters
    PASS=0
    FILTER=0

    # Filter windows and write output
    awk -v pass_ref="$PASS" -v filter_ref="$FILTER" '
    BEGIN{OFS="\t"}
    NR==1{
        print $0, "status"
        next
    }
    {
        if($5 == 0){
            status="filter"
            filter_ref++
        } else {
            status="pass"
            pass_ref++
        }
        print $0, status
    }
    END{
        # Print counts to stdout
        printf "  Windows passed: %d, Windows filtered: %d\n", pass_ref, filter_ref > "/dev/stderr"
    }' $FILE > $OUTFILE

done

echo "Done. Filtered files written to $OUTDIR"
