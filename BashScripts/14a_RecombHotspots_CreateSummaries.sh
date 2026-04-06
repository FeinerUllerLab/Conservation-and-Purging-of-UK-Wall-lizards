#!/bin/bash

# Directories
INDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/RawMaps/"
OUTDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step1_WindowSummary/"
mkdir -p $OUTDIR

# Sliding window parameters
WINDOW=2000   # 2 kb
STEP=1000     # 1 kb

# Chromosome lengths file
CHR_LEN_FILE="Chr_Length.txt"

# Loop over each chromosome
while read -r chr dummy len
do
    echo "Processing chromosome $chr (length $len bp)"

    # Input rmap file
    RMAP_FILE="${INDIR}/${chr}_VO01_n60_p2_b50_w100.rmap"
    if [ ! -f "$RMAP_FILE" ]; then
        echo "  WARNING: rmap file not found: $RMAP_FILE"
        continue
    fi

    # Output summary file
    OUTFILE="${OUTDIR}/${chr}_window_summary.txt"
    echo -e "chr\tstart\tend\tmean_rate\tsites" > $OUTFILE

    # Process using awk
awk -v W=$WINDOW -v S=$STEP -v CHR=$chr -v CHR_LEN=$len '
BEGIN{
    maxwin=int((CHR_LEN-W)/S)

    for(w=0; w<=maxwin; w++){
        win_start[w]=w*S
        win_end[w]=win_start[w]+W
        sum_rate[w]=0
        total_len[w]=0
        sites[w]=0   # explicitly initialize
    }
}
{
    seg_start=$1
    seg_end=$2
    rate=$3

    # determine only relevant windows (IMPORTANT SPEED + CORRECTNESS)
    first_w = int(seg_start / S)
    last_w  = int(seg_end   / S)

    for(w=first_w; w<=last_w; w++){

        # skip invalid windows
        if(w < 0 || w > maxwin) continue

        overlap_start = (seg_start > win_start[w] ? seg_start : win_start[w])
        overlap_end   = (seg_end   < win_end[w]   ? seg_end   : win_end[w])

        if(overlap_start < overlap_end){
            ol_len = overlap_end - overlap_start

            sum_rate[w] += rate * ol_len
            total_len[w] += ol_len


if(seg_start >= win_start[w] && seg_start < win_end[w]){
    sites[w]++
}

        }
    }
}
END{
    for(w=0; w<=maxwin; w++){
        mean_rate = (total_len[w] > 0 ? sum_rate[w]/total_len[w] : 0)

        # enforce true zero if no overlap
        if(total_len[w] == 0){
            sites[w] = 0
        }

        print CHR, win_start[w], win_end[w], mean_rate, sites[w]
    }
}
' $RMAP_FILE >> $OUTFILE

done < $CHR_LEN_FILE

echo "Step 1 completed. Summary files in $OUTDIR"
