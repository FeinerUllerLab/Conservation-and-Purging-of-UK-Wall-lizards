#!/bin/bash

INDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step3_Hotspots/"
OUTDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step4_MergedHotspots/"
mkdir -p $OUTDIR

echo "Merging neighboring/overlapping hotspot windows..."

for FILE in ${INDIR}/*_hotspots.txt
do
    BASENAME=$(basename "$FILE")
    OUTFILE="${OUTDIR}/${BASENAME/.txt/_merged.txt}"
    echo "Processing $BASENAME"

    awk '
    BEGIN{OFS="\t"}
    NR==1{next}  # skip header

    {
        chr=$1
        start=$2
        end=$3
        mean_rate=$4
        sites=$5
        status=$6
        hotspot=$7

        # only process hotspots
        if(hotspot!="hotspot") next

        # start new block if needed
        if(!inblock){
            block_chr=chr
            block_start=start
            block_end=end
            sum_rate=mean_rate
            count_rate=1
            window_count=1
            inblock=1
        } else {
            # check if current window is adjacent or overlapping previous
            if(start <= block_end){ 
                # merge into current block
                block_end=(end>block_end?end:block_end)
                sum_rate+=mean_rate
                count_rate++
                window_count++
            } else {
                # print previous block
                print block_chr, block_start, block_end, sum_rate/count_rate, window_count
                # start new block
                block_chr=chr
                block_start=start
                block_end=end
                sum_rate=mean_rate
                count_rate=1
                window_count=1
            }
        }
    }
    END{
        # print last block
        if(inblock) print block_chr, block_start, block_end, sum_rate/count_rate, window_count
    }' $FILE > $OUTFILE

    # Count number of merged blocks and windows merged
    merged_blocks=$(wc -l < $OUTFILE)
    total_windows=$(awk '$7=="hotspot"{c++}END{print c+0}' $FILE)

    echo "Chromosome $BASENAME: merged $total_windows windows into $merged_blocks blocks"

done

echo "All chromosomes processed. Merged hotspot files in $OUTDIR"

### now create a bed file with all hotspots

MERGED_DIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step4_MergedHotspots/"
OUTFILE="/home/feiner/Projects/UKwallies/RecombHotSpots/All_Merged_Hotspots.bed"

echo -e "chr\tstart\tend\tmean_rate\twindows_merged" > $OUTFILE

for FILE in ${MERGED_DIR}/*_merged.txt
do
    # Skip header line and append to the BED file
    tail -n +1 $FILE >> $OUTFILE
done

echo "All chromosomes merged into $OUTFILE"
