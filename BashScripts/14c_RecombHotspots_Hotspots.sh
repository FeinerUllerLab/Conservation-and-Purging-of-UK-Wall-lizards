#!/bin/bash

INDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step2_Filtered/"
OUTDIR="/home/feiner/Projects/UKwallies/RecombHotSpots/Step3_Hotspots/"
mkdir -p $OUTDIR

WINDOW=2000          # window size in bp
BG=40000             # background window size (±40 kb)
THRESHOLD=5          # hotspot threshold

echo "Calculating hotspots and counting hotspot/background windows..."

for FILE in ${INDIR}/*_filtered.txt
do
    BASENAME=$(basename "$FILE")
    OUTFILE="${OUTDIR}/${BASENAME/.txt/_hotspots.txt}"
    echo "Processing $BASENAME"

    awk -v W=$WINDOW -v BG=$BG -v TH=$THRESHOLD '
    BEGIN{OFS="\t"}
    NR==1{
        # header
        print $0, "hotspot"
        next
    }
    {
        # store data in arrays for background calculation
        chr[NR]=$1
        start[NR]=$2
        end[NR]=$3
        mean_rate[NR]=$4
        sites[NR]=$5
        status[NR]=$6
    }
    END{
        n=NR
        bg_win=int(BG/W)   # how many windows left/right for background
        hotspot_count=0
        background_count=0

        for(i=1;i<=n;i++){
            # skip filtered windows
            if(status[i]=="filter"){
                label="NA"
            } else {
                # compute background mean excluding the current window
                sum=0
                count=0
                for(j=i-bg_win;j<=i+bg_win;j++){
                    if(j>0 && j<=n && j!=i && status[j]=="pass"){
                        sum+=mean_rate[j]
                        count++
                    }
                }
                bg_rate=(count>0 ? sum/count : 0)

                if(bg_rate>0 && mean_rate[i]>=TH*bg_rate){
                    label="hotspot"
                    hotspot_count++
                } else {
                    label="background"
                    background_count++
                }
            }
            print chr[i], start[i], end[i], mean_rate[i], sites[i], status[i], label
        }

        # print counts to standard output
        printf "Chromosome %s: windows hotspot=%d, background=%d\n", chr[1], hotspot_count, background_count > "/dev/stderr"
    }' $FILE > $OUTFILE

done

echo "All chromosomes processed. Hotspot files in $OUTDIR"
