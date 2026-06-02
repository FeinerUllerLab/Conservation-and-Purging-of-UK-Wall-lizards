#! /bin/bash -l

#SBATCH --job-name=RoH_Het
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --mem=125G
#SBATCH --error=logs/RoH_Het.%J.err
#SBATCH --output=logs/RoH_Het.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=feiner@evolbio.mpg.de
#SBATCH --partition=standard

### To obtain Runs of Homozygosity using bcftools and plink and calculate the FROH for each sample (inbreeding coefficient)

# Note that the script needs to be adapted to run on a cluster by either loading the right modules (on Dardel) or by giving the right path to executables (on Wallace)
# Export path to bcftools:
export PATH=/data/biosoftware/bcftools/bcftools-1.21/:$PATH

# VCF file with variant and invariant sites needs to be declared; this has been filtered to exclude sites with any missing genotypes and is used for RoH estimation

VCF=/home/feiner/Projects/UKwallies/Datasets/AllSites_VCFs/All_West-Europe_final.vcf.gz

# A file with the samples' names needs to be declared (incl. outgroup)
samples_list=/home/feiner/Projects/UKwallies/scripts/W-EUR_samples

# The size cut-off for RoH is declared
RoH_size_num=2000000
RoH_size_nam=2Mb

mkdir -p /home/feiner/Projects/UKwallies/Results_RoH_Het_W-EUR
cd /home/feiner/Projects/UKwallies/Results_RoH_Het_W-EUR

#### Add a fictitious fully homozygous sample (Ppi_Hz) to calculate the length of the homozygous genome
# copy the header
#zcat ${VCF} | head -n 10000 | grep "^#" > head.tmp
# add new sample name
#sed -i '$ s/$/\tPxx_Hz/' head.tmp
# add 0/0 genotypes
#zcat ${VCF} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | bgzip > ${VCF}_pseudo
#rm head.tmp

#### BCFTOOLS
# -e tells it to estimate frequencies based only on the actual samples - not the fictitious one
#bcftools roh ${VCF}_pseudo -e ${samples_list} -G 30 -O r -o roh.pseudo
# quality filters: Phred score > 50 & length > 1 Mb
#cat roh.pseudo | awk '$8 > 50' > roh.pseudo.qual
awk -v minsize="$RoH_size_num" '$6 > minsize' roh.pseudo.qual > roh.pseudo.qual.W-EUR.${RoH_size_nam}

#### Calculate FRoH for each sample: (sum of ROHs > minsize)/(length of homozygous genome) and export relevant data
#Calculate length of the homozygous genome
homlength=$(grep "Pxx_Hz" roh.pseudo.qual.W-EUR.${RoH_size_nam} | awk '{sum += $6} END {print sum}')
out="froh_summary_bcftools_W-EUR_${RoH_size_nam}.txt"
echo -e "Sample\tFRoH\tLength\tnRoH" > "$out"
while read -r sample; do
    [ -z "$sample" ] && continue

    # Sum of ROH lengths for this sample
    sum=$(awk -v s="$sample" '$2 == s {total += $6} END {print total+0}' roh.pseudo.qual.W-EUR.${RoH_size_nam})

    # Count number of ROH segments for this sample
    count=$(awk -v s="$sample" '$2 == s {n++} END {print n+0}' roh.pseudo.qual.W-EUR.${RoH_size_nam})

    # Calculate FROH
    ratio=$(awk -v a="$sum" -v b="$homlength" 'BEGIN {
        if (b > 0) printf "%.6f", a / b;
        else print "NA"
    }')
    echo -e "${sample}\t${ratio}\t${sum}\t${count}" >> "$out"
done < "$samples_list"

#### Heterozygosity - genome-wide (remove outgroup)
#plink --vcf ${VCF} -aec --double-id --make-bed --out final.all.plink
#plink -bfile final.all.plink -aec --double-id --het --out het
#plink -bfile final.all.plink -aec --double-id --missing --out miss

#collect output in one file
#Het_sum_geno="heterozygosity_summary_genomewide_W-EUR.tsv"
#echo -e "Sample\tNumberVariableSites\tObservedHomozygous\tMissingSites\tGenotypedSites" > "$Het_sum_geno"
#paste \
#  <(awk 'NR>1 {print $1 "\t" $5 "\t" $3}' het.het) \
#  <(awk 'NR>1 {print $4 "\t" $5}' miss.imiss) \
#>> "$Het_sum_geno"

echo "Summary written to $Het_sum_geno"

#### Heterozygosity - outside of RoHs for each individual
# Continue with the output of BCFTOOLS to identify RoHs for each individual (create .bed file) and use this to exclude from plink het step

ROH="roh.pseudo.qual.W-EUR.${RoH_size_nam}"

# 1) Create empty BED files for all samples
while read -r sample; do
    [ -z "$sample" ] && continue
    : > "RoH_position_${sample}.bed"
done < "$samples_list"
# 2) Fill BED files where ROH data exist
awk -v samples="$samples_list" ' BEGIN {
    # read sample list
    while ((getline < samples) > 0) {
        gsub(/\r/, "", $0)
        keep[$0] = 1
    }
    close(samples)
}
NR <= 3 { next } # skip header lines
# if sample in column 2 is in the list
($2 in keep) {
    outfile = "RoH_position_" $2 ".bed"
    # BED columns: chr, start, end, sample
    print $3, $4, $5, $2 >> outfile
}
' "$ROH"

### now run Plink for each sample
export outgroup
cat ${samples_list} | parallel -j 6 \
'plink -bfile final.all.plink -aec --double-id --exclude range RoH_position_{}.bed --het --out het_{}'
cat ${samples_list} | parallel -j 6 \
'plink -bfile final.all.plink -aec --double-id --exclude range RoH_position_{}.bed --missing --out miss_{}'

# Then collect sample-specific data from all generated Plink outputs
Het_sum_nonRoH="heterozygosity_summary_nonRoH_W-EUR_${RoH_size_nam}.tsv"
echo -e "Sample\tNumberVariableSites\tObservedHomozygous\tMissingSites\tGenotypedSites" > "$Het_sum_nonRoH"

while read -r sample; do
    [ -z "$sample" ] && continue

    MISS_FILE="miss_${sample}.imiss"
    HET_FILE="het_${sample}.het"

    if [ ! -f "$MISS_FILE" ] || [ ! -f "$HET_FILE" ]; then
        echo "WARNING: Missing files for sample $sample, skipping"
        continue
    fi

    # NumberVariableSites: from het file, column 5 ("N(NM)"), row matching sample
    VAR_SITES=$(awk -v col=5 -v samp="$sample" 'NR>1 && $1==samp {print $col}' "$HET_FILE")

    # ObservedHomozygous: from het file, column 3 ("O(HOM)"), row matching sample
    OBS_HOM=$(awk -v col=3 -v samp="$sample" 'NR>1 && $1==samp {print $col}' "$HET_FILE")

    # MissingSites: from miss file, column 4 ("N_MISS"), row matching sample
    N_MISS=$(awk -v col=4 -v samp="$sample" 'NR>1 && $1==samp {print $col}' "$MISS_FILE")

    # GenotypedSites: from miss file, column 5 ("N_GENO"), row matching sample
    N_GENO=$(awk -v col=5 -v samp="$sample" 'NR>1 && $1==samp {print $col}' "$MISS_FILE")

    echo -e "${sample}\t${VAR_SITES}\t${OBS_HOM}\t${N_MISS}\t${N_GENO}" >> "$Het_sum_nonRoH"

done < "$samples_list"

echo "Summary written to $Het_sum_nonRoH"

#remove intermediate files (careful with with if trouble-shooting is needed)
#rm het_*
#rm miss_*
#rm RoH_position_*.bed

echo "All done."
