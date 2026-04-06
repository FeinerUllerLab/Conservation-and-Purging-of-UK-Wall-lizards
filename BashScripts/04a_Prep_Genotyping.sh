#!/bin/bash -l 

# Objective: Generate necessary list files for Genotyping 

# Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping

# Define sample group files
Italian="Italian"  # List of samples for group 1
French="French"  # List of samples for group 2
chromosomes="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/Chromosome_autosomal.txt"  # Chromosome names

# Check if the chromosome file exists
if [[ ! -f $chromosomes ]]; then
    echo "Error: Chromosome file $chromosomes not found!" >&2
    exit 1
fi

# Create GVCF lists for each group and chromosome
for group in Italian; do
    groupfile=${!group}  # Choose one of the groups and run (Italian or French)

    # Check if the sample file exists
    if [[ ! -f $groupfile ]]; then
        echo "Error: Sample file $groupfile not found for group $group!" >&2
        exit 1
    fi

    while read -r chr; do
        gvcf_list="${savedir}/GVCFs_${chr}_${group}.list"
        > $gvcf_list  # Clear or create the list file

        while read -r sample; do
            gvcf="${datadir}/${sample}/ROH_HaplotypeCaller/${sample}_ROH_${chr}.g.vcf.gz"
            if [[ -f $gvcf ]]; then
                realpath "$gvcf" >> $gvcf_list
            else
                echo "Warning: GVCF not found for $sample, $chr" >&2
            fi
        done < "$groupfile"
    done < "$chromosomes"
done

echo "GVCF lists created successfully!"
