#!/bin/bash

# Objective: Annotate SNP dataset using SnpEff

# Activate the mamaba environment (Be sure that it is active when queueing the job)
mamba activate snpeff
module load bcftools/1.20

# Define directories and genome name
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Pmuralis
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Purging
Genome_Database="PodMur_1.0.99"

date 

# Rename chromosome names to match the snpeff database
bcftools annotate --rename-chrs Purging_Chr_names.txt ${datadir}/All_Origins_polarised.vcf.gz -Oz -o ${savedir}/All_Polarised_Purging_tmp.vcf.gz 
bcftools index --tbi ${savedir}/All_Polarised_Purging_tmp.vcf.gz 

# Run SnpEff on the origin subsets
snpEff eff -v ${Genome_Database} ${savedir}/All_Polarised_Purging_tmp.vcf.gz   \
 -stats ${savedir}/All_Polarised_Annotated.html -csvStats ${savedir}/All_Polarised_Annotated.csv > ${savedir}/All_Polarised_annotated.vcf

bgzip ${savedir}/All_Polarised_annotated.vcf
bcftools index --tbi ${savedir}/All_Polarised_annotated.vcf.gz 

rm ${savedir}/All_Polarised_Purging_tmp.*

echo "SnpEff annotation complete."

date