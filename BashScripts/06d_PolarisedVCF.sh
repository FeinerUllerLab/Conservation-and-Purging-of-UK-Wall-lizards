#!/bin/bash -l

# Objective: Produced a filtered VCF with only ancestral REF alleles for annotation and downstream analyses 

# Load modules
module load bcftools/1.20

# Set variables
INPUT_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Final_dataset/All_Origins_SNP_final.vcf.gz"
BED_FILE="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr/High_confidence_ref_ancestral.bed"
OUTPUT_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Pmuralis/All_Origins_polarised.vcf.gz"

echo "Starting VCF filtering"

# Filter VCF using bcftools
bcftools view -T ${BED_FILE} ${INPUT_VCF} -Oz -o ${OUTPUT_VCF}

# Index the filtered VCF
bcftools index --tbi ${OUTPUT_VCF}

# Stats on the final VCF 
bcftools stats ${OUTPUT_VCF} > ${OUTPUT_VCF%.vcf.gz}.stats

echo "VCF filtering complete!"
echo "Filtered VCF created: ${OUTPUT_VCF}"
