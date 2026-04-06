#!/bin/bash -l

# Objective: Extract SNPs from the All-sites dataset and generate an only-variants dataset
# Modules
module load bcftools/1.20

# Define characters  
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data/Filtering_steps 
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Final_dataset
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Steps_Stats
statsdir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Steps_Stats

# Merge both populaions VCFs 
bcftools merge -Oz -o ${datadir}/Merged.vcf.gz ${datadir}/All_French_filtered.vcf.gz ${datadir}/All_Italian_filtered.vcf.gz
bcftools index ${savedir}/Merged.vcf.gz

# Extract only SNPs (remove indels and multi-allelic sites)
bcftools view -m2 -M2 -v snps -Oz -o ${savedir}/Merged_snps.vcf.gz ${savedir}/Merged.vcf.gz
bcftools index ${savedir}/Merged_snps.vcf.gz

echo "SNPs extracted"

# Remove the intermediate Merged VCF to save space
rm ${savedir}/Merged.vcf.gz ${savedir}/Merged.vcf.gz.csi  

# SNPs filtering. 
#Hard-filtering pass
bcftools view ${savedir}/Merged_snps.vcf.gz -f 'PASS,.' -Ob -o ${savedir}/Merged_snps_pass.bcf
bcftools index ${savedir}/Merged_snps_pass.bcf

#for minor allele frequency:

bcftools view -q 0.05:minor ${savedir}/Merged_snps_pass.bcf -Ob -o ${savedir}/Merged_snps_pass_maf0.05.bcf
bcftools index ${savedir}/Merged_snps_pass_maf0.05.bcf

#for genotype quality (GQ) below 20:
bcftools filter -S . -e 'FMT/GQ<20' ${savedir}/Merged_snps_pass_maf0.05.bcf -Ob -o ${savedir}/Merged_snps_pass_maf0.05_GQ20.bcf
bcftools index ${savedir}/Merged_snps_pass_maf0.05_GQ20.bcf

#delete SNPs with any missing genotype data
bcftools filter -e 'F_MISSING > 0' ${savedir}/Merged_snps_pass_maf0.05_GQ20.bcf -Ov -o ${savedir}/All_Origins_SNP_final.vcf
bcftools stats ${savedir}/All_Origins_SNP_final.vcf > ${savedir}/All_Origins_SNP_final.stats

echo "Filters applied"




