#!/bin/bash -l

# Objective: Perfom a PCA using the filtered only-variants VCF file. 

#Module 
module load plink/2.00a5.14
module load bcftools/1.20

#Define directories: 
infodir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Final_dataset
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/PCA
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/PCA

#Change names of chromosomes to valid values: 
bcftools annotate --rename-chrs ${infodir}/New_Chr_names.txt ${datadir}/All_Origins_SNP_final.vcf.gz -Oz -o ${savedir}/All_Origins_SNP_final_tmp.vcf.gz 
bcftools index ${savedir}/All_Origins_SNP_final_tmp.vcf.gz 

# Convert VCF to PLINK binary format
plink --vcf ${savedir}/All_Origins_SNP_final_tmp.vcf.gz  --make-bed --out ${savedir}/All_Origins_SNP_final_plink --allow-extra-chr

#Run the PCA
plink -bfile ${savedir}/All_Origins_SNP_final_plink --aec --pca --out ${savedir}/All_Origins_PCA_FINAL

rm ${savedir}/All_Origins_SNP_final_tmp.*