#!/bin/bash -l

# Objective: Merge the variants of all chromosomes and apply filters

#Load modules 
module load bcftools/1.20
module load gatk/4.5.0.0
PICARD_HOME=/pdc/software/eb/software/picard/2.25.5

#Define working directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data/Filtering_steps
filtdir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data/Filtering_steps
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Diversity_Data
reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna

date 

# Merging all the chromosomes. 
# West-Europe 
java -jar $PICARD_HOME/picard.jar GatherVcfs \
-O ${filtdir}/All_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014743.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014744.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014745.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014746.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014747.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014748.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014749.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014750.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014751.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014752.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014753.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014754.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014755.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014756.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014757.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014758.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014759.1_West-Europe_filtered.vcf.gz \
-I ${datadir}/CM014760.1_West-Europe_filtered.vcf.gz 

# Index
gatk IndexFeatureFile -I ${filtdir}/All_West-Europe_filtered.vcf.gz 
bcftools index ${filtdir}/All_West-Europe_filtered.vcf.gz

echo "West-Europe Whole-genome  VCF completed" 

# Central-Italy 
java -jar $PICARD_HOME/picard.jar GatherVcfs \
-O ${filtdir}/All_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014743.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014744.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014745.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014746.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014747.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014748.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014749.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014750.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014751.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014752.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014753.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014754.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014755.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014756.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014757.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014758.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014759.1_Central-Italy_filtered.vcf.gz \
-I ${datadir}/CM014760.1_Central-Italy_filtered.vcf.gz 

# Index
gatk IndexFeatureFile -I ${filtdir}/All_Central-Italy_filtered.vcf.gz 
bcftools index ${filtdir}/All_Central-Italy_filtered.vcf.gz

echo "Central-Italy Whole-genome  VCF completed"

# keep only 'pass' variants and index again (Hard-filtering):
#West-Europe 
bcftools view ${filtdir}/All_West-Europe_filtered.vcf.gz -f 'PASS,.' -Ob -o ${filtdir}/All_West-Europe_filtered_pass.bcf
bcftools index ${filtdir}/All_West-Europe_filtered_pass.bcf
#Central-Italy
bcftools view ${filtdir}/All_Central-Italy_filtered.vcf.gz -f 'PASS,.' -Ob -o ${filtdir}/All_Central-Italy_filtered_pass.bcf
bcftools index ${filtdir}/All_Central-Italy_filtered_pass.bcf

echo "Final filters starting"

# Further filtering (Missing genotype data must be specially hard to not bias the Gene diversity analysis)

# For genotype quality (GQ) below 20:
bcftools filter -S . -e 'FMT/GQ<20' ${filtdir}/All_West-Europe_filtered_pass.bcf -Ob -o ${filtdir}/All_West-Europe_filtered_pass_GQ20.bcf
bcftools index ${filtdir}/All_West-Europe_filtered_pass_GQ20.bcf

bcftools filter -S . -e 'FMT/GQ<20' ${filtdir}/All_Central-Italy_filtered_pass.bcf -Ob -o ${filtdir}/All_Central-Italy_filtered_pass_GQ20.bcf
bcftools index ${filtdir}/All_Central-Italy_filtered_pass_GQ20.bcf

echo "Genotype quality applied"

# Remove multialellic sites
bcftools view -M2 ${filtdir}/All_West-Europe_filtered_pass_GQ20.bcf -Ob -o ${filtdir}/All_West-Europe_filtered_pass_GQ20_biallelic.bcf
bcftools index ${filtdir}/All_West-Europe_filtered_pass_GQ20_biallelic.bcf

bcftools view -M2 ${filtdir}/All_Central-Italy_filtered_pass_GQ20.bcf -Ob -o ${filtdir}/All_Central-Italy_filtered_pass_GQ20_biallelic.bcf
bcftools index ${filtdir}/All_Central-Italy_filtered_pass_GQ20_biallelic.bcf

echo "Multiallelic sites removed"


# Delete variants with ANY missing genotype data
bcftools filter -e 'F_MISSING > 0' ${filtdir}/All_West-Europe_filtered_pass_GQ20_biallelic.bcf -Ov -o ${savedir}/All_West-Europe_final.vcf
bcftools stats ${savedir}/All_West-Europe_final.vcf > ${savedir}/All_West-Europe_final.stats

bcftools filter -e 'F_MISSING > 0' ${filtdir}/All_Central-Italy_filtered_pass_GQ20_biallelic.bcf -Ov -o ${savedir}/All_Central-Italy_final.vcf
bcftools stats ${savedir}/All_Central-Italy_final.vcf > ${savedir}/All_Central-Italy_final.stats

echo "Filtering complete" 
date 
