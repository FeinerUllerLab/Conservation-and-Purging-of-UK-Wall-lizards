#!/bin/bash -l

# Objective: Array for Genotyping and hard-filtering in Central-Italy-origin samples

#Modules
module load gatk/4.5.0.0
module load bcftools/1.20

# Define the variables 
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/ROH_JointGenotyping
reference=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna
chrom_list="./Chromosome_autosomal.txt"  

# Array the chromosomes
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chrom_list")
group="Central-Italy"

# The list of GVCF was created previously depending on the origin of the samples (Central-Italy or West-Europe)
date 
echo "${chr}_in process" 

#Create the temporary directories: 
 
mkdir -p ${savedir}/tmp_dir_comb_${chr}_${group}
mkdir -p ${savedir}/tmp_dir_geno_${chr}_${group}


#Combination of GVFs 
gatk --java-options "-Xmx8g" CombineGVCFs \
-R ${reference} \
-V ${savedir}/GVCFs_${chr}_${group}.list \
-O ${savedir}/${chr}_${group}.g.vcf.gz \
--tmp-dir ${savedir}/tmp_dir_comb_${chr}_${group}

echo "Combining GVCFs completed"

#Genotyping. 
gatk --java-options "-Xmx8g" GenotypeGVCFs \
-R ${reference} \
-V ${savedir}/${chr}_${group}.g.vcf.gz \
-O ${savedir}/${chr}_${group}.vcf.gz \
--include-non-variant-sites \
--tmp-dir ${savedir}/tmp_dir_geno_${chr}_${group}

echo "Genotyping completed"

#Remove temporary directories/files 
rm ${savedir}/GVCFs_${chr}_${group}.list
rm -rf ${savedir}/tmp_dir_comb_${chr}_${group}
rm -rf ${savedir}/tmp_dir_geno_${chr}_${group}
# rm ${savedir}/${chr}.g.vcf.gz*

#Select for all sites but indel 
bcftools view -e 'STRLEN(ALT) > 1 || STRLEN(REF) > 1' ${savedir}/${chr}_${group}.vcf.gz -Oz -o ${savedir}/${chr}_${group}_allsites.vcf.gz 
bcftools index ${savedir}/${chr}_${group}_allsites.vcf.gz --tbi

echo "All but indels done"

#Select Indel Variants
gatk --java-options "-Xmx8g" SelectVariants \
-V ${savedir}/${chr}_${group}.vcf.gz \
-O ${savedir}/${chr}_${group}_INDEL.vcf.gz \
-select-type INDEL
gatk IndexFeatureFile -I ${savedir}/${chr}_${group}_INDEL.vcf.gz

echo "Indels found" 

# Extract the "SNP" sites (with ALT alleles). 
bcftools view ${savedir}/${chr}_${group}_allsites.vcf.gz --min-alleles 2 -Oz -o ${savedir}/${chr}_${group}_snps.vcf.gz
bcftools index ${savedir}/${chr}_${group}_snps.vcf.gz --tbi

echo "SNPs found"

# Extract the invariant sites (no ALT sites). 
bcftools view ${savedir}/${chr}_${group}_allsites.vcf.gz --max-alleles 1 -Oz -o ${savedir}/${chr}_${group}_invariant.vcf.gz
bcftools index ${savedir}/${chr}_${group}_invariant.vcf.gz --tbi

echo "Invariant sites found" 


echo "Applying filters"

# Hard-Filtering for SNPs - Include QUAL and QD expressions /  Dp estimated by #samples (18 Central-Italy) x Desired covergae (11x to 44x)

gatk --java-options "-Xmx8g" VariantFiltration \
-R ${reference} \
-V "${savedir}/${chr}_${group}_snps.vcf.gz" \
--mask "${savedir}/${chr}_${group}_INDEL.vcf.gz" \
--mask-extension 10 \
--mask-name InDel \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "FS > 60.0" --filter-name "FisherStrand" \
--filter-expression "MQ < 40.0" --filter-name "MapQual40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--filter-expression "DP < 198"  --filter-name "minDepth" \
--filter-expression "DP > 792"  --filter-name "maxDepth" \
-O ${savedir}/${chr}_${group}_snps_filtered.vcf.gz

echo "SNPs filtered"

# Hard-Filtering for invariant sites - Remove QUAL and QD expressions /  Dp estimated by #samples (18 Central-Italy) x Desired covergae (11x to 44x)

gatk --java-options "-Xmx8g" VariantFiltration \
-R ${reference} \
-V "${savedir}/${chr}_${group}_invariant.vcf.gz" \
--mask "${savedir}/${chr}_${group}_INDEL.vcf.gz" \
--mask-extension 10 \
--mask-name InDel \
--filter-expression "FS > 60.0" --filter-name "FisherStrand" \
--filter-expression "MQ < 40.0" --filter-name "MapQual40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--filter-expression "DP < 198"  --filter-name "minDepth" \
--filter-expression "DP > 792"  --filter-name "maxDepth" \
-O ${savedir}/${chr}_${group}_invariant_filtered.vcf.gz

echo "Invariant sites filtered"

# Index after the Hard-filtering
gatk IndexFeatureFile -I ${savedir}/${chr}_${group}_snps_filtered.vcf.gz
gatk IndexFeatureFile -I ${savedir}/${chr}_${group}_invariant_filtered.vcf.gz

# Merging Variant and invariant sites
bcftools concat \
 --allow-overlaps \
 ${savedir}/${chr}_${group}_snps_filtered.vcf.gz \
 ${savedir}/${chr}_${group}_invariant_filtered.vcf.gz \
 -O z \
 -o ${savedir}/${chr}_${group}_filtered.vcf.gz

gatk IndexFeatureFile -I ${savedir}/${chr}_${group}_filtered.vcf.gz

echo "HardFiltering and merging of ${chr} completed"
date