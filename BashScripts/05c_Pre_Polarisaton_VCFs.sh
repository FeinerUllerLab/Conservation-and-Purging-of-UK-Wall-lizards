#!/bin/bash -l

# Objectve: Meging of ingroup and outgroup VCF for polarisation.

# Note that the outgroup.vcf was cretef following steps 1 to 5 (During filteting, multiallelic sites were kept). Mapping was done to the P. muralis reference genome. 

# Module
module load bcftools/1.20

# Paths 
INGROUP="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Final_dataset/All_Origins_SNP_final.vcf.gz"
OUTGROUPS_VCF="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/Polarisation/All_Polarisation_final.vcf.gz"
OUTDIR="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Pol_VCFs"
REF_FASTA="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/RefGenome/GCA_004329235.1_PodMur_1.0_genomic.fna"  

mkdir -p "${OUTDIR}"

# Intermediate files
INGROUP_NOMISS="${OUTDIR}/ingroup.nomiss.vcf.gz"
INGROUP_NORM="${OUTDIR}/ingroup.norm.vcf.gz"
OUTGROUP_NORM="${OUTDIR}/outgroup.norm.vcf.gz"
MERGED_ALL="${OUTDIR}/merged_all.vcf.gz"
MERGED_SITES="${OUTDIR}/merged_sites_in_ingroup.vcf.gz"
FINAL="${OUTDIR}/Pre_polarisation.vcf.gz"

# Step 1: Remove missing genotypes in ingroup
bcftools filter -e 'F_MISSING > 0' "${INGROUP}" -Oz -o "${INGROUP_NOMISS}"
bcftools index -t "${INGROUP_NOMISS}"

# Step 2: Normalize alleles
bcftools norm -f "${REF_FASTA}" -c s "${INGROUP_NOMISS}" -Oz -o "${INGROUP_NORM}"
bcftools norm -f "${REF_FASTA}" -c s "${OUTGROUPS_VCF}" -Oz -o "${OUTGROUP_NORM}"
bcftools index -t "${INGROUP_NORM}"
bcftools index -t "${OUTGROUP_NORM}"

# Step 3: Merge normalized VCFs
bcftools merge -m all -Oz -o "${MERGED_ALL}" "${INGROUP_NORM}" "${OUTGROUP_NORM}"
bcftools index -t "${MERGED_ALL}"

# Step 4: Keep only ingroup sites
bcftools view -T "${INGROUP_NORM}" "${MERGED_ALL}" -Oz -o "${MERGED_SITES}"
bcftools index -t "${MERGED_SITES}"

# Step 5: Remove symbolic alleles (*) and non-ACGT sites
bcftools view -e 'ALT="*"' "${MERGED_SITES}" -Oz -o "${FINAL}"
bcftools index -t "${FINAL}"

# Step 6: QC stats
bcftools stats "${MERGED_SITES}" > "${MERGED_SITES%.vcf.gz}.stats"
bcftools stats "${FINAL}" > "${FINAL%.vcf.gz}.stats"
echo "Check ${FINAL%.vcf.gz}.stats for missingness and site counts"
