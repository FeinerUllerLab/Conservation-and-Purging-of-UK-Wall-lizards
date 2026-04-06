#!/bin/bash

#SBATCH --job-name=SnpEff
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --mem=125G
#SBATCH --error=logs/SnpEff_%A_%a.err
#SBATCH --output=logs/SnpEff_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=feiner@evolbio.mpg.de
#SBATCH --partition=standard

# Objective: Annotate SNP dataset using SnpEff

# Activate the mamaba environment (Be sure that it is active when queueing the job)
#mamba activate snpeff
#module load bcftools/1.20

module load java/x64/24u2


# Define directories and genome name
datadir=/home/feiner/Projects/UKwallies/Datasets/Pol_VCFs/
mkdir -p /home/feiner/Projects/UKwallies/Purging
savedir=/home/feiner/Projects/UKwallies/Purging
Genome_Database="PodMur_1.0#1"

date 

# Rename chromosome names to match the snpeff database
bcftools annotate --rename-chrs Purging_Chr_names_NCBI.txt ${datadir}/All_Origins_polarised.vcf.gz -Oz -o ${savedir}/All_Polarised_Purging_tmp.vcf.gz 

bcftools index --tbi ${savedir}/All_Polarised_Purging_tmp.vcf.gz 

# Run SnpEff on the origin subsets
java -jar /data/biosoftware/SnpEff/snpEff/snpEff.jar ann -dataDir /mnt/beegfs/databases/snpEff_data \
-stats ${savedir}/All_Polarised_Annotated_NCBI.html -csvStats ${savedir}/All_Polarised_Annotated_NCBI.csv \
${Genome_Database} ${savedir}/All_Polarised_Purging_tmp.vcf.gz > ${savedir}/All_Polarised_annotated_NCBI.vcf

bgzip ${savedir}/All_Polarised_annotated_NCBI.vcf
bcftools index --tbi ${savedir}/All_Polarised_annotated_NCBI.vcf.gz 

rm ${savedir}/All_Polarised_Purging_tmp.*

echo "SnpEff annotation complete."

date
