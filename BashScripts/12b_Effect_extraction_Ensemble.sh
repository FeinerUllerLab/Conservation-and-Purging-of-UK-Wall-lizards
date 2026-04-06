#!/bin/bash -l

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

# Objective: Separate variants based on the annotation impacts using SnpSift

# Modules and environments (make sure that the enviroment is loaded while the job is running) 
#conda activate snpeff 
#module load bcftools/1.20

module load java/x64/24u2

# Define directories and genome name
datadir=/home/feiner/Projects/UKwallies/Purging
mkdir -p /home/feiner/Projects/UKwallies/Purging/Polarised_Impacts_Ensemble
savedir=/home/feiner/Projects/UKwallies/Purging/Polarised_Impacts_Ensemble

date 

# Filter annotoations based on impact 
java -jar /data/biosoftware/SnpEff/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'HIGH')" ${datadir}/All_Polarised_annotated_Ensemble.vcf.gz > ${savedir}/All_Polarised_high.vcf

bgzip ${savedir}/All_Polarised_high.vcf
bcftools index --tbi ${savedir}/All_Polarised_high.vcf.gz

java -jar /data/biosoftware/SnpEff/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'MODERATE')" ${datadir}/All_Polarised_annotated_Ensemble.vcf.gz > ${savedir}/All_Polarised_moderate.vcf

bgzip ${savedir}/All_Polarised_moderate.vcf
bcftools index --tbi ${savedir}/All_Polarised_moderate.vcf.gz

java -jar /data/biosoftware/SnpEff/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'LOW')" ${datadir}/All_Polarised_annotated_Ensemble.vcf.gz > ${savedir}/All_Polarised_low.vcf

bgzip ${savedir}/All_Polarised_low.vcf
bcftools index --tbi ${savedir}/All_Polarised_low.vcf.gz

java -jar /data/biosoftware/SnpEff/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'MODIFIER')" ${datadir}/All_Polarised_annotated_Ensemble.vcf.gz > ${savedir}/All_Polarised_modifier.vcf

bgzip ${savedir}/All_Polarised_modifier.vcf
bcftools index --tbi ${savedir}/All_Polarised_modifier.vcf.gz

echo "VCF subsets done"

# Extract Intergenic regions
java -jar /data/biosoftware/SnpEff/snpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'intergenic_region')" ${datadir}/All_Polarised_annotated_Ensemble.vcf.gz > ${savedir}/All_Polarised_intergenic.vcf

bgzip ${savedir}/All_Polarised_intergenic.vcf
bcftools index --tbi ${savedir}/All_Polarised_intergenic.vcf.gz
