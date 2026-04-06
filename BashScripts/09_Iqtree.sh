#!/bin/bash -l

# Objective: Perform an Iqtree analysis for the only-variant data
    #OBS: It requres two independent runs, One to generate a varsites file and one to run the actrual analysys

# Modules
module load PDC/23.12
module load iqtree
module load python

# Get vcf2phylip
git clone https://github.com/edgardomortiz/vcf2phylip.git
Path=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts/vcf2phylip

# Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling/FinalVCFs/OnlySNPs_Data/Final_dataset
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Tree 
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Tree 

#Transform VCF to PHYLIP format. 
python3 $Path/vcf2phylip.py -i ${datadir}/All_Orifgins_final.vcf.gz --output-folder ${savedir}

#Run Iqtree to create Varsites file: 
iqtree2 -s ${savedir}/All_Origins_SNP_final.min4.phy -m GTR+ASC -st DNA -nt AUTO

# Run after (different job)
# Run iqtree with bootstraping  
iqtree2 -s ${savedir}/All_Origins_SNP_final.min4.phy.varsites.phy -m GTR+ASC -st DNA -nt AUTO -cptime 400 -B 1000

#End 
