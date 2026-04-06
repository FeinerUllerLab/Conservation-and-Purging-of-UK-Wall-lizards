#!/bin/bash -l

# Objective: Evaluate population structure for K 1 to 10 uses the .bed only-snps data generated in step 06

# Modules
module load  bioinfo-tools
module load ADMIXTURE/1.3.0

#Define directories
datadir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen
mkdir -p /cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Admixture 
savedir=/cfs/klemming/projects/snic/snic2022-23-124/Santiago/PopGen/Admixture 

#Run ADMIXTURE for 1 to 10 ks
for K in {1..10}; do
    echo "Running ADMIXTURE for K=${K}"
    admixture --cv -j8 ${datadir}/All_Origins_SNP_final_plink.bed $K | tee ${savedir}/Admixture_K${K}.out
done

# Gather CV errors into a summary file
grep -h "CV error" ${savedir}/Admixture_K*.out > ${savedir}/CV_errors_summary.txt
echo "Cross-validation completed. Results saved in ${savedir}/CV_errors_summary.txt"

#End 
