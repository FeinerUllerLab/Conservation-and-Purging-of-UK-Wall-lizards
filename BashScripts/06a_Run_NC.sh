#!/bin/bash -l

# Objective: Run the nucletide count scrip on the filtired VCF that contains ingroup and outgroups.
    # The output are est-sfs compatible nucleotide counts per chromosome 

# It uses a VCF with both ingroup and outgroup samples to get organized nucleotide counts data for est-sfs. 

module load python/3.12.3
module load PDCOLD/23.12
module load pysam/0.22.1-cpeGNU-23.12

scriptdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts"
inputdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Pol_VCFs"
outdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr"
mkdir -p ${outdir}

date 

# Run the script - argument order matters 
python3 "${scriptdir}/13_Nucliotide_count.py" \
    "${inputdir}/Pre_polarization.vcf.gz" \
    "ALL" \
    "${outdir}" \
    "${scriptdir}/Samples_Used" \
    "${scriptdir}/Polarization_Out1" \
    "${scriptdir}/Polarization_Out2" \
    "${scriptdir}/Polarization_Out3"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Python script failed" >&2
    exit 1
fi

echo "Successfully generated chromosome count files in ${outdir}"
date
