#!/bin/bash

# Objective: Run ests-sfs per chormosome using the kimura model. 
    # Other model options are available (JC - Rate 6)

# Environment 
mamba activate est-sfst #Before running

# Directories
scriptdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Scripts"
dir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation"
countdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr"
outdir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Polarisation/Per_chr/Kimura"
mkdir -p ${outdir}

# Path to chromosome list file
chromosome_list="${scriptdir}/Chromosome_autosomal.txt"

# Read the specific chromosome for this array task
chromosome=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$chromosome_list" | tr -d '\r')

# Input file for this chromosome
input_file="${countdir}/nucleotide_counts_${chromosome}.tsv"

# Check if the input file exists
if [[ ! -f "$input_file" ]]; then
    echo "ERROR: Input file not found: $input_file" >&2
    exit 1
fi

# Run est-sfs for this chromosome with chromosome-specific output names
echo "Processing chromosome $chromosome with input file $input_file"
est-sfs ${dir}/config-kimura.txt $input_file ${dir}/seedfile.txt \
    ${outdir}/PmuralisKimura-Out_${chromosome}.txt \
    ${outdir}/PmuralisKimura-pval_${chromosome}.txt

# Check if the command was successful
if [[ $? -ne 0 ]]; then
    echo "ERROR: est-sfs failed for chromosome $chromosome" >&2
    exit 1
fi

echo "Successfully processed chromosome $chromosome"