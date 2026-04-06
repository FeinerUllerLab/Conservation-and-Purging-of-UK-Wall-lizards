#! /bin/bash -l

# Objective: Trimming and preparing for mapping

module load trimmomatic/0.39

# Define the source and target directories
src_dir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/Rawdata"
mkdir /cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling
target_dir="/cfs/klemming/projects/snic/snic2022-23-124/Santiago/SNP_calling"

# File containing the list of samples to be used
dir_list_file="Samples_Used2"

# Loop over the directories listed in the file
while IFS= read -r dir_name; do
  # Construct the full path to the source directory
  dir="$src_dir/$dir_name"

  # Check if the directory exists before processing
  if [ -d "$dir" ]; then
    # Create a new folder with the same name in the target directory and also create subdirectories
    mkdir -p "$target_dir/$dir_name/trimReads"

    # Change into the source subdirectory
    cd "$dir" || { echo "Failed to cd into $dir"; continue; }

    # Count the number of files ending with *fq.gz*
    fq_files=(*fq.gz)
    fq_count=${#fq_files[@]}

    if [[ $fq_count -eq 2 ]]; then
      # If the number of fq.gz files is 2, use these files in trimmomatic
      input_files=("${fq_files[@]}")
    else
      # If not, use the file names containing *combined* in trimmomatic
      input_files=(*combined*)
    fi

    # Check if input_files array is non-empty
    if [[ ${#input_files[@]} -eq 0 ]]; then
      echo "No suitable input files found in $dir, skipping."
      continue
    fi

    # Run Trimmomatic on the selected files
    TRIMMOMATIC_ROOT=/pdc/software/eb/software/trimmomatic/0.39
    trimmomatic_jar="$TRIMMOMATIC_ROOT/trimmomatic.jar"
    adapter_file="$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa"

    java -jar "$trimmomatic_jar" PE -phred33 \
    "${input_files[0]}" "${input_files[1]}" \
    "$target_dir/$dir_name/trimReads/${dir_name}_trim_1.fq.gz" "$target_dir/$dir_name/trimReads/${dir_name}_trim_unpair1.fq.gz" \
    "$target_dir/$dir_name/trimReads/${dir_name}_trim_2.fq.gz" "$target_dir/$dir_name/trimReads/${dir_name}_trim_unpair2.fq.gz" \
    ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70

    # Print a message that the job is done and the directory name
    echo "Job done for directory: $dir_name"
    date

    # Return to the parent directory
    cd "$src_dir" || { echo "Failed to return to $src_dir"; exit 1; }
  else
    echo "Directory $dir does not exist, skipping."
  fi
done < "$dir_list_file"
