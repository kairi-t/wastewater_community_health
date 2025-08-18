i#!/bin/bash

# Initialize an array
r1_files=()

# Fill the array using a while loop and read
while IFS= read -r file; do
  r1_files+=("$file")
done < <(ls *_clean_R1.fastq.gz)

# Loop through each R1 file
for r1_file in "${r1_files[@]}"; do
  # Extract the sample ID by removing the read id part
  sampleID=${r1_file%_clean_R1.fastq.gz}

  # Define the R2 file name
  r2_file="${sampleID}_clean_R2.fastq.gz"

  # Check if corresponding R2 file exists
  if [[ -f "$r2_file" ]]; then
    # Concatenate R1 and R2 files into a merged file
    cat "$r1_file" "$r2_file" > "${sampleID}_merged.fq.gz"
  else
    echo "Missing R2 file for $sampleID, skipping..."
  fi
done

echo "Merging complete!"
