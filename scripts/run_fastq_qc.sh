# Programmed by: Kairi Tanaka
# Programmed on: 10-05-2024
# Programmed to: change the compression from bz2 to gz
#!/bin/bash

# Activate conda env containing fastQC and fastp
source ~/.condarc
conda init
conda activate bioinfo

# Loop over all directories starting with 'sample'
for sample_name in Li_*; do
    if [ -d "$sample_name" ]; then
        echo "Processing $sample_name"

        # Define the paths to the forward and reverse read files
        fwd_file="$sample_name/${sample_name}.1.fq.gz"
        rev_file="$sample_name/${sample_name}.2.fq.gz"

        # Sanity check: Ensure both files exist
        if [ ! -f "$fwd_file" ]; then
            echo "Error: Forward read file $fwd_file does not exist!"
            continue
        fi
        if [ ! -f "$rev_file" ]; then
            echo "Error: Reverse read file $rev_file does not exist!"
            continue
        fi

        # Run FastQC
        echo "Running QC using FastQC"
        fastqc "$fwd_file" "$rev_file"

        # Create folder to save FastQC results
        mkdir "$sample_name/${sample_name}_fastqc/"
        mv "$sample_name/${sample_name}.?_fastqc.*" "$sample_name/${sample_name}_fastqc/"

        # Run fastp
        echo "Running fastp"
        fastp -i "$fwd_file" -I "$rev_file"\
        -o "${sample_name}/${sample_name}_clean_1.fq.gz" -O "${sample_name}/${sample_name}_clean_2.fq.gz"\
        --detect_adapter_for_pe\
        -R "${sample_name}/${sample_name}_fastp_report" -h "${sample_name}/${sample_name}_fastp_report.html"
        -w 8 # use 8 cores; change if necessary

        # Create folder to save fastp results
        mkdir "$sample_name/${sample_name}_fastp/"
        mv "$sample_name/*_fastp_*" "$sample_name/${sample_name}_fastp/"

    else
        echo "Skipping $sample_name as it is not a directory"
    fi
done

# Deactivate conda env
conda deactivate
