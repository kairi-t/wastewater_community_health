# Programmed by: Kairi Tanaka
# Programmed on: 10-05-2024
# Programmed to: change the compression from bz2 to gz
#!/bin/bash

# Loop over all directories starting with 'sample'
for sample_dir in Li_*; do
    if [ -d "$sample_dir" ]; then
        echo "Processing $sample_dir"

        # Define the paths to the forward and reverse read files
        fwd_file="$sample_dir/${sample_dir}.1.fq.bz2"
        rev_file="$sample_dir/${sample_dir}.2.fq.bz2"

        # Sanity check: Ensure both files exist
        if [ ! -f "$fwd_file" ]; then
            echo "Error: Forward read file $fwd_file does not exist!"
            continue
        fi
        if [ ! -f "$rev_file" ]; then
            echo "Error: Reverse read file $rev_file does not exist!"
            continue
        fi

        # Convert the forward read file from bz2 to gz
        echo "Converting $fwd_file to gz"
        bzcat "$fwd_file" | gzip -c > "${sample_dir}/${sample_dir}.1.fq.gz"

        # Convert the reverse read file from bz2 to gz
        echo "Converting $rev_file to gz"
        bzcat "$rev_file" | gzip -c > "${sample_dir}/${sample_dir}.2.fq.gz"

    else
        echo "Skipping $sample_dir as it is not a directory"
    fi
done
