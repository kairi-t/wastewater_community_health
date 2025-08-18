# Programmed by: Kairi Tanaka
# Programmed on: 03-11-2024
# Programmed to: run fastp on the new sewage project data
# Modified by: Kairi Tanaka
# Modified on: 05-30-2024
# Modified to: adjust the code for less strict quality filtering and turn off merging

#!bin/bash/

# Set the forward and reverse extensions
fwd_ext=".1.fq.gz"
rev_ext=".2.fq.gz"

echo "Preparetion complete..."
echo "Initiating fastp for-loop"

# Loop through files
for file in $(ls *.1.fq.gz); do

  filename=$(echo "$file" | cut -d '.' -f 1)

  echo "Working on "
  echo $filename
  echo "..."

  # Extract sample name from filename
  fwd_sample="$filename$fwd_ext"
  echo $fwd_sample

  rev_sample="$filename$rev_ext"
  echo $rev_sample

  # Check if reverse read file exists
  if [ ! -f "$rev_sample" ]; then
    echo "Error: Reverse file not found: $rev_sample"
    continue
  fi

  # Run fastp
  fastp -i "$fwd_sample" -I "$rev_sample"\
        -o "clean_16s_fastq/${filename}_clean_1.fq.gz"\
        -O "clean_16s_fastq/${filename}_clean_2.fq.gz"\
        --detect_adapter_for_pe\
        --trim_poly_g --poly_g_min_len 5\
        --trim_poly_x --poly_x_min_len 5\
        --verbose\
        --overrepresentation_analysis\
        -R "${filename}_fastp_report"\
        -h "fastp_reports_16s/${filename}_fastp_report.html"\
        -w 10
done

echo "Running fastp completed!"
