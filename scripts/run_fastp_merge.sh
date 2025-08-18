conda init
module load conda

conda activate /usr/local/usrapps/methanomics/ktanaka2/bioinfo
fwd_ext="_R1.fastq.bz2"
rev_ext="_R2.fastq.bz2"

echo "Preparation complete..."
echo "Running fastp..."

for file in $(ls *_1.fastq.bz2); do
  filename=$(echo "$file" | cut -d '_' -f 1,2)

  echo "Working on "
  echo $filename
  echo "..."

  fwd_sample="$filename$fwd_ext"
  echo $fwd_sample

  rev_sample="$filename$rev_ext"
  echo $rev_sample

  if [ ! -f "$rev_sample"]; then
    echo "Error: Reverse file not found: $rev_file"
    continue
  fi

  fastp -i "$fwd_sample" -I "$rev_sample"\
	-o "${filename}_clean_R1.fastq.bz2" -O "${filename}_clean_R2.fastq.bz2"\
	--detect_adapter_for_pe\
	--trim_poly_g --poly_g_min_len 7\
	--trim_poly_x --poly_x_min_len 7\
	--verbose\
	--overrepresentation_analysis\
	-R "${filename}_fastp_report" -h "${filename}_fastp_report.html"\
	-w 12

  echo "fastp completed!"
done
