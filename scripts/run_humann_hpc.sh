#BSUB -n 24
#!/bin/bash
#BSUB -J humann_run
#BSUB -W 80:00
#BSUB -R span[hosts=1]
#BSUB -o stdout.%J
#BSUB -e stderr.%J

source ~/.bashrc

conda init
module load conda
conda activate /usr/local/usrapps/methanomics/ktanaka2/humann

mkdir sams
mkdir bowtie2
mkdir profiles

for f in *_merged.fq.gz; do
  echo "Running MetaPhlAn on ${f}"

  bn=$(echo "${f}" | cut -d '_' -f 1,2 )

  metaphlan ${f} --input_type fastq\
	    -t rel_ab_w_read_stats\
            -s sams/${bn}.sam.bz2\
  	    --bowtie2out bowtie2/${bn}.bowtie2.bz2\
  	    -o profiles/${bn}_profiled.tsv\
  	    --offline -t rel_ab_w_read_stats\
  	    --bowtie2db /rs1/researchers/z/zlyu3/biobakery_databases/metaphlan_databases/\
  	    --index mpa_vJun23_CHOCOPhlAnSGB_202307\
  	    --nproc 24;

  humann --input ${f}\
	 --output humann_out\
	 --threads 24 --memory-use maximum\
	 --taxonomic-profile profiles/${bn}_profiled.tsv\
	 --nucleotide-database /rs1/researchers/z/zlyu3/biobakery_databases/humann_databases/data/chocophlan\
	 --protein-database    /rs1/researchers/z/zlyu3/biobakery_databases/humann_databases/data/uniref\
	 --pathways-database   /rs1/researchers/z/zlyu3/biobakery_databases/humann_databases/data/pathways/metacyc_pathways\
	 --verbose\
	 --metaphlan-options "-t rel_ab_w_read_stats --offline --bowtie2db /rs1/researchers/z/zlyu3/biobakery_databases/metaphlan_databases/ --index mpa_vJun23_CHOCOPhlAnSGB_202307";

done
