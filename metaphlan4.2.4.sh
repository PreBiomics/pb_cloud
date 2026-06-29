#!/bin/bash
sample_name=$1
mkdir -p /output/tmp/mpa422/tmp
echo "Getting mpa database content"
ls /databases/
zcat /input/${sample_name}/*.fastq.gz | metaphlan --input_type fastq --tmp_dir /output/tmp/mpa422/tmp --index mpa_vJan25_CHOCOPhlAnSGB_202503 --bowtie2db /databases/ -o /output/tmp/mpa422/${sample_name}/${sample_name}_profile.tsv --force --no_map
sed -i "2d" /output/tmp/mpa422/${sample_name}/${sample_name}_profile.tsv
