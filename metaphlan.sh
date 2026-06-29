#!/bin/bash
sample_name=$1
mkdir -p /output/tmp/mpa411/tmp
echo "Getting mpa database content"
ls /databases/mpa_database
zcat /input/${sample_name}/*.fastq.gz | metaphlan --input_type fastq --tmp_dir /output/tmp/mpa411/tmp --index mpa_vJan21_TOY_CHOCOPhlAnSGB_202103 --bowtie2db /databases/mpa_database -o /output/tmp/mpa411/${sample_name}/${sample_name}_profile.tsv --force --no_map
sed -i -e 's/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103/mpa_vJun23_CHOCOPhlAnSGB_202307/g' /output/tmp/mpa411/${sample_name}/${sample_name}_profile.tsv
