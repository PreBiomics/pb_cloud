#!/bin/bash
sample_name=$1
mkdir -p /output/tmp/mpa411/tmp
zcat /data/preprocessing/${sample_name}/*.fastq.gz | metaphlan --input_type fastq --tmp_dir /output/mpa411/tmp --index mpa_vJan21_TOY_CHOCOPhlAnSGB_202103 --bowtie2db /databases/mpa_database -o /output/mpa411/${sample_name}/${sample_name}_profile.tsv --force --no_map
sed -i -e 's/mpa_vJan21_TOY_CHOCOPhlAnSGB_202103/mpa_vJun23_CHOCOPhlAnSGB_202307/g' /data/mpa411/${sample_name}/${sample_name}_profile.tsv
