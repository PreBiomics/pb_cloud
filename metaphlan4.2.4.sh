#!/bin/bash
sample_name=$1
nproc=$2
input=$3
mkdir -p /data/mpa424/${sample_name}
bzcat /data/${input}/${sample_name}/*.fastq.gz | metaphlan --input_type fastq --tmp_dir /tmp --index mpa_vJan25_CHOCOPhlAnSGB_202503 --bowtie2db /databases/mpa_database -o /data/mpa424/${sample_name}/${sample_name}_profile.tsv --force --no_map --skip_unclassified_estimation --nproc $nproc
