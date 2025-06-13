#!/bin/bash
sample_name=$1
mkdir -p /data/hnn39/${sample_name}
mkdir -p /rubineto/hnn39_${sample_name}

bzcat /data/preprocessing/${sample_name}/*.fastq.bz2 > /rubineto/hnn39_${sample_name}/${sample_name}.fastq
humann --input /rubineto/hnn39_${sample_name}/${sample_name}.fastq --output /rubineto/hnn39_${sample_name} --input-format fastq --nucleotide-database /databases/hnn_database/chocophlan --protein-database /databases/hnn_database/uniref --search-mode uniref90 --verbose --remove-temp-output --bypass-translated-search --taxonomic-profile /data/mpa411/${sample_name}/${sample_name}_profile.tsv
rm /rubineto/hnn39_${sample_name}/${sample_name}.fastq
mv /rubineto/hnn39_${sample_name}/${sample_name}_genefamilies.tsv /data/hnn39/${sample_name}/${sample_name}_genefamilies.tsv
mv /rubineto/hnn39_${sample_name}/${sample_name}_pathabundance.tsv /data/hnn39/${sample_name}/${sample_name}_pathabundance.tsv
mv /rubineto/hnn39_${sample_name}/${sample_name}_pathcoverage.tsv /data/hnn39/${sample_name}/${sample_name}_pathcoverage.tsv
rm -rf rm /rubineto/hnn39_${sample_name}/
