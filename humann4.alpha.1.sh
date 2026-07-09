#!/bin/bash
set -euo pipefail
sample_name=${1}
nproc=${2}

mkdir -p /output/tmp/hnn4/tmp/${sample_name}/

zcat /input/${sample_name}/*.fastq.gz > /output/tmp/hnn4/tmp/${sample_name}/${sample_name}.fastq

humann --input /output/tmp/hnn4/tmp/${sample_name}/${sample_name}.fastq \
	--output /output/tmp/hnn4/${sample_name} --input-format fastq \
	--nucleotide-database /databases/chocophlan --protein-database /databases/uniref \
	--threads ${nproc} --verbose --remove-temp-output \
	--metaphlan-options --index mpa_vOct22_CHOCOPhlAnSGB_202403 --offline \
	--bowtie2db /databases/mpa_database -t rel_ab_w_read_stats \
	--utility-database /databases/utility_mapping --remove-temp-output
  
humann_renorm_table -i /output/tmp/hnn4/${sample_name}/${sample_name}_2_genefamilies.tsv -u relab -o /output/tmp/hnn4/${sample_name}/${sample_name}_2_genefamilies_relab.tsv
humann_renorm_table -i /output/tmp/hnn4/${sample_name}/${sample_name}_4_pathabundance.tsv -u relab -o /output/tmp/hnn4/${sample_name}/${sample_name}_4_pathabundance_relab.tsv

rm -rf /output/tmp/hnn4/tmp/${sample_name}/
