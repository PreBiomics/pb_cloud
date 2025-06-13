#!/bin/bash
sample_name=$1
kneaddata --input1 /data/demultiplexed_data/${sample_name}/${sample_name}_R1.fastq.gz --input2 /data/demultiplexed_data/${sample_name}/${sample_name}_R2.fastq.gz --output-prefix ${sample_name} -db /databases/knd_database/hg39 --output /data/preprocessing/${sample_name}
rm /data/preprocessing/${sample_name}/*_hg39_* /data/preprocessing/${sample_name}/*.trimmed.* /data/preprocessing/${sample_name}/*.repeats.* 
mv /data/preprocessing/${sample_name}/${sample_name}_paired_1.fastq /data/preprocessing/${sample_name}/${sample_name}_R1.fastq
mv /data/preprocessing/${sample_name}/${sample_name}_paired_2.fastq /data/preprocessing/${sample_name}/${sample_name}_R2.fastq
cat /data/preprocessing/${sample_name}/${sample_name}_unmatched_2.fastq >> /data/preprocessing/${sample_name}/${sample_name}_unmatched_1.fastq
mv /data/preprocessing/${sample_name}/${sample_name}_unmatched_1.fastq /data/preprocessing/${sample_name}/${sample_name}_UN.fastq
rm /data/preprocessing/${sample_name}/${sample_name}_unmatched_2.fastq
bzip2 /data/preprocessing/${sample_name}/${sample_name}_R1.fastq
bzip2 /data/preprocessing/${sample_name}/${sample_name}_R2.fastq
bzip2 /data/preprocessing/${sample_name}/${sample_name}_UN.fastq
