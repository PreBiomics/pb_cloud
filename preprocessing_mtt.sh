#!/bin/bash
set -euo pipefail

export PATH=/output/pb_cloud/preprocessing:$PATH

sample_name=${1}
nproc=${2}
host=${3}
mapper=${4}
qc=${5}

preprocessing.py \
	-m ${mapper} \
	-i /input/${sample_name}/ \
	-s ${sample_name} \
	-f R1 \
	-r R2 \
	-n ${nproc} \
	-t ${nproc} \
	${host} \
	-d /databases/kraken2_databases \
	-x /databases/bowtie2_indexes \
	-p \
	-o /output/tmp/preprocessing/${sample_name} \
	-q ${qc} \
	--verbose
