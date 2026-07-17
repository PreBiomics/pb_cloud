#!/bin/bash
set -euo pipefail

export PATH=/output/pb_cloud/preprocessing:$PATH

sample_name=${1}
nproc=${2}
host=${3}
mapper=${4}
qc=${5}

bowtie2 -h
SplitUnmappedSAM.py -h
