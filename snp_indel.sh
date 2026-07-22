#!/bin/bash

set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage:"
    echo "$0 <sample_name> <reference_name> <threads>"
    exit 1
fi

SAMPLE=$1
R1=/input/${SAMPLE}/${SAMPLE}_R1.fastq.gz
R2=/input/${SAMPLE}/${SAMPLE}_R2.fastq.gz
REFERENCE=/databases/$2
OUTDIR=/output/tmp/wgs/${SAMPLE}
THREADS=$3

mkdir -p "${OUTDIR}"

SAM="${OUTDIR}/${SAMPLE}.sam"
BAM="${OUTDIR}/${SAMPLE}.bam"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"

GVCF="${OUTDIR}/${SAMPLE}.g.vcf.gz"
RAW_VCF="${OUTDIR}/${SAMPLE}.raw.vcf.gz"
FILTERED_VCF="${OUTDIR}/${SAMPLE}.filtered.vcf.gz"

echo "[$(date)] Running BWA-MEM2..."

bwa-mem2 mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" \
    > "${SAM}"

echo "[$(date)] Converting SAM to BAM..."

samtools view \
    -@ ${THREADS} \
    -bS "${SAM}" \
    -o "${BAM}"

rm "${SAM}"

echo "[$(date)] Sorting BAM..."

samtools sort \
    -@ ${THREADS} \
    -o "${SORTED_BAM}" \
    "${BAM}"

rm "${BAM}"


echo "[$(date)] Indexing BAM..."

samtools index \
    -@ ${THREADS} \
    "${SORTED_BAM}"
