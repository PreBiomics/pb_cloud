#!/bin/bash

set -euo pipefail

if [ "$#" -lt 6 ]; then
    echo "Usage:"
    echo "$0 <R1.fastq.gz> <R2.fastq.gz> <reference.fa> <sample_name> <output_dir> <threads>"
    exit 1
fi

R1=$1
R2=$2
REFERENCE=$3
SAMPLE=$4
OUTDIR=$5
THREADS=${6}

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


echo "[$(date)] Running GATK HaplotypeCaller..."

gatk HaplotypeCaller \
    -R "${REFERENCE}" \
    -I "${SORTED_BAM}" \
    -O "${GVCF}" \
    -ERC GVCF


echo "[$(date)] Converting GVCF to VCF..."

gatk GenotypeGVCFs \
    -R "${REFERENCE}" \
    -V "${GVCF}" \
    -O "${RAW_VCF}"


echo "[$(date)] Filtering variants..."

bcftools filter \
    -i 'QUAL>=30 && DP>=10' \
    "${RAW_VCF}" \
    -Oz \
    -o "${FILTERED_VCF}"

echo "[$(date)] Indexing final VCF..."

bcftools index \
    "${FILTERED_VCF}"


echo "===================================="
echo "Pipeline finished"
echo "Final VCF:"
echo "${FILTERED_VCF}"
echo "===================================="
