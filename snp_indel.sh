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
QC="${OUTDIR}/QC"
THREADS=$3
SCATTER_COUNT=64

mkdir -p "${OUTDIR}"
mkdir -p "${QC}"

BAM="${OUTDIR}/${SAMPLE}.bam"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"


echo "[$(date)] Running BWA-MEM2..."

bwa-mem2 mem \
    -t ${THREADS} \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" | samtools view \
    -@ "${THREADS}" \
    -b \
    -o "${BAM}" -

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

echo "[$(date)] Collecting QC metrics"

# Alignment statistics
samtools flagstat "${SORTED_BAM}" \
    > "${QC}/${SAMPLE}.flagstat.txt"

samtools stats "${SORTED_BAM}" \
    > "${QC}/${SAMPLE}.stats.txt"

# Coverage statistics
samtools coverage "${SORTED_BAM}" \
    > "${QC}/${SAMPLE}.coverage.txt"

# Mean depth
samtools depth -a "${SORTED_BAM}" \
    | awk '{sum+=$3;cnt++} END {print sum/cnt}' \
    > "${QC}/${SAMPLE}.mean_depth.txt"

echo "[$(date)] BaseRecalibrator"

gatk BaseRecalibrator \
    -R "${REFERENCE}.fa" \
    -I "${SORTED_BAM}" \
    --known-sites /databases/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-sites /databases/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites /databases/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O "${OUTDIR}/${SAMPLE}.recal.table"

echo "[$(date)] ApplyBQSR"

gatk ApplyBQSR \
    -R "${REFERENCE}.fa" \
    -I "${SORTED_BAM}" \
    --bqsr-recal-file "${OUTDIR}/${SAMPLE}.recal.table" \
    -O "${OUTDIR}/${SAMPLE}.recal.bam"

samtools index "${OUTDIR}/${SAMPLE}.recal.bam"

rm "${SORTED_BAM}"
rm "${SORTED_BAM}.bai"

echo "[$(date)] HaplotypeCaller"

##############################################################################
# HaplotypeCaller (parallelized)
###############################################################################

echo "[$(date)] Preparing intervals"

INTERVAL_DIR="${OUTDIR}/${2}_intervals_${SCATTER_COUNT}"

if [ ! -d "${INTERVAL_DIR}" ]; then
    mkdir -p "${INTERVAL_DIR}"

    gatk SplitIntervals \
        -R "${REFERENCE}.fa" \
        -L "${REFERENCE}.fa" \
        --scatter-count ${SCATTER_COUNT} \
        -O "${INTERVAL_DIR}"
fi

GVCF_DIR="${OUTDIR}/gvcf_parts"
mkdir -p "${GVCF_DIR}"

echo "[$(date)] Running HaplotypeCaller on ${SCATTER_COUNT} intervals"

find "${INTERVAL_DIR}" -name "*.interval_list" | sort | \
parallel -j ${THREADS} \
'
INTERVAL={}
NAME=$(basename ${INTERVAL} .interval_list)

gatk HaplotypeCaller \
    -R '"${REFERENCE}.fa"' \
    -I '"${OUTDIR}/${SAMPLE}.recal.bam"' \
    -L ${INTERVAL} \
    -ERC GVCF \
    --native-pair-hmm-threads 1 \
    -O '"${GVCF_DIR}"'/${NAME}.g.vcf.gz
'

echo "[$(date)] Gathering GVCFs"

INPUTS=()

while IFS= read -r file
do
    INPUTS+=("-I" "$file")
done < <(find "${GVCF_DIR}" -name "*.g.vcf.gz" | sort)

gatk GatherVcfs \
    "${INPUTS[@]}" \
    -O "${OUTDIR}/${SAMPLE}.vcf.gz"

rm -rf "${GVCF_DIR}"

###############################################################################
# Variant selection
###############################################################################

echo "[$(date)] Select SNP and INDELs"

gatk SelectVariants \
        -R "${REFERENCE}.fa" \
        -V ${OUTDIR}/${SAMPLE}.vcf.gz \
        --select-type-to-include SNP \
        -O "${OUTDIR}/${SAMPLE}.SNP.vcf.gz"

gatk SelectVariants \
        -R "${REFERENCE}.fa" \
        -V ${OUTDIR}/${SAMPLE}.vcf.gz \
        --select-type-to-include INDEL \
        -O "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz"

rm "${OUTDIR}/${SAMPLE}.recal.table"
rm "${OUTDIR}/${SAMPLE}.recal.bam.bai"
rm "${OUTDIR}/${SAMPLE}.recal.bai"
rm "${OUTDIR}/${SAMPLE}.SNP.vcf.gz.tbi"
rm "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz.tbi"
rm "${OUTDIR}/${SAMPLE}.vcf.gz.tbi"
rm -rf "${OUTDIR}/${2}_intervals_*"

echo ""
echo "Pipeline completed successfully."
echo ""
echo "Output:"
echo "  QC metrics : ${QC}"
echo "  Recal BAM  : ${OUTDIR}/${SAMPLE}.recal.bam"
echo "  gVCF       : ${OUTDIR}/${SAMPLE}.vcf.gz"
echo "  SNP gVCF   : ${OUTDIR}/${SAMPLE}.SNP.vcf.gz"
echo "  INDEL gVCF : ${OUTDIR}/${SAMPLE}.INDEL.vcf.gz"
