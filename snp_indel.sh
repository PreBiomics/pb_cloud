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
RECALL_BAM="${OUTDIR}/${SAMPLE}.recal.bam"

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
# samtools flagstat "${SORTED_BAM}" \
#     > "${QC}/${SAMPLE}.flagstat.txt"

# samtools stats "${SORTED_BAM}" \
#     > "${QC}/${SAMPLE}.stats.txt"

# Coverage statistics
# samtools coverage "${SORTED_BAM}" \
#     > "${QC}/${SAMPLE}.coverage.txt"

# Mean depth
# samtools depth -a "${SORTED_BAM}" \
#     | awk '{sum+=$3;cnt++} END {print sum/cnt}' \
#     > "${QC}/${SAMPLE}.mean_depth.txt"

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
    -O $RECALL_BAM
samtools index $RECALL_BAM

rm "${SORTED_BAM}"
rm "${SORTED_BAM}.bai"

echo "[$(date)] HaplotypeCaller"

################################################################################
# Interval generation (only once per reference)
################################################################################

echo
echo "[$(date)] Preparing intervals"

INTERVAL_DIR="${OUTDIR}/${2}_intervals"

if [[ ! -d "${INTERVAL_DIR}" ]]; then
    mkdir -p "${INTERVAL_DIR}"
    TARGET=100000000
    COUNT=0
    while read -r CHR LEN _; do
        START=1
        while (( START <= LEN )); do
            END=$((START + TARGET - 1))
            if (( END > LEN )); then
                END=$LEN
            fi
            printf "%s:%d-%d\n" \
                "$CHR" \
                "$START" \
                "$END" \
                > "${INTERVAL_DIR}/$(printf "%03d" ${COUNT}).interval_list"
            COUNT=$((COUNT + 1))
            START=$((END + 1))
        done
    done < "${REFERENCE}.fa.fai"
fi

################################################################################
# Parallel HaplotypeCaller
################################################################################

echo
echo "[$(date)] Running HaplotypeCaller"
GVCF_DIR="${OUTDIR}/gvcf_parts"
mkdir -p "${GVCF_DIR}"

export REFERENCE
export RECAL_BAM
export GVCF_DIR

parallel \
    --jobs "${THREADS}" \
    --halt soon,fail=1 \
'
FILE={}
NAME=$(basename "${FILE}" .interval_list)

gatk HaplotypeCaller \
    -R '"${REFERENCE}.fa"' \
    -I '"${RECAL_BAM}"' \
    -L "${FILE}" \
    -ERC GVCF \
    --native-pair-hmm-threads 1 \
    -O '"${GVCF_DIR}"'/"${NAME}".g.vcf.gz
' ::: "${INTERVAL_DIR}"/*.interval_list

################################################################################
# Gather GVCFs
################################################################################

echo
echo "[$(date)] Gathering GVCFs"

INPUTS=()
for FILE in "${GVCF_DIR}"/*.g.vcf.gz
do
    INPUTS+=("-I")
    INPUTS+=("${FILE}")
done

gatk GatherVcfs \
    "${INPUTS[@]}" \
    -O "${OUTDIR}/${SAMPLE}.vcf.gz"

rm -rf "${GVCF_DIR}"


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
rm "$INTERVAL_FILE"

echo ""
echo "Pipeline completed successfully."
echo ""
echo "Output:"
echo "  QC metrics : ${QC}"
echo "  Recal BAM  : ${OUTDIR}/${SAMPLE}.recal.bam"
echo "  gVCF       : ${OUTDIR}/${SAMPLE}.vcf.gz"
echo "  SNP gVCF   : ${OUTDIR}/${SAMPLE}.SNP.vcf.gz"
echo "  INDEL gVCF : ${OUTDIR}/${SAMPLE}.INDEL.vcf.gz"
