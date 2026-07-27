#!/bin/bash
set -euo pipefail

export PATH=/output/pb_cloud/wgs:$PATH

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
SHARDS=$((THREADS * 2))

mkdir -p "${QC}"

BAM="${OUTDIR}/${SAMPLE}.bam"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
RECAL_BAM="${OUTDIR}/${SAMPLE}.recal.bam"

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

FLAGSTAT_THREADS=1
MOSDEPTH_THREADS=$(( THREADS - FLAGSTAT_THREADS ))
(( MOSDEPTH_THREADS < 1 )) && MOSDEPTH_THREADS=1

samtools flagstat \
    -@ "${FLAGSTAT_THREADS}" \
    "${SORTED_BAM}" \
    > "${QC}/${SAMPLE}.flagstat.txt" &

mosdepth \
    --threads "${MOSDEPTH_THREADS}" \
    --no-per-base \
    "${QC}/${SAMPLE}" \
    "${SORTED_BAM}" &

wait

grep "mapped (" "${QC}/${SAMPLE}.flagstat.txt" \
| sed -E 's/.*\(([0-9.]+)%.*/\1/' \
> "${QC}/${SAMPLE}.mapping_rate.txt"

awk '
$1=="total"{
    printf("%.2f\n",$4)
}
' "${QC}/${SAMPLE}.mosdepth.summary.txt" \
> "${QC}/${SAMPLE}.mean_depth.txt"

cp \
"${QC}/${SAMPLE}.mean_depth.txt" \
"${QC}/${SAMPLE}.sequencing_depth.txt"
MEAN=$(cat "${QC}/${SAMPLE}.mean_depth.txt")
awk -v mean="${MEAN}" '
BEGIN{
    threshold = 0.2 * mean
}
$1 >= threshold{
    pct = $2
}
END{
    printf("%.2f\n", pct * 100)
}
' "${QC}/${SAMPLE}.mosdepth.global.dist.txt" \
> "${QC}/${SAMPLE}.coverage_uniformity.txt"

rm -f "${QC}/${SAMPLE}.mosdepth.global.dist.txt"
rm -f "${QC}/${SAMPLE}.mosdepth.region.dist.txt"
rm -f "${QC}/${SAMPLE}.regions.bed.gz"
rm -f "${QC}/${SAMPLE}.regions.bed.gz.csi"

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
    -O $RECAL_BAM
samtools index $RECAL_BAM

rm "${SORTED_BAM}"
rm "${SORTED_BAM}.bai"

echo "[$(date)] HaplotypeCaller"

echo
echo "[$(date)] Preparing intervals"

SHARD_DIR="${OUTDIR}/${2}_shards"
mkdir -p "${SHARD_DIR}"

if [[ ! -f "${SHARD_DIR}/shard_31.list" ]]; then
    rm -f "${SHARD_DIR}"/*.list
    TOTAL=$(awk '
    {
        sum += $2
    }
    END{
        printf "%.0f\n", sum
    }
    ' "${REFERENCE}.fa.fai")
    TARGET=$((TOTAL / SHARDS))
    SHARD=0
    CURRENT=0
    while read -r CHR LEN REST
    do
        START=1
        while (( START <= LEN ))
        do
            REMAIN=$((TARGET - CURRENT))
            LEFT=$((LEN - START + 1))
            # Last shard gets everything that's left
            if (( SHARD == SHARDS - 1 )); then
                echo "${CHR}:${START}-${LEN}" \
                    >> "${SHARD_DIR}/shard_$(printf "%02d" "${SHARD}").list"
                break
            fi
            if (( LEFT <= REMAIN )); then
                echo "${CHR}:${START}-${LEN}" \
                    >> "${SHARD_DIR}/shard_$(printf "%02d" "${SHARD}").list"
                CURRENT=$((CURRENT + LEFT))
                START=$((LEN + 1))
            else
                END=$((START + REMAIN - 1))
                echo "${CHR}:${START}-${END}" \
                    >> "${SHARD_DIR}/shard_$(printf "%02d" "${SHARD}").list"
                START=$((END + 1))
                SHARD=$((SHARD + 1))
                CURRENT=0
            fi
        done
    done < "${REFERENCE}.fa.fai"
fi

echo
echo "[$(date)] Running HaplotypeCaller"

GVCF_DIR="${OUTDIR}/gvcf_parts"
mkdir -p "${GVCF_DIR}"

parallel \
    -j ${THREADS} \
    --halt soon,fail=1 \
    hc_worker.sh \
        "${REFERENCE}" \
        "${RECAL_BAM}" \
        {} \
        "${GVCF_DIR}"/{/.}.g.vcf.gz \
    ::: "${SHARD_DIR}"/*.list

rm -rf "${SHARD_DIR}"

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

gatk IndexFeatureFile \
    -I "${OUTDIR}/${SAMPLE}.vcf.gz"

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

bcftools index -n \
    "${OUTDIR}/${SAMPLE}.SNP.vcf.gz" \
    > "${QC}/${SAMPLE}.num_snps.txt"

bcftools index -n \
    "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz" \
    > "${QC}/${SAMPLE}.num_indels.txt"

rm "${OUTDIR}/${SAMPLE}.recal.table"
rm "${OUTDIR}/${SAMPLE}.recal.bam.bai"
rm "${OUTDIR}/${SAMPLE}.recal.bai"
rm "${OUTDIR}/${SAMPLE}.SNP.vcf.gz.tbi"
rm "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz.tbi"
rm "${OUTDIR}/${SAMPLE}.vcf.gz.tbi"

echo ""
echo "Pipeline completed successfully."
echo ""
echo "Output:"
echo "  QC metrics : ${QC}"
echo "  Recal BAM  : ${OUTDIR}/${SAMPLE}.recal.bam"
echo "  gVCF       : ${OUTDIR}/${SAMPLE}.vcf.gz"
echo "  SNP gVCF   : ${OUTDIR}/${SAMPLE}.SNP.vcf.gz"
echo "  INDEL gVCF : ${OUTDIR}/${SAMPLE}.INDEL.vcf.gz"
