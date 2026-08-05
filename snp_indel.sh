#!/bin/bash
set -euo pipefail

#export PATH=/output/pb_cloud/wgs:$PATH
export PATH=/script/wgs:$PATH

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

FINAL_BAM="${OUTDIR}/${SAMPLE}.bam"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
RECAL_BAM="${OUTDIR}/${SAMPLE}.recal.bam"
TMP_BAM="${OUTDIR}/${SAMPLE}.tmp.bam"

R1_FASTP=/input/${SAMPLE}/${SAMPLE}_R1.fastp.fastq.gz
R2_FASTP=/input/${SAMPLE}/${SAMPLE}_R2.fastp.fastq.gz
FASTP_HTML="/input/${SAMPLE}/${SAMPLE}.html"
FASTP_JSON="/input/${SAMPLE}/${SAMPLE}.json"
FLAGSTAT="${QC}/${SAMPLE}.flagstat.txt"
MOSDEPTH_PREFIX="${QC}/${SAMPLE}"
REPORT="${OUTDIR}/${SAMPLE}.report.tsv"

echo
echo "[$(date)] Running FASTp..."

fastp -i "${R1}" -o "${R1_FASTP}" -I "${R2}" -O "${R2_FASTP}" -z 4 -e 20 -l 75 -w "${THREADS}" -q 20 \
    --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT \
    --trim_poly_g --trim_poly_x -y -n 2 -h "${FASTP_HTML}" -j "${FASTP_JSON}"

mv "${R1_FASTP}" "${R1}"
mv "${R2_FASTP}" "${R2}"
rm -f "${FASTP_HTML}"

echo
echo "[$(date)] Running BWA-MEM2..."

bwa-mem2 mem \
    -t "${THREADS}" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    "${REFERENCE}" \
    "${R1}" \
    "${R2}" | \
samtools view \
    -@ "${THREADS}" \
    -b - | \
samtools sort \
    -@ "${THREADS}" \
    -o "${SORTED_BAM}" \
    -

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

awk '
$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ {
    print $1 "\t0\t" $2
}
' "${REFERENCE}.fa.fai" > "${QC}/primary_genome.bed"

mosdepth \
  -b "${QC}/primary_genome.bed" \
  --no-per-base \
  --threads "${MOSDEPTH_THREADS}" \
  "${QC}/${SAMPLE}" \
  "${SORTED_BAM}" &

wait

RAW_READS=$(jq -r '.summary.before_filtering.total_reads' "${FASTP_JSON}")
QC_READS=$(jq -r '.summary.after_filtering.total_reads' "${FASTP_JSON}")
Q30_RATE=$(jq -r '
    (.summary.after_filtering.q30_rate * 100)
    | tostring
' "${FASTP_JSON}")

MAPPING_RATE=$(
awk '
/mapped \(/{
    sub(/^.*\(/,"")
    sub(/%.*/,"")
    print
    exit
}
' "${FLAGSTAT}"
)

MEAN_DEPTH=$(
awk '
$1=="total"{
    print $4
}
' "${MOSDEPTH_PREFIX}.mosdepth.summary.txt"
)

########################################
# Coverage uniformity
# (% genome covered at >=20% of mean depth)
########################################

THRESHOLD=$(awk -v d="${MEAN_DEPTH}" 'BEGIN{printf "%.0f", d*0.2}')

UNIFORMITY=$(
awk -v t="$THRESHOLD" '
$1=="total" {
    if ($2 >= t)
        val = $3 * 100
}
END {
    printf "%.2f\n", val
}' "${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt"
)

echo
echo "[$(date)] Preparing intervals"

SHARD_DIR="${OUTDIR}/${2}_shards"
mkdir -p "${SHARD_DIR}"

EXPECTED=$(printf "shard_%02d.list" $((SHARDS-1)))

if [[ ! -f "${SHARD_DIR}/${EXPECTED}" ]]
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
echo "[$(date)] BaseRecalibrator"

BQSR_DIR="${OUTDIR}/bqsr_parts"
BQSR_BAM_DIR="${OUTDIR}/recal_parts"
INTERVAL_DIR="${OUTDIR}/intervals"

mkdir -p "${BQSR_DIR}" "${BQSR_BAM_DIR}" "${INTERVAL_DIR}"

PRIMARY=(
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18
    chr19 chr20 chr21 chr22 chrX chrY chrM
)

JOBFILE=$(mktemp)

# Create one interval list per primary chromosome
for c in "${PRIMARY[@]}"; do
    LIST="${INTERVAL_DIR}/${c}.list"
    echo "${c}" > "${LIST}"
    echo "${LIST}" >> "${JOBFILE}"
done

# Create interval list for all remaining contigs
OTHER_LIST="${INTERVAL_DIR}/other.list"
> "${OTHER_LIST}"

while read -r CONTIG LEN REST
do
    if [[ ! " ${PRIMARY[*]} " =~ " ${CONTIG} " ]]; then
        echo "${CONTIG}" >> "${OTHER_LIST}"
    fi
done < "${REFERENCE}.fa.fai"

# Only add other.list if it contains anything
[[ -s "${OTHER_LIST}" ]] && echo "${OTHER_LIST}" >> "${JOBFILE}"

parallel -j "${THREADS}" '
LIST={}
NAME=$(basename "$LIST" .list)

gatk BaseRecalibrator \
    -R "'"${REFERENCE}.fa"'" \
    -I "'"${SORTED_BAM}"'" \
    --known-sites /databases/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --known-sites /databases/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites /databases/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -L "$LIST" \
    -O "'"${BQSR_DIR}"'/${NAME}.table"
' :::: "${JOBFILE}"

echo "[$(date)] GatherBQSRReports"

TABLES=()

while read -r LIST
do
    NAME=$(basename "${LIST}" .list)
    TABLE="${BQSR_DIR}/${NAME}.table"

    [[ -f "${TABLE}" ]] && TABLES+=(-I "${TABLE}")
done < "${JOBFILE}"

gatk GatherBQSRReports \
    "${TABLES[@]}" \
    -O "${OUTDIR}/${SAMPLE}.recal.table"

rm -rf "${BQSR_DIR}"

echo "[$(date)] ApplyBQSR"

parallel -j "${THREADS}" '
LIST={}
NAME=$(basename "$LIST" .list)

gatk ApplyBQSR \
    -R "'"${REFERENCE}.fa"'" \
    -I "'"${SORTED_BAM}"'" \
    --bqsr-recal-file "'"${OUTDIR}/${SAMPLE}.recal.table"'" \
    -L "$LIST" \
    -O "'"${BQSR_BAM_DIR}"'/${NAME}.bam"

' :::: "${JOBFILE}"

rm -f "${SORTED_BAM}"
rm -f "${SORTED_BAM}.bai"

echo "[$(date)] GatherBamFiles"

BAMS=()

while read -r LIST
do
    NAME=$(basename "${LIST}" .list)
    BAM="${BQSR_BAM_DIR}/${NAME}.bam"

    [[ -f "${BAM}" ]] && BAMS+=(-I "${BAM}")
done < "${JOBFILE}"

gatk GatherBamFiles \
    "${BAMS[@]}" \
    -O "${RECAL_BAM}"

rm -rf "${BQSR_BAM_DIR}"
rm -rf "${INTERVAL_DIR}"
rm -f "${JOBFILE}"

samtools index "${RECAL_BAM}"

echo "[$(date)] HaplotypeCaller"

echo
echo "[$(date)] Running HaplotypeCaller"

GVCF_DIR="${OUTDIR}/gvcf_parts"
VCF_DIR="${OUTDIR}/vcf_parts"
mkdir -p "${GVCF_DIR}"
mkdir -p "${VCF_DIR}"

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
    -O "${OUTDIR}/${SAMPLE}.g.vcf.gz"

echo "[$(date)] Genotyping shards"

find "${GVCF_DIR}" -name "shard_*.g.vcf.gz" \
| sort \
| parallel -j "${THREADS}" '
    BASE=$(basename {} .g.vcf.gz)
    gatk IndexFeatureFile \
        -I {}
    gatk GenotypeGVCFs \
        -R "'"${REFERENCE}.fa"'" \
        -V {} \
        -O "'"${VCF_DIR}"'/${BASE}.vcf.gz"
    gatk IndexFeatureFile \
        -I "'"${VCF_DIR}"'/${BASE}.vcf.gz"
'

rm -rf "${GVCF_DIR}"

echo "[$(date)] Gathering VCFs"

VCF_INPUTS=()

for VCF in $(find "${VCF_DIR}" -name "shard_*.vcf.gz" | sort)
do
    VCF_INPUTS+=("-I" "$VCF")
done

gatk GatherVcfs \
    "${VCF_INPUTS[@]}" \
    -O "${OUTDIR}/${SAMPLE}.vcf.gz"
    
rm -rf "${VCF_DIR}"

gatk IndexFeatureFile \
    -I "${OUTDIR}/${SAMPLE}.vcf.gz"

echo "[$(date)] Select SNP and INDELs"

gatk SelectVariants \
        -R "${REFERENCE}.fa" \
        -V ${OUTDIR}/${SAMPLE}.vcf.gz \
        --select-type-to-include SNP \
        -O "${OUTDIR}/${SAMPLE}.SNP.vcf.gz" &

gatk SelectVariants \
        -R "${REFERENCE}.fa" \
        -V ${OUTDIR}/${SAMPLE}.vcf.gz \
        --select-type-to-include INDEL \
        -O "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz" &
wait 

SNPS=$(bcftools view -H "${OUTDIR}/${SAMPLE}.SNP.vcf.gz" | wc -l)
INDELS=$(bcftools view -H "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz" | wc -l)

{
printf "Metric\tValue\n"
printf "Sample\t%s\n" "${SAMPLE}"
printf "Raw_reads\t%s\n" "${RAW_READS}"
printf "QC_passed_reads\t%s\n" "${QC_READS}"
printf "Q30_percent\t%.2f\n" "${Q30_RATE}"
printf "Mean_depth\t%.2f\n" "${MEAN_DEPTH}"
printf "Coverage_uniformity_percent\t%.2f\n" "${UNIFORMITY}"
printf "Mapping_rate_percent\t%.2f\n" "${MAPPING_RATE}"
printf "SNPs\t%s\n" "${SNPS}"
printf "INDELs\t%s\n" "${INDELS}"
} > "${REPORT}"

echo "QC report written to ${REPORT}"

rm -f "${OUTDIR}/${SAMPLE}.recal.table"
rm -f "${OUTDIR}/${SAMPLE}.recal.bam.bai"
rm -f "${OUTDIR}/${SAMPLE}.recal.bai"
rm -f "${OUTDIR}/${SAMPLE}.SNP.vcf.gz.tbi"
rm -f "${OUTDIR}/${SAMPLE}.INDEL.vcf.gz.tbi"
rm -f "${OUTDIR}/${SAMPLE}.vcf.gz"
rm -f "${OUTDIR}/${SAMPLE}.vcf.gz.tbi"
rm -rf ${QC}
mv ${RECAL_BAM} ${FINAL_BAM}

echo ""
echo "Pipeline completed successfully."
echo ""
echo "Outputs:"
echo "  QC metrics : ${REPORT}"
echo "  Recal BAM  : ${OUTDIR}/${SAMPLE}.bam"
echo "  gVCF       : ${OUTDIR}/${SAMPLE}.g.vcf.gz"
echo "  SNP gVCF   : ${OUTDIR}/${SAMPLE}.SNP.vcf.gz"
echo "  INDEL gVCF : ${OUTDIR}/${SAMPLE}.INDEL.vcf.gz"
