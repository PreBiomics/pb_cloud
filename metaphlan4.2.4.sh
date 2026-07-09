#!/bin/bash
set -euo pipefail
sample_name=\$1
nproc=\$2
index=mpa_vJan25_CHOCOPhlAnSGB_202503

mkdir -p "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp"
mkdir -p "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}"

echo "Running default profiling..."
zcat "/input/\${sample_name}"/*.fastq.gz | metaphlan \\
  --input_type fastq --skip_unclassified_estimation \\
  --tmp_dir /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp \\
  --index ${index} \\
  --db_dir /databases/ \\
  --mapout "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}.mapout.bz2" \\
  -o "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}_profile.tsv" \\
  --nproc \${nproc} --offline
sed -i "2d" "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}_profile.tsv"

echo "Running unclassified profiling..."
metaphlan "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}.mapout.bz2" \\
  --input_type mapout \\
  --tmp_dir /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp \\
  --index ${index} \\
  --db_dir /databases/ \\
  -o "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}_unclassified.tsv" \\
  --nproc 1 --offline
sed -i "2d" "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}_profile.tsv"

rm "/output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/\${sample_name}/\${sample_name}.mapout.bz2"
