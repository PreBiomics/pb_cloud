#!/bin/bash
set -euo pipefail
sample_name=${1}
index=mpa_vJan25_CHOCOPhlAnSGB_202503

mkdir -p /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp
mkdir -p /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}

echo "Running statq 0.015..."
metaphlan /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}.mapout.bz2 \
  --input_type mapout --skip_unclassified_estimation \
  --tmp_dir /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp \
  --index ${index} \
  --db_dir /databases/ --stat_q 0.015 \
  -o /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_015.tsv \
  --nproc 1 --offline
sed -i 2d /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_015.tsv

echo "Running statq 0.01..."
metaphlan /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}.mapout.bz2 \
  --input_type mapout --skip_unclassified_estimation \
  --tmp_dir /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp \
  --index ${index} --stat_q 0.01 \
  --db_dir /databases/ \
  -o /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_010.tsv \
  --nproc 1 --offline
sed -i 2d /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_010.tsv

echo "Running statq 0.005..."
metaphlan /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}.mapout.bz2 \
  --input_type mapout --skip_unclassified_estimation \
  --tmp_dir /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/tmp \
  --index ${index} --stat_q 0.005 \
  --db_dir /databases/ \
  -o /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_005.tsv \
  --nproc 1 --offline
sed -i 2d /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}_005.tsv

if ${rm_mapout}; then
  rm /output/tmp/mpa-4.2.4_mpa_vJan25_CHOCOPhlAnSGB_202503/${sample_name}/${sample_name}.mapout.bz2
fi
