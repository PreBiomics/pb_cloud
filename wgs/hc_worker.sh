#!/bin/bash

REFERENCE=$1
BAM=$2
LIST=$3
OUT=$4

ARGS=()

while read L
do
	    ARGS+=("-L")
	        ARGS+=("$L")
	done < "$LIST"

	gatk HaplotypeCaller \
		    -R "${REFERENCE}.fa" \
		        -I "$BAM" \
			    "${ARGS[@]}" \
			        -ERC GVCF \
				    --native-pair-hmm-threads 1 \
				        -O "$OUT"
