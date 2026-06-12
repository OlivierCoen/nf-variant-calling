#!/usr/bin/env bash

set -euo pipefail

VCF=$1
DESIGN=$2
GENOTYPES=$3
OUTDIR=$4
STRICT_MODE=$5

mkdir -p $OUTDIR

prefix=$(basename $VCF .vcf.gz).snp_indel
filtered_vcf=${OUTDIR}/${prefix}.genotype_filtered.vcf
scores=${OUTDIR}/${prefix}.distance_scores.txt

if [ "$STRICT_MODE" = "true" ]; then
    strict_flag="--strict"
else
    strict_flag=""
fi
bin/filter_by_genotype.py \
    --vcf $VCF \
    --genotypes $GENOTYPES \
    --out $filtered_vcf \
    $strict_flag

if [ -f $scores ]; then
    rm $scores
fi
bin/compute_distance_to_expected_allele_frequency.py \
    --vcf $filtered_vcf \
    --genotypes $GENOTYPES \
    --out $scores

bin/aggregate_data.py \
    --vcf $filtered_vcf \
    --pvalues $scores \
    --design $DESIGN \
    --prefix ${OUTDIR}/$prefix \
    --window-size 20000

cp -r modules/local/dash_app/apps ${OUTDIR}/dash_app
mkdir -p ${OUTDIR}/dash_app/variants/data ${OUTDIR}/dash_app/windows/data

mv ${OUTDIR}/${prefix}.grouped_variants.parquet ${OUTDIR}/dash_app/windows/data
mv ${OUTDIR}/${prefix}.formated_variants.parquet ${OUTDIR}/dash_app/variants/data

