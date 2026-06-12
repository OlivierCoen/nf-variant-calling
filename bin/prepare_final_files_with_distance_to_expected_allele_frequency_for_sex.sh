#!/usr/bin/env bash

set -euo pipefail

VCF=$1
DESIGN=$2
GENOTYPES_1=$3
GENOTYPES_2=$4
OUTDIR=$5
STRICT_MODE=$6

mkdir -p $OUTDIR

prefix=$(basename $VCF .vcf.gz).snp_indel
filtered_vcf_1=${OUTDIR}/${prefix}.genotype_filtered.1.vcf
filtered_vcf_2=${OUTDIR}/${prefix}.genotype_filtered.2.vcf
scores_1=${OUTDIR}/${prefix}.distance_scores.1.txt
scores_2=${OUTDIR}/${prefix}.distance_scores.2.txt

if [ "$STRICT_MODE" = "true" ]; then
    strict_flag="--strict"
else
    strict_flag=""
fi
bin/filter_by_genotype.py \
    --vcf $VCF \
    --genotypes $GENOTYPES_1 \
    --out $filtered_vcf_1 \
    $strict_flag

bin/filter_by_genotype.py \
    --vcf $VCF \
    --genotypes $GENOTYPES_2 \
    --out $filtered_vcf_2 \
    $strict_flag

bin/compute_distance_to_expected_allele_frequency.py \
    --vcf $filtered_vcf_1 \
    --genotypes $GENOTYPES_1 \
    --out $scores_1

bin/compute_distance_to_expected_allele_frequency.py \
    --vcf $filtered_vcf_2 \
    --genotypes $GENOTYPES_2 \
    --out $scores_2

bin/aggregate_data2.py \
    --vcf1 $filtered_vcf_1 \
    --vcf2 $filtered_vcf_2 \
    --pvalues1 $scores_1 \
    --pvalues2 $scores_2 \
    --design $DESIGN \
    --prefix ${OUTDIR}/$prefix \
    --window-size 20000

cp -r modules/local/dash_app/apps ${OUTDIR}/dash_app
mkdir -p ${OUTDIR}/dash_app/variants/data ${OUTDIR}/dash_app/windows/data

mv ${OUTDIR}/${prefix}.grouped_variants.parquet ${OUTDIR}/dash_app/windows/data
mv ${OUTDIR}/${prefix}.formated_variants.parquet ${OUTDIR}/dash_app/variants/data

