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
sorted_vcf_1=${OUTDIR}/${prefix}.genotype_filtered.1.sorted.vcf.gz
sorted_vcf_2=${OUTDIR}/${prefix}.genotype_filtered.2.sorted.vcf.gz
filtered_vcf=${OUTDIR}/${prefix}.genotype_filtered.vcf
density_scores=${OUTDIR}/${prefix}.density_scores.txt

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

micromamba run -n vcf \
    bcftools sort $filtered_vcf_1 -o $sorted_vcf_1 -O z --write-index

micromamba run -n vcf \
    bcftools sort $filtered_vcf_2 -o  $sorted_vcf_2 -O z --write-index

micromamba run -n vcf \
    bcftools concat $sorted_vcf_1 $sorted_vcf_2 \
        -o $filtered_vcf -a

if [ -f $density_scores ]; then
    rm $density_scores
fi
bin/compute_snp_density.py \
    --vcf $filtered_vcf \
    --out $density_scores \
    --window-size 20000

bin/aggregate_data.py \
    --vcf $filtered_vcf \
    --pvalues $density_scores \
    --design $DESIGN \
    --prefix ${OUTDIR}/$prefix \
    --window-size 20000

cp -r modules/local/dash_app/apps ${OUTDIR}/dash_app
mkdir -p ${OUTDIR}/dash_app/variants/data ${OUTDIR}/dash_app/windows/data

mv ${OUTDIR}/${prefix}.grouped_variants.parquet ${OUTDIR}/dash_app/windows/data
mv ${OUTDIR}/${prefix}.formated_variants.parquet ${OUTDIR}/dash_app/variants/data

