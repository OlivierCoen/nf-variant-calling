#!/usr/bin/env bash

set -euo pipefail

VCF=$1
DESIGN=$2
OUTDIR=$3

mkdir -p $OUTDIR

prefix=$(basename $VCF .vcf.gz).snp_indel
RO_counts=${OUTDIR}/${prefix}.RO_counts.parquet
AO_counts=${OUTDIR}/${prefix}.AO_counts.parquet
pvalues=${OUTDIR}/${prefix}.pvalues.txt

bin/separate_vcf_data.py \
    --vcf $VCF

mv RO_counts.parquet $RO_counts
mv AO_counts.parquet $AO_counts

if [ -f $pvalues ]; then
    rm $pvalues
fi
micromamba run -n stat_test bin/compute_statistical_test.R \
    --method fet \
    --RO $RO_counts \
    --AO $AO_counts \
    --design $DESIGN \
    --out $pvalues

bin/aggregate_data.py \
    --vcf $filtered_vcf \
    --pvalues $pvalues \
    --design $DESIGN \
    --prefix ${OUTDIR}/$prefix \
    --window-size 20000

cp -r modules/local/dash_app ${OUTDIR}/dash_app
mkdir -p ${OUTDIR}/dash_app/apps/variants/data ${OUTDIR}/dash_app/apps/windows/data

mv ${OUTDIR}/${prefix}.grouped_variants.parquet ${OUTDIR}/dash_app/apps/windows/data
mv ${OUTDIR}/${prefix}.formated_variants.parquet ${OUTDIR}/dash_app/apps/variants/data

