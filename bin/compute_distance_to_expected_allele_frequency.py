#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
from tqdm import tqdm
from common import parse_vcf_data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GENOTYPES_TO_ALLELE_FREQUENCIES = {
    "0/0": 0,
    "0/1": 0.5,
    "1/0": 0.5,
    "1/1": 1,
}

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Get reference and alternative counts")
    parser.add_argument(
        "--vcf",
        type=Path,
        dest="vcf_file",
        required=True,
        help="Path to VCF file",
    )
    parser.add_argument(
        "--genotypes",
        type=Path,
        dest="genotype_file",
        required=True,
        help="Path to genotypes file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="Path to output file",
    )
    return parser.parse_args()


def get_position_in_format(vcf_lf: pl.LazyFrame, info: str) -> int:
    fmt_series = vcf_lf.select("FORMAT").collect().to_series()
    if fmt_series.unique().len() > 1:
        raise ValueError(f"More than one format found: {fmt_series.unique()}")
    fmt = fmt_series.unique().item()
    return fmt.split(":").index(info)


def extract_first_alternative_allele_frequency(
    vcf_lf: pl.LazyFrame, sample_cols: list[str], idx: int
) -> pl.DataFrame:
    return vcf_lf.select(
        pl.col(sample_cols).str.split(":").list[idx].replace(".", None).cast(pl.Float64)
    ).collect()


def compute_distance():
    pass

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")

    vcf_lf = parse_vcf_data(args.vcf_file)
    
    VAF1_index = get_position_in_format(vcf_lf, "VAF1")
 
    sample_cols = vcf_lf.collect_schema().names()[9:]
    VAF1_df = extract_first_alternative_allele_frequency(vcf_lf, sample_cols, VAF1_index)
    
    sample_genotypes_df = pl.read_csv(args.genotype_file)
    sample_genotypes_df = sample_genotypes_df.with_columns(
        genotypes=pl.col("genotypes").replace("het", "0/1:1/0").replace("hom", "0/0:1/1").str.split(":")
    )
    sample_to_genotypes = { d["sample"]: d["genotypes"] for d in sample_genotypes_df.to_dicts() }
    expected_frequencies = { 
        sample: list(set([GENOTYPES_TO_ALLELE_FREQUENCIES[g] for g in genotypes])) 
        for sample, genotypes in sample_to_genotypes.items() 
    }
    
    logger.info(f"Computing distance to expected allele frequency for {len(sample_cols)} samples")
    for sample in tqdm(sample_cols):
        freqs = expected_frequencies[sample]
        VAF1_df = VAF1_df.with_columns(
            pl.col(sample).map_elements(
                    lambda x: min((x - freq) ** 2 for freq in freqs)
                ).alias(f"{sample}_dist")
        )

    distance_cols = [f"{sample}_dist" for sample in sample_cols]
    dist_df = VAF1_df.select(distance_cols)
    del VAF1_df

    logger.info("Computing number of non-null allele frequencies per SNP")
    computations = dist_df.with_columns(
        pl.fold(
            acc=pl.lit(0),
            function=lambda acc, x: acc + x.is_not_null().cast(pl.Int32),
            exprs=pl.all(),
        ).alias("not_null_count")
    )

    logger.info("Computing distance scores")
    # mean and not sum in order to compensate for missing values
    distance_scores = dist_df.sum_horizontal(ignore_nulls=True).alias("sum_of_squares")
    
    logger.info("Computing final scores")
    
    K = 6
    final_scores = (
        computations.hstack([distance_scores])
        .select(
            (pl.col("sum_of_squares").sqrt() / pl.col("not_null_count").cast(pl.Float32)) ** K
        )
        .to_series()
        .replace(0, 1E-40)
        .to_list()
    )
    
    with open(args.outfile, 'w') as fout:
        fout.writelines([f'{score:.10f}\n' for score in final_scores])


if __name__ == "__main__":
    main()
