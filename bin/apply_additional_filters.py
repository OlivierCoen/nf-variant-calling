#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

import polars as pl
from common import parse_vcf_data, parse_vcf_header

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate and filter VCF files")
    parser.add_argument(
        "--vcf",
        type=Path,
        dest="vcf_file",
        required=True,
        help="Path to VCF file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="Path to output VCF file",
    )
    parser.add_argument(
        "--min-depth-quantile",
        dest="min_depth_quantile",
        type=float,
        default=0.1,
        help="Minimum depth quantile",
    )
    parser.add_argument(
        "--max-depth-quantile",
        dest="max_depth_quantile",
        type=float,
        default=0.9,
        help="Maximum depth quantile",
    )
    return parser.parse_args()


def get_snps_with_too_low_depth(depth_df: pl.DataFrame, quantile: float):
    too_low_df = depth_df.select(pl.all().lt(pl.all().quantile(quantile)))
    return too_low_df.select(
        (pl.sum_horizontal(pl.all() == 0) == 1).alias("not_too_low")
    )


def get_snps_with_too_high_depth(depth_df: pl.DataFrame, quantile: float):
    too_high_df = depth_df.select(pl.all().gt(pl.all().quantile(quantile)))
    return too_high_df.select(
        (pl.sum_horizontal(pl.all() == 0) == 1).alias("not_too_high")
    )


def get_filter_mask(
    vcf_lf: pl.LazyFrame,
    min_depth_quantile: float,
    max_depth_quantile: float,
) -> pl.Series:
    depth_df = vcf_lf.select(
        pl.col("INFO").str.extract(r"DP=(\d+);", 1).cast(pl.Int64).alias("total_depth")
    ).collect()
    nb_snps = len(depth_df)

    not_too_low_depth_filter = get_snps_with_too_low_depth(
        depth_df, quantile=min_depth_quantile
    )
    logger.info(
        f"{nb_snps - not_too_low_depth_filter.sum().item()} SNPs show too low depth"
    )

    not_too_high_depth_filter = get_snps_with_too_high_depth(
        depth_df, quantile=max_depth_quantile
    )

    logger.info(
        f"{nb_snps - not_too_high_depth_filter.sum().item()} SNPs show too high depth"
    )

    return (
        pl.concat(
            [not_too_low_depth_filter, not_too_high_depth_filter],
            how="horizontal",
        )
        .select(mask=pl.all_horizontal(pl.all()))
        .to_series()
    )


def range_starts(lengths: pl.Series) -> pl.Series:
    """
    This function is used to convert relative positions within segments to absolute genome-wide positions.
    The output is a vector where each element represents the starting position of the corresponding segment in a contiguous sequence.
    """
    return pl.concat(
        [pl.Series([1], dtype=lengths.dtype), lengths.head(lengths.len() - 1)]
    ).cum_sum()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")
    vcf_lf = parse_vcf_data(args.vcf_file)
    vcf_header_lines = parse_vcf_header(args.vcf_file)
    nb_original_snps = vcf_lf.select(pl.len()).collect().item()

    if nb_original_snps == 0:
        logger.info("No SNPs found in VCF file")
        sys.exit(0)

    logger.info(f"Parsed {nb_original_snps} SNPs")

    logger.info("Separating read count columns")

    filter_mask = get_filter_mask(
        vcf_lf, args.min_depth_quantile, args.max_depth_quantile
    )
    logger.info(
        f"Kept {filter_mask.sum()} SNPs out of {nb_original_snps} ({filter_mask.sum() / nb_original_snps:.2%})"
    )

    logger.info("Applying filters")
    vcf_lf = vcf_lf.filter(filter_mask)

    logger.info(f"Exporting fitlered VCF to {args.outfile}")

    if vcf_lf.select(pl.len()).collect().item() == 0:
        logger.info("No variants left after filtering")
        sys.exit(0)

    with open(args.outfile, "a") as fout:
        fout.writelines(vcf_header_lines)
        vcf_lf.collect().rename({"CHROM": "#CHROM"}).write_csv(fout, separator="\t")


if __name__ == "__main__":
    main()
