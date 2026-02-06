#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
from common import parse_vcf_data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

RO_OUTFILE = "RO_counts.csv"
AO_OUTFILE = "AO_counts.csv"


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
    return parser.parse_args()


def get_position_in_format(vcf_lf: pl.LazyFrame, info: str) -> int:
    fmt_series = vcf_lf.select("FORMAT").collect().to_series()
    if fmt_series.unique().len() > 1:
        raise ValueError(f"More than one format found: {fmt_series.unique()}")
    fmt = fmt_series.unique().item()
    return fmt.split(":").index(info)


def extract_counts(
    vcf_lf: pl.LazyFrame, sample_cols: list[str], idx: int
) -> pl.LazyFrame:
    return vcf_lf.select(
        pl.col(sample_cols).str.split(":").list[idx].replace(".", None).cast(pl.Int64)
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")

    vcf_lf = parse_vcf_data(args.vcf_file)

    RO_index = get_position_in_format(vcf_lf, "RO")
    AO_index = get_position_in_format(vcf_lf, "AO")

    sample_cols = vcf_lf.collect_schema().names()[9:]
    RO_lf = extract_counts(vcf_lf, sample_cols, RO_index)
    AO_lf = extract_counts(vcf_lf, sample_cols, AO_index)

    RO_lf.collect().write_csv(RO_OUTFILE)
    AO_lf.collect().write_csv(AO_OUTFILE)


if __name__ == "__main__":
    main()
