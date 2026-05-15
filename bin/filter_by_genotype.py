#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl

from common import parse_vcf_data, get_position_in_format, parse_vcf_header

pl.Config.set_streaming_chunk_size(int(1e6))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

OUTFILE_SUFFIX = ".genotype_filtered.vcf"


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
        "--genotypes",
        dest="genotype_file",
        type=Path,
        required=True,
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()



def add_windows(lf: pl.LazyFrame, window_size: int) -> pl.LazyFrame:
    return lf.with_columns(
        ((pl.col("POS") // window_size) * window_size + int(window_size / 2)).alias(
            "window"
        )
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

    header = parse_vcf_header(args.vcf_file)

    sample_genotypes_df = pl.read_csv(args.genotype_file)
    sample_genotypes_df = sample_genotypes_df.with_columns(
        genotypes=pl.col("genotypes").replace("het", "0/1:1/0").replace("hom", "0/0:1/1").str.split(":")
    )
    sample_to_genotypes = { d["sample"]: d["genotypes"] for d in sample_genotypes_df.to_dicts() }

    if not args.strict:
        sample_to_genotypes = {
            k: lst + ["."] for k, lst in sample_to_genotypes.items()
        }

    GT_index = get_position_in_format(vcf_lf, "GT")
    samples = vcf_lf.collect_schema().names()[9:]

    logger.info("Filtering genotypes")
    for sample in samples:
        expected_genotypes = sample_to_genotypes[sample]
        vcf_lf = vcf_lf.filter(
            pl.col(sample).str.split(":").list[GT_index].is_in(expected_genotypes)
        )

    outfile = args.vcf_file.name.rstrip(".gz").replace(".vcf", OUTFILE_SUFFIX)
    logger.info(f"Writing filtered data to {outfile}")
    with open(outfile, 'w') as fout:
        fout.writelines(header)
        vcf_lf.rename({"CHROM": "#CHROM"}).sink_csv(fout, separator="\t")


if __name__ == "__main__":
    main()
