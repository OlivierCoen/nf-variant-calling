#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
from common import parse_vcf_data
from tqdm import tqdm

pl.Config.set_streaming_chunk_size(1e6)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

QUANTILE = 0.05


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
        "--pvalues",
        type=Path,
        dest="pvalue_file",
        required=True,
        help="Path to file containing p-values",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="Path to output parquet file",
    )
    parser.add_argument(
        "--window-size",
        dest="window_size",
        type=int,
        required=True,
        help="Window size",
    )
    return parser.parse_args()


def parse_pvalues(cmh_pvalues_file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        cmh_pvalues_file, has_header=False, new_columns=["pvalue"], null_values=["NA"]
    )


def add_windows(lf: pl.LazyFrame, window_size: int) -> pl.LazyFrame:
    return lf.with_columns(
        ((pl.col("POS") // window_size) * window_size + int(window_size / 2)).alias(
            "window"
        )
    )


def get_length(lf: pl.LazyFrame) -> int:
    return lf.select(pl.len()).collect().item(0, 0)


def get_chromosomes(lf: pl.LazyFrame) -> list[str]:
    return lf.select("CHROM").unique().collect().to_series().to_list()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")
    vcf_lf = parse_vcf_data(args.vcf_file)
    vcf_lf.sink_parquet("vcf.parquet")
    logger.info("Parsing p-values")
    pvalues_lf = parse_pvalues(args.pvalue_file)

    logger.info("Associating SNPs to windows")
    vcf_lf = add_windows(vcf_lf, args.window_size)

    nb_variants = get_length(vcf_lf)
    nb_pvalues = get_length(pvalues_lf)
    if nb_variants != nb_pvalues:
        raise ValueError(
            f"Number of variants ({nb_variants}) and number of pvalues ({nb_pvalues}) do not match."
        )

    vcf_lf = pl.concat([vcf_lf, pvalues_lf], how="horizontal")

    logger.info("Getting list of contigs")
    contigs = get_chromosomes(vcf_lf)

    logger.info(
        f"Computing quantile {QUANTILE} of pvalue for each pair of contig & window"
    )
    quantile_files = []
    for i, contig in tqdm(enumerate(contigs), total=len(contigs)):
        file = f"quantile_{i}.parquet"
        lf = (
            vcf_lf.filter(pl.col("CHROM") == contig)
            .group_by("window")
            .agg(pl.col("pvalue").quantile(QUANTILE).alias("pvalue_quantile"))
        )
        lf.sink_parquet(file)
        quantile_files.append(file)

    quantile_lfs = [pl.scan_parquet(file) for file in quantile_files]
    quantile_lf = pl.concat(quantile_lfs, how="vertical")

    logger.info(f"Saving data to {args.outfile}")
    quantile_lf.sink_parquet(args.outfile)

    logger.info("Removing temporary files")
    for file in quantile_files:
        Path(file).unlink()

    # quantile_lf.show_graph(engine="streaming", plan_stage="physical")


if __name__ == "__main__":
    main()
