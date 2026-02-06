#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from common import parse_vcf_data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGZISE = (12, 6)
WINDOW_SIZE = int(1e6)


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
        "--filtered-vcf",
        type=Path,
        dest="filtered_vcf_file",
        required=True,
        help="Path to filtered VCF file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="Path to output VCF file",
    )
    parser.add_argument(
        "--fai",
        type=Path,
        dest="fai_file",
        required=True,
        help="Fai index file of the reference genome",
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


def get_genome_positions(
    vcf_file: Path, fai_file: Path, window_size: int
) -> pl.DataFrame:

    chrom_df = (
        pl.scan_csv(fai_file, has_header=False, separator="\t")
        .rename({"column_1": "chrom", "column_2": "len"})
        .with_columns(pl.col("chrom").str.strip_chars())
        .select(["chrom", "len"])
        .collect()
    )

    scaff_starts = range_starts(chrom_df["len"])
    chroms_with_starts = chrom_df.with_columns(scaffold_start=scaff_starts)

    vcf_lf = parse_vcf_data(vcf_file)

    return (
        vcf_lf.select(["CHROM", "POS"])
        .collect()
        .join(chroms_with_starts, left_on="CHROM", right_on="chrom", how="left")
        .with_columns(
            (pl.col("POS") + pl.col("scaffold_start") - 1).alias("genome_position")
        )
        .with_columns(
            ((pl.col("genome_position") // window_size) * window_size)
            .cast(pl.Int64)
            .alias("window")
        )
    )


def get_nb_snps_per_window(genome_positions: pl.DataFrame):
    return genome_positions.group_by("window").agg(pl.len().alias("N")).sort("window")


def plot_effect_of_filters_on_snps(
    genome_positions: pl.DataFrame,
    filtered_genome_positions: pl.DataFrame,
    outfile: Path,
):
    n_snps = get_nb_snps_per_window(genome_positions)
    n_snps_filtered = get_nb_snps_per_window(filtered_genome_positions)

    n_snps = n_snps.join(
        n_snps_filtered, on="window", how="left", suffix="_filtered"
    ).fill_null(0)

    fig, ax = plt.subplots(figsize=FIGZISE)
    plt.scatter(n_snps["N"], n_snps["N_filtered"])
    plt.xlabel("SNP count (before filtering)")
    plt.ylabel("SNP count (after filtering)")
    plt.title("Effect of filters on SNP density per 1 Mb window")

    plt.savefig(outfile, bbox_inches="tight")
    # closes current figure window to save memory
    plt.close()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")
    genome_positions = get_genome_positions(args.vcf_file, args.fai_file, WINDOW_SIZE)
    filtered_genome_positions = get_genome_positions(
        args.filtered_vcf_file, args.fai_file, WINDOW_SIZE
    )

    logger.info("Plotting effect of filters on SNP density")
    plot_effect_of_filters_on_snps(
        genome_positions, filtered_genome_positions, args.outfile
    )

    logger.info("Done")


if __name__ == "__main__":
    main()
