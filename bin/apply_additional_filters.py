#!/usr/bin/env python3

import argparse
import gzip
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGZISE = (12, 6)


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
        "--fai",
        type=Path,
        dest="fai_file",
        required=True,
        help="Fai index file of the reference genome",
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
    parser.add_argument(
        "--window-size",
        dest="window_size",
        type=int,
        default=1e6,
        help="Window size",
    )
    return parser.parse_args()


def parse_vcf_header(vcf_file: Path) -> list:
    header_lines = []
    with gzip.open(vcf_file, "rb") as fin:
        for line in fin.readlines():
            line = line.decode("utf-8")
            if line.startswith("##"):
                header_lines.append(line)
            else:
                break
    return header_lines


def parse_vcf_data(vcf_file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        vcf_file, separator="\t", has_header=True, comment_prefix="##"
    ).rename({"#CHROM": "CHROM"})


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
    vcf_lf: pl.LazyFrame, fai_file: Path, window_size: int
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


def get_nb_filtered_snps_per_window(genome_positions: pl.DataFrame, mask: pl.Series):
    """Count SNPs per window (after filtering)"""
    return (
        genome_positions.filter(mask)
        .group_by("window")
        .agg(pl.len().alias("N"))
        .sort("window")
    )


def plot_effect_of_filters_on_snps(filters: pl.Series, genome_positions: pl.DataFrame):
    n_snps = get_nb_snps_per_window(genome_positions)
    n_snps_filtered = get_nb_filtered_snps_per_window(genome_positions, filters)

    n_snps = n_snps.join(
        n_snps_filtered, on="window", how="left", suffix="_filtered"
    ).fill_null(0)

    fig, ax = plt.subplots(figsize=FIGZISE)
    plt.scatter(n_snps["N"], n_snps["N_filtered"])
    plt.xlabel("SNP count (before filtering)")
    plt.ylabel("SNP count (after filtering)")
    plt.title("Effect of filters on SNP density per 1 Mb window")

    outfile = "effect_of_filters_on_snps.png"
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
    vcf_lf = parse_vcf_data(args.vcf_file)
    vcf_header_lines = parse_vcf_header(args.vcf_file)
    nb_original_snps = vcf_lf.select(pl.len()).collect().item()
    logger.info(f"Parsed {nb_original_snps} SNPs")

    logger.info("Separating read count columns")

    filter_mask = get_filter_mask(
        vcf_lf, args.min_depth_quantile, args.max_depth_quantile
    )
    logger.info(
        f"Kept {filter_mask.sum()} SNPs out of {nb_original_snps} ({filter_mask.sum() / nb_original_snps:.2%})"
    )

    logger.info("Plotting effect of filters on SNP density")
    genome_positions = get_genome_positions(vcf_lf, args.fai_file, args.window_size)
    plot_effect_of_filters_on_snps(filter_mask, genome_positions)

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
