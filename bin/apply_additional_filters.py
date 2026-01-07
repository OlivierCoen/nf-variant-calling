#!/usr/bin/env python3


import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

FIGZISE = (12, 6)

OUTFILE_SUFFIX = "filtered.vcf"


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
        "--fai",
        type=Path,
        dest="fai_file",
        required=True,
        help="Fai index file of the reference genome",
    )
    parser.add_argument(
        "--max-depth-quantile",
        dest="max_depth_quantile",
        type=float,
        default=0.9,
        help="Maximum depth quantile",
    )
    parser.add_argument(
        "--min-freq",
        dest="min_allele_frequency",
        type=float,
        default=0.1,
        help="Minimum allele frequency",
    )
    parser.add_argument(
        "--max-freq",
        dest="max_allele_frequency",
        type=float,
        default=0.99,
        help="Maximum allele frequency",
    )
    parser.add_argument(
        "--min-qual",
        dest="min_quality",
        type=float,
        default=10,
        help="Minimum quality score",
    )
    parser.add_argument(
        "--window-size",
        dest="window_size",
        type=int,
        default=1e6,
        help="Window size",
    )
    return parser.parse_args()


def clean_duplicated_sample_name(sample_name: str) -> str:
    """
    Clean string like LU26_M_LU26_M to get only LU26_M.
    These duplicated strings are caused by nf-core/sarek
    """
    lst = sample_name.split("_")
    if lst[:2] == lst[2:]:
        return "_".join(lst[:2])
    return sample_name


def parse_vcf(vcf_file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        vcf_file, separator="\t", has_header=True, null_values="."
    ).rename(lambda col: col.split("]")[1])


def get_read_count_columns(vcf_lf: pl.LazyFrame, suffix: str) -> pl.DataFrame:
    return (
        vcf_lf.select(pl.col(f"^.*:{suffix}$"))
        .rename(lambda col: col.split(":")[0])
        .rename(clean_duplicated_sample_name)
        # .fill_null(0)
        .collect()
    )


def mask_snps_with_too_low_depth(depth_df: pl.DataFrame, threshold: int):
    return depth_df.select(not_too_low=pl.any_horizontal(pl.all() >= threshold))


def mask_snps_with_too_high_depth(depth_df: pl.DataFrame, quantile: float):
    # Check rows with no values exceeding the 90th percentile
    too_high = depth_df.select(pl.all().gt(pl.all().quantile(quantile)))
    return too_high.select(not_too_high=(pl.sum_horizontal(pl.all() == 0) == 1))


def mask_too_rare_snps(
    AO_df: pl.DataFrame,
    depth_df: pl.DataFrame,
    min_frequency: float,
    max_frequency: float,
):
    freq_df = AO_df / depth_df
    return freq_df.select(
        not_too_rare=(
            pl.sum_horizontal(pl.all().is_between(min_frequency, max_frequency)) == 1
        ).fill_null(False)
    )


def mask_snps_with_low_quality(vcf_df: pl.DataFrame, min_quality: float):
    return vcf_df.select(not_low_quality=(pl.col("QUAL") >= min_quality))


def get_filters(
    vcf_df: pl.DataFrame,
    AO_df: pl.DataFrame,
    depth_df: pl.DataFrame,
    args,
):
    too_low_depth_filter = mask_snps_with_too_low_depth(
        depth_df, threshold=args.min_depth_threshold
    )
    logger.info(
        f"{len(too_low_depth_filter) - too_low_depth_filter.sum().item()} SNPs show too low depth"
    )

    too_high_depth_filter = mask_snps_with_too_high_depth(
        depth_df, quantile=args.max_depth_quantile
    )
    logger.info(
        f"{len(too_low_depth_filter) - too_high_depth_filter.sum().item()} SNPs show too high depth"
    )

    too_rare_filter = mask_too_rare_snps(
        AO_df,
        depth_df,
        min_frequency=args.min_allele_frequency,
        max_frequency=args.max_allele_frequency,
    )
    logger.info(
        f"{len(too_low_depth_filter) - too_rare_filter.sum().item()} SNPs are too rare"
    )

    low_quality_filter = mask_snps_with_low_quality(
        vcf_df,
        min_quality=args.min_quality,
    )
    logger.info(
        f"{len(too_low_depth_filter) - low_quality_filter.sum().item()} SNPs have low quality"
    )

    return (
        pl.concat(
            [
                too_low_depth_filter,
                too_high_depth_filter,
                too_rare_filter,
                low_quality_filter,
            ],
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


def get_genome_positions(vcf_df: pl.DataFrame, args):
    chrom_df = (
        pl.scan_csv(args.fai_file, has_header=False, separator="\t")
        .rename({"column_1": "chrom", "column_2": "len"})
        .with_columns(pl.col("chrom").str.strip_chars())
        .select(["chrom", "len"])
        .collect()
    )

    scaff_starts = range_starts(chrom_df["len"])
    chroms_with_starts = chrom_df.with_columns(scaffold_start=scaff_starts)

    return (
        vcf_df.join(chroms_with_starts, left_on="CHROM", right_on="chrom", how="left")
        .with_columns(
            (pl.col("POS") + pl.col("scaffold_start") - 1).alias("genome_position")
        )
        .with_columns(
            ((pl.col("genome_position") // args.window_size) * args.window_size)
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


def plot_effect_of_filters_on_snps(vcf_df, filters, args):
    genome_positions = get_genome_positions(vcf_df, args)
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
    vcf_lf = parse_vcf(args.vcf_file)

    logger.info("Separating read count columns")
    AO_df = get_read_count_columns(vcf_lf, suffix="AO")
    RO_df = get_read_count_columns(vcf_lf, suffix="RO")
    # vcf_df = vcf_lf.select(pl.exclude(["REF", "ALT", "^.*:AO$", "^.*:RO$"])).collect()
    vcf_df = vcf_lf.select(pl.exclude(["^.*:AO$", "^.*:RO$"])).collect()

    depth_df = AO_df + RO_df

    filters = get_filters(vcf_df, AO_df, depth_df, args)
    logger.info(f"Kept {filters.sum()} SNPs out of {len(vcf_df)}")

    logger.info("Plotting effect of filters on SNP density")
    plot_effect_of_filters_on_snps(vcf_df, filters, args)

    logger.info("Applying filters")
    vcf_df = vcf_df.filter(filters)

    outfile = Path(args.vcf_file).with_suffix(OUTFILE_SUFFIX)
    logger.info(f"Exporting VCF to {outfile}")
    vcf_df.write_csv(outfile)


if __name__ == "__main__":
    main()
