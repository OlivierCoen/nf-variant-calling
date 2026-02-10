#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
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
        "--variants",
        type=Path,
        dest="variant_file",
        required=True,
        help="Path to variant file",
    )
    parser.add_argument(
        "--pvalues",
        type=Path,
        dest="pvalue_file",
        required=True,
        help="Path to file containing p-values",
    )
    parser.add_argument(
        "--RO",
        type=Path,
        dest="RO_file",
        required=True,
        help="Path to file containing reference counts",
    )
    parser.add_argument(
        "--AO",
        type=Path,
        dest="AO_file",
        required=True,
        help="Path to file containing reference counts",
    )
    parser.add_argument(
        "--design",
        type=Path,
        dest="design_file",
        required=True,
        help="Path to design file",
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


def get_allele_count_expr(populations: list[str]):
    for pop in populations:
        yield (
            f"{pop}_F: "
            + pl.col(f"{pop}_F_RO_mean").round(2).cast(pl.String)
            + " (±"
            + pl.col(f"{pop}_F_RO_std").round(2).cast(pl.String)
            + ") / "
            + pl.col(f"{pop}_F_AO_mean").round(2).cast(pl.String)
            + " (±"
            + pl.col(f"{pop}_F_AO_std").round(2).cast(pl.String)
            + f") <br>{pop}_M: "
            + pl.col(f"{pop}_M_RO_mean").round(2).cast(pl.String)
            + " (±"
            + pl.col(f"{pop}_M_RO_std").round(2).cast(pl.String)
            + ") / "
            + pl.col(f"{pop}_M_AO_mean").round(2).cast(pl.String)
            + " (±"
            + pl.col(f"{pop}_M_AO_std").round(2).cast(pl.String)
            + ")"
        ).alias(pop)


def add_total_depth(lf: pl.LazyFrame) -> pl.Series:
    return lf.with_columns(
        pl.col("INFO").str.extract(r"DP=(\d+)").cast(pl.Int64).alias("total_depth")
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

    try:
        window_size = int(float(args.window_size))
        logger.info(f"Window size: {window_size}")
    except ValueError:
        raise TypeError(f"Could not cast {args.window_size} to integer.")

    logger.info("Parsing VCF file")
    variant_lf = pl.scan_parquet(args.variant_file)

    logger.info("Parsing p-values")
    pvalues_lf = parse_pvalues(args.pvalue_file)

    logger.info("Parsing allele counts")
    RO_lf = pl.scan_parquet(args.RO_file).select(pl.all().name.suffix("_RO"))
    AO_lf = pl.scan_parquet(args.AO_file).select(pl.all().name.suffix("_AO"))

    RO_cols = RO_lf.collect_schema().names()
    AO_cols = AO_lf.collect_schema().names()

    design_df = pl.read_csv(args.design_file)
    populations = design_df["population"].unique().to_list()

    logger.info("Associating SNPs to windows")
    variant_lf = add_windows(variant_lf, window_size)

    nb_variants = get_length(variant_lf)
    nb_pvalues = get_length(pvalues_lf)
    if nb_variants != nb_pvalues:
        raise ValueError(
            f"Number of variants ({nb_variants}) and number of pvalues ({nb_pvalues}) do not match."
        )

    lf = pl.concat([variant_lf, pvalues_lf, RO_lf, AO_lf], how="horizontal")

    logger.info("Getting list of contigs")
    contigs = get_chromosomes(lf)

    logger.info("Computing total depth")
    lf = add_total_depth(lf)

    logger.info(
        f"Computing quantile {QUANTILE} of pvalue for each pair of contig & window"
    )
    chrom_files = []
    for i, contig in tqdm(enumerate(contigs), total=len(contigs)):
        contig_lf = (
            lf.filter(pl.col("CHROM") == contig)
            .group_by(  # inserting CHROM in the group_by in order to get it back after the agg
                ["CHROM", "window"]
            )
            .agg(
                pl.col("pvalue").quantile(QUANTILE),
                pl.col("QUAL").mean().alias("quality"),
                pl.col("total_depth").mean(),
                pl.col(RO_cols + AO_cols).mean().name.suffix("_mean"),
                pl.col(RO_cols + AO_cols).std().name.suffix("_std"),
            )
            .with_columns(get_allele_count_expr(populations))
            .with_columns(
                pl.concat_str(populations, separator="<br>").alias("allele_counts")
            )
            .rename({"CHROM": "chromosome", "window": "position"})
            .select(
                [
                    "chromosome",
                    "position",
                    "pvalue",
                    "quality",
                    "total_depth",
                    "allele_counts",
                ]
            )
        )

        chrom_file = f"chrom_{i}.parquet"
        contig_lf.sink_parquet(chrom_file)
        chrom_files.append(chrom_file)

    chrom_lfs = [pl.scan_parquet(file) for file in chrom_files]
    total_lf = pl.concat(chrom_lfs, how="vertical")

    logger.info(f"Saving data to {args.outfile}")
    total_lf.sink_parquet(args.outfile)

    logger.info("Removing temporary files")
    for file in chrom_files:
        Path(file).unlink()

    # quantile_lf.show_graph(engine="streaming", plan_stage="physical")


if __name__ == "__main__":
    main()
