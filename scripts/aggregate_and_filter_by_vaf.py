#!/usr/bin/env python3

import argparse
import logging
from multiprocessing import Value
from pathlib import Path

import polars as pl

pl.Config.set_streaming_chunk_size(1e6)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GROUPED_VARIANTS_OUTFILE_SUFFIX = "grouped_variants.parquet"
VARIANTS_OUTFILE_SUFFIX = "formated_variants.parquet"

QUANTILE = 0.05

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Get reference and alternative counts")
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
        "--prefix",
        required=True,
        help="Prefix for output files",
    )
    return parser.parse_args()


def parse_pvalues(cmh_pvalues_file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        cmh_pvalues_file, has_header=False, new_columns=["pvalue"], null_values=["NA"]
    )


def get_allele_count_expr(populations: list[str]):
    for pop in populations:
        yield (
            f"{pop}_F: "
            + pl.col(f"{pop}_F_RO").cast(pl.String)
            + " / "
            + pl.col(f"{pop}_F_AO").cast(pl.String)
            + f"<br>{pop}_M: "
            + pl.col(f"{pop}_M_RO").cast(pl.String)
            + " / "
            + pl.col(f"{pop}_M_AO").cast(pl.String)
        ).alias(pop)


def add_total_depth(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.with_columns(
        pl.col("INFO").str.extract(r"DP=(\d+)").cast(pl.Int64).alias("total_depth")
    )


def get_length(lf: pl.LazyFrame) -> int:
    return lf.select(pl.len()).collect().item(0, 0)


def get_position_in_format(vcf_lf: pl.LazyFrame, info: str) -> int:
    fmt_series = vcf_lf["FORMAT"]
    if fmt_series.unique().len() > 1:
        raise ValueError(f"More than one format found: {fmt_series.unique()}")
    fmt = fmt_series.unique().item()
    return fmt.split(":").index(info)


def extract_vaf(vcf_lf: pl.LazyFrame, sample_cols: list[str], idx: int) -> pl.LazyFrame:
    return vcf_lf.select(
        pl.col("POS"),
        pl.col(sample_cols)
        .str.split(":")
        .list[idx]
        .replace(".", None)
        .cast(pl.Float64),
        # .cast(pl.Float64),
    )


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")
    variant_lf = pl.scan_parquet(args.variant_file)

    logger.info("Parsing p-values")
    pvalues_lf = parse_pvalues(args.pvalue_file)

    logger.info("Parsing allele counts")
    design_df = pl.read_csv(args.design_file)
    populations = design_df["population"].unique().to_list()

    logger.info("Parsing allele counts")
    RO_lf = pl.scan_parquet(args.RO_file).select(pl.all().name.suffix("_RO"))
    AO_lf = pl.scan_parquet(args.AO_file).select(pl.all().name.suffix("_AO"))

    nb_variants = get_length(variant_lf)
    nb_pvalues = get_length(pvalues_lf)
    if nb_variants != nb_pvalues:
        raise ValueError(
            f"Number of variants ({nb_variants}) and number of pvalues ({nb_pvalues}) do not match."
        )

    df = pl.concat(
        [
            variant_lf,
            pvalues_lf,
            RO_lf,
            AO_lf,
        ],
        how="horizontal",
    ).collect()

    logger.info("Getting list of contigs")

    logger.info("Computing total depth")
    df = add_total_depth(df)

    VAF_index = get_position_in_format(df, "VAF")

    sample_cols = df.collect_schema().names()[9:15]
    vaf_lf = extract_vaf(df, sample_cols, VAF_index)

    valid_positions = (
        vaf_lf.filter(
            [
                (pl.col(col) <= 0.1) | (pl.col(col) >= 0.9)
                for col in sample_cols
                if col.endswith("_F")
            ]
        )
        .filter(
            [
                (pl.col(col) >= 0.4) & (pl.col(col) <= 0.6)
                for col in sample_cols
                if col.endswith("_M")
            ]
        )
        .select("POS")
        .to_series()
        .to_list()
    )

    filtered_lf = df.filter(pl.col("POS").is_in(valid_positions))

    variant_lf = (
        filtered_lf.with_columns(get_allele_count_expr(populations))
        .with_columns(
            pl.concat_str(populations, separator="<br>").alias("allele_counts")
        )
        .rename({"CHROM": "chromosome", "POS": "position", "QUAL": "quality"})
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
    print(len(variant_lf))
    logger.info("Saving data")
    variant_lf.write_parquet(f"{args.prefix}.{VARIANTS_OUTFILE_SUFFIX}")


if __name__ == "__main__":
    main()
