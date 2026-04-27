#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
from tqdm import tqdm

from common import parse_vcf_data, get_position_in_format

pl.Config.set_streaming_chunk_size(int(1e6))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

GROUPED_VARIANTS_OUTFILE_SUFFIX = "grouped_variants.parquet"
VARIANTS_OUTFILE_SUFFIX = "formated_variants.parquet"

QUANTILE = 0.05

logger.info(f"Polars pool size: {pl.thread_pool_size()}")


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
        "--design",
        type=Path,
        dest="design_file",
        required=True,
        help="Path to design",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for output files",
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


def format_count_per_sample(phenotype_to_samples: dict):
    expr = ( pl.lit("Counts (REF / ALT / ...):<br>") )
    for phenotype, samples in phenotype_to_samples.items():
        expr += (pl.lit(f"  Phenotype: {phenotype}:<br>"))
        for sample in samples:
            ad_col = f"{sample}_AD"
            expr += (
                pl.lit(f"    {sample}: ")
                + pl.when(
                    pl.col(ad_col).is_not_null()
                ).then(
                    pl.col(ad_col).list.join(" / ")
                ).otherwise(
                    pl.lit("")
                )
                + pl.lit("<br>")
            )
    yield expr.alias("allele_counts")


def get_window_allele_count_expr(populations: list[str]):
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
    vcf_lf = parse_vcf_data(args.vcf_file)

    logger.info("Parsing p-values")
    pvalues_lf = parse_pvalues(args.pvalue_file)

    samples = vcf_lf.collect_schema().names()[9:]
    design_df = pl.read_csv(args.design_file)
    design_df = design_df.filter(pl.col("sample").is_in(samples))
    phenotype_to_samples = {
        d["phenotype"]: d['sample']
        for d in design_df.group_by("phenotype").agg("sample").to_dicts()
    }

    AD_index = get_position_in_format(vcf_lf, "AD")
    
    logger.info("Extracting allele depths")
    for sample in samples:
        vcf_lf = vcf_lf.with_columns(
            pl.col(sample)
            .str.split(":").list[AD_index]
            .replace(".", None)
            .str.split(",")
            .alias(f"{sample}_AD")
        )

    logger.info("Associating SNPs to windows")
    vcf_lf = add_windows(vcf_lf, window_size)

    nb_variants = get_length(vcf_lf)
    nb_pvalues = get_length(pvalues_lf)
    if nb_variants != nb_pvalues:
        raise ValueError(
            f"Number of variants ({nb_variants}) and number of pvalues ({nb_pvalues}) do not match."
        )

    lf = pl.concat([vcf_lf, pvalues_lf], how="horizontal")


    logger.info("Computing total depth")
    lf = add_total_depth(lf)

    logger.info(
        f"Computing quantile {QUANTILE} of pvalue for each pair of contig & window"
    )

    logger.info("Saving variants")
    variant_lf = (
        lf
        .with_columns(format_count_per_sample(phenotype_to_samples))
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
    variant_lf.sink_parquet(f"{args.prefix}.{VARIANTS_OUTFILE_SUFFIX}")

    logger.info("Saving windows")
    window_lf = (
        lf.group_by(["CHROM", "window"])
        .agg(
            pl.col("pvalue").quantile(QUANTILE),
            pl.col("QUAL").mean().alias("quality"),
            pl.col("total_depth").mean(),
        )
        #.with_columns(get_window_allele_count_expr(populations))
        #.with_columns(pl.concat_str(populations, separator="<br>").alias("allele_counts"))
        .rename({"CHROM": "chromosome", "window": "position"})
        .select(
            [
                "chromosome",
                "position",
                "pvalue",
                "quality",
                "total_depth",
                #"allele_counts",
            ]
        )
    )
    window_lf.sink_parquet(f"{args.prefix}.{GROUPED_VARIANTS_OUTFILE_SUFFIX}")




if __name__ == "__main__":
    main()
