#!/usr/bin/env python3

from pathlib import Path
import argparse
import polars as pl
import pybigtools

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
        "--scores",
        type=Path,
        dest="score_file",
        required=True,
        help="Path to file containing scores",
    )
    parser.add_argument(
        "--is-neg-log",
        dest="is_neg_log",
        action="store_true",
        help="Whether scores are negative log-transformed",
    )
    parser.add_argument(
        "--fai",
        type=Path,
        dest="fai_file",
        required=True,
        help="Path to FAI file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="output_file",
        required=True,
        help="Path to output file",
    )
    return parser.parse_args()


def parse_vcf_data(vcf_file: Path) -> pl.DataFrame:
    return pl.read_csv(
        vcf_file, separator="\t", has_header=True, comment_prefix="##", low_memory=True
    )


def parse_scores(score_file: Path) -> pl.DataFrame:
    return pl.read_csv(
        score_file, has_header=False, new_columns=["score"], null_values=["NA"]
    )
    
#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():

    args = parse_args()
    
    vcf_df = parse_vcf_data(args.vcf_file)
    score_df = parse_scores(args.score_file)

    df = vcf_df.hstack(score_df)
    
    bw_df = (
        df.filter(pl.col("score").is_not_null())
        .select(
            pl.col("#CHROM").alias("chromosome"),
            (pl.col("POS") - 1).alias("start"),
            pl.col("POS").alias("end"),
            pl.col("score")
        )
        .unique()
        .sort(by=["chromosome", "start"])
    )

    if not args.is_neg_log:
        bw_df = bw_df.with_columns(
            pl.when(pl.col("score") == 0)
            .then(pl.lit(40))
            .otherwise(-pl.col("score").log10())
            .alias("score")
        )
    
    #bed_df.write_csv(BED_FILE, separator="\t", include_header=True, float_precision=6)
    
    chrom_size_df = (
        pl.scan_csv(args.fai_file, separator="\t", has_header=False)
        .rename({"column_1": "chromosome", "column_2": "size"})
        .select(["chromosome", "size"])
        .collect()
    )
    
    chrom2sizes = {d["chromosome"]: d["size"] for d in chrom_size_df.to_dicts()}

    # cannot make a context manager here
    bw = pybigtools.open(args.output_file, "w")
    vals = bw_df.select(["chromosome", "start", "end", "score"]).rows()
    bw.write(chrom2sizes, iter(vals))
    bw.close()
    

if __name__ == "__main__":
    main()
