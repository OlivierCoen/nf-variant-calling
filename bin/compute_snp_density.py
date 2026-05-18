#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

import polars as pl
from tqdm import tqdm
from common import parse_vcf_data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
    parser.add_argument(
        "--window-size",
        dest="window_size",
        type=int,
        required=True,
        help="Window size",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="outfile",
        required=True,
        help="Path to output file",
    )
    return parser.parse_args()

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():
    args = parse_args()

    logger.info("Parsing VCF file")

    vcf_lf = parse_vcf_data(args.vcf_file)

    contigs = (
        vcf_lf.select("CHROM")
        .collect()
        .unique(maintain_order=True)
        .to_series().to_list()
    )

    logger.info(f"Computing SNP density for {len(contigs)} contigs")
    for contig in tqdm(contigs):
        contig_lf = vcf_lf.filter(pl.col("CHROM") == contig)
        max_position = contig_lf.select("POS").collect().max().item()
        positions = list(range(1, max_position + 1))
        
        all_positions_df = pl.DataFrame({"position": positions})
        snp_positions_df = contig_lf.with_columns(pl.lit(1).alias("count")).select(
            pl.col("POS").alias("position"), 
            pl.col("count")
        ).collect()
    
        #weights = [2 * (1 + i) / args.window_size for i in range(args.window_size // 2)]
        #weights = weights + list(reversed(weights))
        
        scores = (
            all_positions_df
                .join(snp_positions_df, on="position", how="left")
                .with_columns(pl.col("count").alias("original_count"))
                .with_columns(pl.col("count").fill_null(0))
                .with_columns(
                    pl.col("count").rolling_mean(
                        window_size=args.window_size, 
                        min_samples=1, 
                        center=True, 
                        #weights=weights
                    ).alias("density")
                )
                .filter(pl.col("original_count").is_not_null())
                .select((1 - 300 * pl.col("density")).alias("score"))
                #.drop_nulls()
                .to_series().to_list()
        )

        scores = [max(1e-10, score) for score in scores]

        with open(args.outfile, 'a') as fout:
            fout.writelines([f'{score:.10f}\n' for score in scores])


if __name__ == "__main__":
    main()
