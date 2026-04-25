import argparse
import logging
from pathlib import Path
import polars as pl

from scipy.stats import fisher_exact
import numpy as np
from scipy.stats import MonteCarloMethod

from common import parse_vcf_data

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

def parse_args():
    parser = argparse.ArgumentParser()
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
    return parser.parse_args()


def get_position_in_format(vcf_lf: pl.LazyFrame, info: str) -> int:
    fmt_series = vcf_lf.select("FORMAT").collect().to_series()
    if fmt_series.unique().len() > 1:
        raise ValueError(f"More than one format found: {fmt_series.unique()}")
    fmt = fmt_series.unique().item()
    return fmt.split(":").index(info)


fisher_exact([[8, 2, 3], [1, 5, 4]], method=method)

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

def main():
    args = parse_args()
    
    logger.info("Parsing allele counts")
    RO_lf = pl.scan_parquet(args.RO_file).select(pl.all().name.suffix("_RO"))
    AO_lf = pl.scan_parquet(args.AO_file).select(pl.all().name.suffix("_AO"))

    RO_cols = RO_lf.collect_schema().names()
    AO_cols = AO_lf.collect_schema().names()
    
    design_df = pl.read_csv(args.design_file)
    
    rng = np.random.default_rng(seed=42)
    method = MonteCarloMethod(rng=rng)
    
    