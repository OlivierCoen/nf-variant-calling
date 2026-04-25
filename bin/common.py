import gzip
from pathlib import Path

import polars as pl


def parse_vcf_header(vcf_file: Path) -> list:
    header_lines = []
    with gzip.open(vcf_file, "rb") as fin:
        for line in fin:
            line = line.decode("utf-8")
            if line.startswith("##"):
                header_lines.append(line)
            else:
                break
    return header_lines


def parse_vcf_data(vcf_file: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        vcf_file, separator="\t", has_header=True, comment_prefix="##", low_memory=True
    ).rename({"#CHROM": "CHROM"})


def get_position_in_format(vcf_lf: pl.LazyFrame, info: str) -> int:
    fmt_series = vcf_lf.select("FORMAT").collect().to_series()
    if fmt_series.unique().len() > 1:
        raise ValueError(f"More than one format found: {fmt_series.unique()}")
    fmt = fmt_series.unique().item()
    return fmt.split(":").index(info)


def extract_counts(
    vcf_lf: pl.LazyFrame, sample_cols: list[str], idx: int
) -> pl.LazyFrame:
    return vcf_lf.select(
        pl.col(sample_cols).str.split(":").list[idx].replace(".", None).cast(pl.Int64)
    )
