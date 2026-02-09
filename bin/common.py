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
