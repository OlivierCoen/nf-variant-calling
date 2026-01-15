from functools import lru_cache
from pathlib import Path

import pandas as pd
import polars as pl
from src.utils import config


@lru_cache(maxsize=None)
class DataManager:
    def __init__(self):
        # parse initial data
        raw_data = dict(
            snp_indel=dict(
                variants=self.parse_vcf_data(config.SNP_INDEL["vcf_file"]),
                pvalues=self.parse_pvalues(config.SNP_INDEL["pvalues_file"]),
            ),
            sv=dict(
                variants=self.parse_vcf_data(config.SV["vcf_file"]),
                pvalues=self.parse_pvalues(config.SV["pvalues_file"]),
            ),
        )
        # compute other static values
        self.data = {}
        for key in raw_data:
            self.data[key] = self.prepare_data(
                raw_data[key]["variants"], raw_data[key]["pvalues"]
            )

    @staticmethod
    def parse_vcf_data(vcf_file: str) -> pl.LazyFrame | None:
        if not Path(vcf_file).is_file():
            return None
        return pl.scan_csv(
            vcf_file, separator="\t", has_header=True, comment_prefix="##"
        ).rename({"#CHROM": "CHROM"})

    @staticmethod
    def parse_pvalues(cmh_pvalues_file: str) -> list[float | None] | None:
        if not Path(cmh_pvalues_file).is_file():
            return None
        pvalues = []
        with open(cmh_pvalues_file, "r") as fin:
            for line in fin.readlines():
                val = line.strip()
                if val == "NA":
                    pvalues.append(None)
                else:
                    pvalues.append(float(val))
        return pvalues

    @staticmethod
    def get_position_in_format(vcf_lf: pl.LazyFrame, count_type: str) -> int:
        fmt_series = vcf_lf.select("FORMAT").collect().to_series()
        if fmt_series.unique().len() > 1:
            raise ValueError(f"More than one format found: {fmt_series.unique()}")
        fmt = fmt_series.unique().item()
        return fmt.split(":").index(count_type)

    @staticmethod
    def extract_counts(vcf_lf: pl.LazyFrame, count_type: str) -> pl.DataFrame:
        sample_cols = vcf_lf.collect_schema().names()[9:]
        count_type_position = DataManager.get_position_in_format(vcf_lf, count_type)

        count_df = vcf_lf.select(
            pl.col(sample_cols)
            .str.split(":")
            .list[count_type_position]
            .replace(".", None)
            .cast(pl.String())
        ).collect()

        return count_df.select(pl.col(col) for col in sample_cols)

    @staticmethod
    def get_annotation(vcf_lf: pl.LazyFrame) -> pl.Series:
        RO_df = DataManager.extract_counts(vcf_lf, "RO")
        AO_df = DataManager.extract_counts(vcf_lf, "AO")
        df = RO_df + "/" + AO_df
        for col in df.columns:
            df = df.with_columns(
                pl.when(pl.col(col).is_not_null())
                .then(pl.lit(col) + ": " + pl.col(col))
                .otherwise(pl.col(col))
                .alias(col)
            )
        return df.select(
            pl.concat_str(pl.all(), separator="\n", ignore_nulls=True)
        ).to_series()

    @staticmethod
    def get_chrom_ame_mapping(vcf_lf: pl.LazyFrame) -> dict[str, int]:
        unique_chroms = (
            vcf_lf.select("CHROM").collect().unique().sort("CHROM").to_series()
        )
        return {chrom: i + 1 for i, chrom in enumerate(unique_chroms)}

    @staticmethod
    def prepare_data(
        vcf_lf: pl.LazyFrame | None, pvalues: list[float | None]
    ) -> pl.LazyFrame | None:
        if vcf_lf is None:
            return None

        chrom_to_index = DataManager.get_chrom_ame_mapping(vcf_lf)
        annot_series = DataManager.get_annotation(vcf_lf)

        return vcf_lf.select(
            pl.col("CHROM").replace(chrom_to_index).cast(pl.Int64).alias("chromosome"),
            pl.col("POS").cast(pl.Int64).alias("position"),
            pl.col("QUAL").cast(pl.Float64).alias("quality"),
            (pl.col("CHROM") + "_" + pl.col("POS").cast(pl.String)).alias("snp"),
            (pl.col("CHROM") + "_" + pl.col("POS").cast(pl.String)).alias("gene"),
            pl.Series(pvalues).cast(pl.Float64).alias("cmh_pvalue"),
            annot_series.cast(pl.String).alias("annotation"),
        ).filter(pl.col("cmh_pvalue").is_not_null())

    def get_manhattanplot_data(
        self, data_type: str, min_quality: float
    ) -> pd.DataFrame:
        if self.data[data_type] is None:
            return pd.DataFrame()
        return (
            self.data[data_type]
            .filter(pl.col("quality") >= min_quality)
            .collect()
            .to_pandas()
        )
