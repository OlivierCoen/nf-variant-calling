import logging
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from src.utils import config

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


@lru_cache(maxsize=None)
class DataManager:
    def __init__(self):
        folder = Path(config.DATA_FOLDER)
        # compute other static values
        self.variants = {}
        self.quantiles = {}
        for variant_type in config.VARIANT_TYPES:
            logger.info(f"Preparing {variant_type} data")

            variant_file = folder / f"{variant_type}.{config.INPUT_FILE_SUFFIX}"
            self.variants[variant_type] = self.prepare_data(variant_file)

            self.quantiles[variant_type] = self.get_quantiles(variant_type)

            logger.info(f"{variant_type} data loaded")

    @staticmethod
    def parse_parquet(file: Path) -> pl.LazyFrame | None:
        logger.info(f"Parsing {file} file")
        if not Path(file).is_file():
            return None
        return pl.scan_parquet(file)

    def get_max(self, value: str, variant_type: str) -> int:
        if self.variants[variant_type] is None:
            return 0
        return (
            self.variants[variant_type]
            .select(pl.col(value).cast(pl.Int64).max().alias("max_val"))
            .collect()
            .to_series()
            .item(0)
        )

    def get_min(self, value: str, variant_type: str) -> int:
        if self.variants[variant_type] is None:
            return 0
        return (
            self.variants[variant_type]
            .select(pl.col(value).cast(pl.Int64).min().alias("min_val"))
            .collect()
            .to_series()
            .item(0)
        )

    def get_sorted_chromosomes(self, variant_type: str) -> pl.LazyFrame:
        if self.variants[variant_type] is None:
            return pl.LazyFrame(schema=["chromosome"])

        return (
            self.variants[variant_type]
            .group_by("chromosome")
            .agg(pl.col("position").max().alias("max_pos"))
            .sort("max_pos", descending=True)
        )

    def get_chromosomes(self, variant_type: str) -> list[str]:
        sorted_chromosomes = self.get_sorted_chromosomes(variant_type)
        if sorted_chromosomes is None:
            return []
        return (
            sorted_chromosomes.select(pl.col("chromosome"))
            .collect()
            .to_series()
            .to_list()
        )

    @staticmethod
    def prepare_data(file: Path) -> pl.LazyFrame | None:
        lf = DataManager.parse_parquet(file)
        if lf is None:
            return None

        return lf.select(
            pl.col(
                [
                    "chromosome",
                    "position",
                    "pvalue",
                    "quality",
                    "total_depth",
                    "allele_counts",
                ]
            ),
            (pl.col("chromosome") + "_" + pl.col("position").cast(pl.String)).alias(
                "snp"
            ),
            (pl.col("chromosome") + "_" + pl.col("position").cast(pl.String)).alias(
                "gene"
            ),
        )

    def get_quantiles(self, variant_type: str) -> dict[float, float]:
        if self.variants[variant_type] is None:
            return {}

        pvalues = (
            self.variants[variant_type]
            .select("pvalue")
            .filter(pl.col("pvalue").is_not_nan())
            .collect()
            .to_series()
            .to_numpy()
        )

        quantile_options = np.arange(
            0, 1.01, 0.01
        )  # 101 elements from 0 to 1 inclusive

        return {
            (i / 100): float(quantile)
            for i, quantile in enumerate(np.quantile(pvalues, quantile_options))
        }

    def get_manhattanplot_data(
        self,
        variant_type: str,
        pvalue_quantile_range: list[float],
        quality_range: list[int],
        depth_range: list[int],
        chromosome: str,
    ) -> pd.DataFrame:
        if self.variants[variant_type] is None:
            return pd.DataFrame()

        min_pvalue = self.quantiles[variant_type][pvalue_quantile_range[0]]
        max_pvalue = self.quantiles[variant_type][pvalue_quantile_range[1]]

        return (
            self.variants[variant_type]
            .filter(pl.col("chromosome") == chromosome)
            .filter(pl.col("pvalue").is_between(min_pvalue, max_pvalue))
            .filter(
                pl.col("quality").is_between(quality_range[0], quality_range[1] + 1)
            )
            .filter(
                pl.col("total_depth").is_between(depth_range[0], depth_range[1] + 1)
            )
            .collect()
            .to_pandas()
        )
