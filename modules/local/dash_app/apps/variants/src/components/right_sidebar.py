import dash_mantine_components as dmc
from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

snp_indel_min_depth = data_manager.get_min("total_depth", "snp_indel")
snp_indel_max_depth = data_manager.get_max("total_depth", "snp_indel")
snp_indel_min_quality = data_manager.get_min("quality", "snp_indel")
snp_indel_max_quality = data_manager.get_max("quality", "snp_indel")

chromosomes = data_manager.get_chromosomes("snp_indel")
snp_indel_num_chromosomes = len(chromosomes)
default_chromosome = chromosomes[0] if chromosomes else None

DEFAULT_NB_CHROMOSOMES = 100

options = dmc.Stack(
    [
        dmc.Select(
            label="Chromosome",
            placeholder="Select chromosome to display",
            id="chromosome",
            value=default_chromosome,
            data=chromosomes,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
        ),
        dmc.Text("Range of pvalue quantiles"),
        dmc.RangeSlider(
            id="pvalue-quantile-range",
            value=[0, 1],
            color="purple",
            min=0,
            max=1,
            step=0.05,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text("Range of quality"),
        dmc.RangeSlider(
            id="quality-range",
            value=[snp_indel_min_quality, snp_indel_max_quality],
            color="teal",
            min=snp_indel_min_quality,
            max=snp_indel_max_quality,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text("Range of depth"),
        dmc.RangeSlider(
            id="depth-range",
            value=[snp_indel_min_depth, snp_indel_max_depth],
            color="blue",
            min=snp_indel_min_depth,
            max=snp_indel_max_depth,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Button("Update", id="update-button"),
    ],
    align="stretch",
    gap="xl",
    style=style.SIDEBAR,
)
