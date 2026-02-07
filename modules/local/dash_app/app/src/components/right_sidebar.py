import dash_mantine_components as dmc

from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

snp_indel_min_depth = data_manager.get_min("total_depth", "snp_indel")
snp_indel_max_depth = data_manager.get_max("total_depth", "snp_indel")
snp_indel_min_quality = data_manager.get_min("quality", "snp_indel")
snp_indel_max_quality = data_manager.get_max("quality", "snp_indel")


filter_stack = dmc.Stack(
    [
        dmc.Text("Minimum quality"),
        dmc.RangeSlider(
            id="snp-indel-range-quality",
            value=[snp_indel_min_depth, snp_indel_max_depth],
            color="teal",
            min=snp_indel_min_depth,
            max=snp_indel_max_depth,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text("Range of total depth"),
        dmc.RangeSlider(
            id="snp-indel-range-depth",
            value=[snp_indel_min_quality, snp_indel_max_quality],
            color="blue",
            min=snp_indel_min_quality,
            max=snp_indel_max_quality,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
    ],
    align="stretch",
    gap="xl",
)

options = dmc.Accordion(
    value="filters",
    children=[
        dmc.AccordionItem(
            [
                dmc.AccordionControl("Variant filtering"),
                dmc.AccordionPanel(filter_stack),
            ],
            value="filters",
        )
    ],
    id="sidebar-items",
    style=style.SIDEBAR,
)
