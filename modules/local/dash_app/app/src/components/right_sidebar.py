import dash_mantine_components as dmc
from src.utils import style
from src.utils.data_management import DataManager

data_manager = DataManager()

snp_indel_min_depth = data_manager.get_min("total_depth", "snp_indel", window=False)
snp_indel_max_depth = data_manager.get_max("total_depth", "snp_indel", window=False)
snp_indel_min_quality = data_manager.get_min("quality", "snp_indel", window=False)
snp_indel_max_quality = data_manager.get_max("quality", "snp_indel", window=False)

snp_indel_min_depth_window = data_manager.get_min(
    "total_depth", "snp_indel", window=True
)
snp_indel_max_depth_window = data_manager.get_max(
    "total_depth", "snp_indel", window=True
)
snp_indel_min_quality_window = data_manager.get_min("quality", "snp_indel", window=True)
snp_indel_max_quality_window = data_manager.get_max("quality", "snp_indel", window=True)

chromosomes = data_manager.get_chromosomes("snp_indel")
snp_indel_num_chromosomes = len(chromosomes)

DEFAULT_NB_CHROMOSOMES = 100


window_filter_stack = dmc.Stack(
    [
        dmc.Text("Number of chromosomes displayed"),
        dmc.Slider(
            id="window-nb-chromosomes",
            value=DEFAULT_NB_CHROMOSOMES,
            color="orange",
            min=1,
            max=snp_indel_num_chromosomes,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text("Range of quality"),
        dmc.RangeSlider(
            id="window-quality-range",
            value=[snp_indel_min_quality_window, snp_indel_max_quality_window],
            color="teal",
            min=snp_indel_min_quality_window,
            max=snp_indel_max_quality_window,
            step=1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
        dmc.Text("Range of total depth"),
        dmc.RangeSlider(
            id="window-depth-range",
            value=[snp_indel_min_depth_window, snp_indel_max_depth_window],
            color="blue",
            min=snp_indel_min_depth_window,
            max=snp_indel_max_depth_window,
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

variant_filter_stack = dmc.Stack(
    [
        dmc.Text("Chromosomes"),
        dmc.MultiSelect(
            id="variant-chromosomes",
            value=chromosomes[0],
            data=chromosomes,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            w=400,
            mb=10,
        ),
        dmc.Text("Range of quality"),
        dmc.RangeSlider(
            id="variant-quality-range",
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
            id="variant-depth-range",
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
    ],
    align="stretch",
    gap="xl",
)

options = dmc.Stack(
    [
        dmc.Accordion(
            value="filters",
            children=[
                dmc.AccordionItem(
                    [
                        dmc.AccordionControl("Sliding windows"),
                        dmc.AccordionPanel(window_filter_stack),
                        dmc.Button("Update", id="window-update-button"),
                    ],
                    value="windows",
                ),
                dmc.AccordionItem(
                    [
                        dmc.AccordionControl("Variants"),
                        dmc.AccordionPanel(variant_filter_stack),
                        dmc.Button("Update", id="variant-update-button"),
                    ],
                    value="variants",
                ),
            ],
            id="sidebar-items",
            style=style.SIDEBAR,
            multiple=True,
        ),
    ],
    align="stretch",
    gap="xl",
)
