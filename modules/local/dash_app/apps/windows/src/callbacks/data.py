import logging

from dash_extensions.enrich import Input, Output, State, Trigger, callback
from src.components.dash_bio.manhattan import ManhattanPlot
from src.utils import config
from src.utils.data_management import DataManager

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

data_manager = DataManager()

##############################################
##############################################
# CALLBACKS
##############################################
##############################################

manhattan_plot_kwargs = dict(
    logp=True,
    suggestiveline_color="#AA00AA",
    genomewideline_color="#AA5500",
    highlight_color="#4F2196",
)


def register_callbacks():

    @callback(
        Output("snp-indel-graph", "figure"),
        Output("sv-graph", "figure"),
        Trigger("update-button", "n_clicks"),
        State("quality-range", "value"),
        State("depth-range", "value"),
        State("nb-chromosomes", "value"),
        # background=True,
    )
    def update_graph_data(
        quality_range: list[int],
        depth_range: list[int],
        nb_chromosomes: int,
    ):
        figures = {}

        for variant_type in config.VARIANT_TYPES:
            if variant_type in data_manager.parsed_variant_types:
                logger.info(f"Updating {variant_type} graph data")
                df = data_manager.get_manhattanplot_data(
                    variant_type, quality_range, depth_range, nb_chromosomes
                )
                if df.empty:
                    fig = {}
                else:
                    fig = ManhattanPlot(
                        dataframe=df,
                        **manhattan_plot_kwargs,
                    )
                figures[variant_type] = fig
            else:
                figures[variant_type] = {}

        return tuple(list(figures.values()))
