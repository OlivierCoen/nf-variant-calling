import logging

from dash_extensions.enrich import Input, Output, State, Trigger, callback
from src.components.dash_bio.manhattan import ManhattanPlot
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
)


def register_callbacks():

    @callback(
        Output("snp-indel-graph", "figure"),
        Output("sv-graph", "figure"),
        Trigger("variant-update-button", "n_clicks"),
        State("variant-quality-range", "value"),
        State("variant-depth-range", "value"),
        State("variant-chromosomes", "value"),
        # background=True,
    )
    def update_variant_graph_data(
        quality_range: list[int],
        depth_range: list[int],
        chromosomes: str | list[str],
    ):
        figures = {}
        if isinstance(chromosomes, str):
            chromosomes = [chromosomes]

        for data_type in ["snp_indel", "sv"]:
            logger.info(f"Updating {data_type} graph data")
            df = data_manager.get_variant_manhattanplot_data(
                data_type, quality_range, depth_range, chromosomes
            )
            if df.empty:
                fig = {}
            else:
                fig = ManhattanPlot(
                    dataframe=df,
                    **manhattan_plot_kwargs,
                )
            figures[data_type] = fig

        return tuple(list(figures.values()))

    @callback(
        Output("snp-indel-window-graph", "figure"),
        Output("sv-window-graph", "figure"),
        Trigger("window-update-button", "n_clicks"),
        State("window-quality-range", "value"),
        State("window-depth-range", "value"),
        State("window-nb-chromosomes", "value"),
        # background=True,
    )
    def update_window_graph_data(
        quality_range: list[int],
        depth_range: list[int],
        nb_chromosomes: int,
    ):
        figures = {}

        for data_type in ["snp_indel", "sv"]:
            logger.info(f"Updating {data_type} graph data")
            df = data_manager.get_window_manhattanplot_data(
                data_type, quality_range, depth_range, nb_chromosomes
            )
            if df.empty:
                fig = {}
            else:
                fig = ManhattanPlot(
                    dataframe=df,
                    **manhattan_plot_kwargs,
                )
            figures[data_type] = fig

        return tuple(list(figures.values()))
