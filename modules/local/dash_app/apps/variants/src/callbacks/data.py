import logging

from dash_extensions.enrich import Input, Output, State, Trigger, callback
from src.components.dash_bio.manhattan import ManhattanPlot
from src.utils.data_management import DataManager
from src.utils import config

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
        State("pvalue-quantile-range", "value"),
        State("quality-range", "value"),
        State("depth-range", "value"),
        State("chromosome", "value"),
        # background=True,
    )
    def update_graph_data(
        pvalue_quantile_range: list[float],
        quality_range: list[int],
        depth_range: list[int],
        chromosome: str,
    ):
        figures = {}

        for variant_type in config.VARIANT_TYPES:
            logger.info(f"Updating {variant_type} graph data")
            df = data_manager.get_manhattanplot_data(
                variant_type,
                pvalue_quantile_range,
                quality_range,
                depth_range,
                chromosome,
            )
            if df.empty:
                fig = {}
            else:
                fig = ManhattanPlot(
                    dataframe=df,
                    **manhattan_plot_kwargs,
                )
            figures[variant_type] = fig

        return tuple(list(figures.values()))
