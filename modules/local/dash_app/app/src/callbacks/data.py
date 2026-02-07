from dash_extensions.enrich import Input, Output, State, Trigger, callback

from src.components.dash_bio.manhattan import ManhattanPlot
from src.utils.data_management import DataManager

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
        Input("snp-indel-range-quality", "value"),
        Input("snp-indel-range-depth", "value"),
    )
    def update_graph_data(quality_range: list[int], depth_range: list[int]):
        figures = {}

        for data_type in ["snp_indel", "sv"]:
            df = data_manager.get_manhattanplot_data(
                data_type, quality_range, depth_range
            )
            print(df.columns)
            if df.empty:
                fig = {}
            else:
                fig = ManhattanPlot(
                    dataframe=df,
                    **manhattan_plot_kwargs,
                )

            figures[data_type] = fig

        return tuple(list(figures.values()))
