import dash_bio
from dash_extensions.enrich import Input, Output, State, Trigger, callback
from src.utils.data_management import DataManager

data_manager = DataManager()

##############################################
##############################################
# CALLBACKS
##############################################
##############################################

manhattan_plot_kwargs = dict(
    chrm="chromosome",
    bp="position",
    p="cmh_pvalue",
    snp="snp",
    gene="gene",
    annotation="annotation",
    logp=True,
    suggestiveline_color="#AA00AA",
    genomewideline_color="#AA5500",
)


def register_callbacks():
    @callback(
        Output("drawer", "opened"),
        Trigger("settings-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def open_drawer():
        return True

    @callback(
        Output("snp-indel-graph", "figure"),
        Output("sv-graph", "figure"),
        Output("snp-indel-graph", "style"),
        Output("sv-graph", "style"),
        Input("snp-indel-min-quality-log", "value"),
        State("snp-indel-graph", "style"),
        State("sv-graph", "style"),
    )
    def update_graph_data(
        min_quality_log: float, snp_indel_graph_style: dict, sv_graph_style: dict
    ):
        figures = {}
        styles = {"snp_indel": snp_indel_graph_style, "sv": sv_graph_style}

        for data_type in ["snp_indel", "sv"]:
            df = data_manager.get_manhattanplot_data(data_type, min_quality_log)

            if df.empty:
                fig = {}
                styles[data_type]["display"] = "none"
            else:
                fig = dash_bio.ManhattanPlot(
                    dataframe=df,
                    **manhattan_plot_kwargs,
                )
                styles[data_type]["display"] = "block"

            figures[data_type] = fig

        return tuple(list(figures.values()) + list(styles.values()))
