from dash_extensions.enrich import dcc
from src.utils import style


def get_graph(graph_id: str):
    return dcc.Graph(id=graph_id, figure={}, style=style.GRAPH)


snp_indel_graph = get_graph("snp-indel-graph")
snp_indel_window_graph = get_graph("snp-indel-window-graph")

sv_graph = get_graph("sv-graph")
sv_window_graph = get_graph("sv-window-graph")
