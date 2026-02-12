import dash_mantine_components as dmc
from dash_iconify import DashIconify
from src.components import graphs
from src.utils import style

gene_icon = DashIconify(icon="material-symbols:genetics", width=20)

sample_icon = DashIconify(icon="ic:baseline-dashboard-customize", width=20)


variant_tabs = dmc.Tabs(
    children=[
        dmc.TabsList(
            children=[
                dmc.TabsTab(
                    dmc.Text("SNPs / Indels", fw=500),
                    value="snp_indel",
                    color="teal",
                    style=style.HEADER_TABLIST_ITEM,
                ),
                dmc.TabsTab(
                    dmc.Text("Structural variants", fw=500),
                    value="sv",
                    color="red",
                    style=style.HEADER_TABLIST_ITEM,
                ),
            ],
            style=style.HEADER_TABLIST,
        ),
        dmc.TabsPanel(
            children=[graphs.snp_indel_graph],
            style=style.TABS_PANEL,
            value="snp_indel",
        ),
        dmc.TabsPanel(
            children=[graphs.sv_graph],
            style=style.TABS_PANEL,
            value="sv",
        ),
    ],
    id="tabs",
    variant="default",
    radius="md",
    orientation="horizontal",
    placement="right",
    value="snp_indel",
    persistence=True,
    persisted_props=["value"],
    persistence_type="session",
    style=style.TAB,
)
