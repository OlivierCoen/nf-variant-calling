import dash_mantine_components as dmc
from src.utils import style

snp_indel_selection_stack = dmc.Stack(
    [
        dmc.Slider(
            id="snp-indel-min-quality-log",
            value=0,
            color="teal",
            min=0,
            max=100,
            step=0.1,
            persistence=True,
            persisted_props=["value"],
            persistence_type="session",
            mb=35,
        ),
    ],
    align="stretch",
    gap="xl",
)

sidebar_stack = dmc.Accordion(
    value="gene_selection",
    children=[
        dmc.AccordionItem(
            [
                dmc.AccordionControl("SNP selection"),
                dmc.AccordionPanel(snp_indel_selection_stack),
            ],
            value="gene_selection",
        )
    ],
    id="sidebar-items",
    style={"marginTop": "20px", "display": "none"},
)

drawer = dmc.Drawer(
    children=[sidebar_stack],
    id="drawer",
    opened=False,
    position="right",
    withCloseButton=True,
    closeOnEscape=True,
    overlayProps=dict(backgroundOpacity=0),
    trapFocus=False,
    zIndex=10000,
    style=style.SIDEBAR,
)
