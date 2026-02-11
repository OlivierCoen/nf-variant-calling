import socket

import dash_mantine_components as dmc
import diskcache
from dash import DiskcacheManager
from dash_extensions.enrich import DashProxy, TriggerTransform, html
from src.callbacks import data
from src.components import right_sidebar, stores, tabs, tooltips
from src.utils import config, style

DEBUG = True
# DEBUG = False


cache = diskcache.Cache("./cache")
background_callback_manager = DiskcacheManager(cache)

app = DashProxy(
    __name__,
    title=config.APP_TITLE,
    prevent_initial_callbacks="initial_duplicate",
    suppress_callback_exceptions=(not DEBUG),
    update_title=config.UPDATE_TITLE,
    external_stylesheets=[dmc.styles.ALL],
    background_callback_manager=background_callback_manager,
    transforms=[TriggerTransform()],
)

# -------------------- LAYOUT --------------------


def serve_layout():
    return dmc.MantineProvider(
        children=[
            html.Div(
                [
                    right_sidebar.options,
                    tabs.variant_tabs,
                    *stores.stores_to_load,
                    *tooltips.tooltips_to_load,
                ],
                id="layout",
                style=style.LAYOUT,
            )
        ]
    )


app.layout = serve_layout

# -------------------- IMPORTING CALLBACKS --------------------

data.register_callbacks()

# -------------------- LAUNCH SERVER --------------------


def find_port(port: int) -> int:
    """Find a port not in use starting at given port"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(("localhost", port)) == 0:
            return find_port(port=port + 1)
        else:
            return port


if __name__ == "__main__":
    # setting prune_errors to False avoids error message pruning
    # in order to get original tracebacks
    # (very useful for debugging)
    prune_errors = False if DEBUG else True
    app.run(
        debug=DEBUG,
        host=config.HOST,
        port=find_port(port=config.PLOTLY_APP_PORT),
        dev_tools_prune_errors=prune_errors,
    )
