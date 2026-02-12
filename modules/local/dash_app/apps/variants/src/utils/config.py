PLOTLY_APP_PORT = 8080
HOST = "127.0.0.1"

LOGGING_FORMAT = "[%(asctime)s] [%(name)s] %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d_%H-%M-%S"

APP_TITLE = "Variant analysis"
UPDATE_TITLE = "Updating ..."

DATA_FOLDER = "data"

VARIANT_TYPES = ["snp_indel", "sv"]
INPUT_FILE_SUFFIX = "formated_variants.parquet"


AG_GRID_DEFAULT_COLUMN_DEF = {
    "filter": True,
    "resizable": True,
    "editable": False,
    "sortable": True,
}

AG_GRID_DEFAULT_OPTIONS = {"pagination": True, "paginationAutoPageSize": True}
