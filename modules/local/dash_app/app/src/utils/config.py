PLOTLY_APP_PORT = 8080
HOST = "0.0.0.0"

LOGGING_FORMAT = "[%(asctime)s] [%(name)s] %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d_%H-%M-%S"

APP_TITLE = "Variant analysis"
UPDATE_TITLE = "Updating ..."

DATA_FOLDER = "data"

SNP_INDEL = {
    "vcf_file": "data/snp_indel.filtered.vcf.gz",
    "pvalues_file": "data/snp_indel.cmh_pvalues.txt",
}


SV = {"vcf_file": "data/sv.filtered.vcf.gz", "pvalues_file": "data/sv.cmh_pvalues.txt"}


AG_GRID_DEFAULT_COLUMN_DEF = {
    "filter": True,
    "resizable": True,
    "editable": False,
    "sortable": True,
}

AG_GRID_DEFAULT_OPTIONS = {"pagination": True, "paginationAutoPageSize": True}
