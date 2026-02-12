import logging

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pandas.api.types import is_numeric_dtype
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

SUGGESTIVE_LINE_LABEL = "suggestive line"
GENOMEWIDE_LINE_LABEL = "genomewide line"


def ManhattanPlot(
    dataframe,
    chrm="chromosome",
    bp="position",
    p="pvalue",
    snp="snp",
    gene="gene",
    quality="quality",
    total_depth="total_depth",
    allele_counts="allele_counts",
    logp=True,
    title="Manhattan Plot",
    showgrid=True,
    xlabel=None,
    ylabel="-log10(p)",
    point_size=5,
    showlegend=True,
    col=None,
    suggestiveline_value=-np.log10(1e-8),
    suggestiveline_color="#636efa",
    suggestiveline_width=1,
    genomewideline_value=-np.log10(5e-8),
    genomewideline_color="#EF553B",
    genomewideline_width=1,
    highlight=True,
    highlight_color="black",
):
    mh = _ManhattanPlot(
        dataframe,
        chrm=chrm,
        bp=bp,
        p=p,
        snp=snp,
        gene=gene,
        quality=quality,
        total_depth=total_depth,
        allele_counts=allele_counts,
        logp=logp,
    )

    return mh.figure(
        title=title,
        showgrid=showgrid,
        xlabel=xlabel,
        ylabel=ylabel,
        point_size=point_size,
        showlegend=showlegend,
        col=col,
        suggestiveline_value=suggestiveline_value,
        suggestiveline_color=suggestiveline_color,
        suggestiveline_width=suggestiveline_width,
        genomewideline_value=genomewideline_value,
        genomewideline_color=genomewideline_color,
        genomewideline_width=genomewideline_width,
        highlight=highlight,
        highlight_color=highlight_color,
    )


class _ManhattanPlot:
    def __init__(
        self,
        x,
        chrm="chromosome",
        bp="position",
        p="pvalue",
        snp="snp",
        gene="gene",
        quality="quality",
        total_depth="total_depth",
        allele_counts="allele_counts",
        logp=True,
    ):
        # checking the validity of the arguments

        # Make sure you have chrm, bp and p columns and that they are of
        # numeric type
        if chrm not in x.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % chrm)

        if bp not in x.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % bp)
        else:
            if not is_numeric_dtype(x[bp].dtype):
                raise TypeError("%s column should be numeric type" % bp)

        if p not in x.columns.values:
            raise KeyError("Column %s not found in 'x' data.frame" % p)
        else:
            if not is_numeric_dtype(x[p].dtype):
                raise TypeError("%s column should be numeric type" % p)

        # Create a new DataFrame with columns named after chrm, bp, and p.
        self.data = x

        # approximating chromosome size byt the max value of its positions
        chrom_size_dict = (
            self.data.groupby("chromosome")["position"].agg("max").to_dict()
        )
        self.data["chrom_size"] = self.data["chromosome"].map(chrom_size_dict)

        # sorting datarame to have chromosomes in descending order and positions in ascending order
        self.data = self.data.sort_values(
            by=["chrom_size", bp], ascending=[False, True]
        ).drop(columns=["chrom_size"])

        # sorting chrom_size_dict
        chrom_size_dict = {
            k: v
            for k, v in sorted(
                chrom_size_dict.items(), key=lambda item: item[1], reverse=True
            )
        }

        # making dict of cumulative sizes (considering all previous chromosomes) for each chromosome
        self.cumulative_previous_sizes = {}
        self.cumulative_sizes = {}
        current_pos = 0
        for chrom, size in chrom_size_dict.items():
            self.cumulative_previous_sizes[chrom] = current_pos
            current_pos += size
            self.cumulative_sizes[chrom] = current_pos

        self.xlabel = ""
        self.ticks = []
        self.ticksLabels = []
        self.chrom_colname = chrm
        self.position_colname = bp
        self.pvalue_colname = p
        self.snp_colname = snp
        self.gene_colname = gene
        self.quality_colname = quality
        self.total_depth_colname = total_depth
        self.allele_counts_colname = allele_counts
        self.logp = logp

        # This section sets up positions and ticks. Ticks should be placed in
        # the middle of a chromosome. The new pos column is added that keeps
        # a running sum of the positions of each successive chromosome.
        # For example:
        # chrm bp pos
        # 1   1  1
        # 1   2  2
        # 2   1  3
        # 2   2  4
        # 3   1  5

        self.chromosomes = list(chrom_size_dict.keys())

        self.data["cumulative_position"] = self.data.apply(
            lambda row: (
                self.cumulative_previous_sizes[row["chromosome"]] + row["position"]
            ),
            axis=1,
        )

        # sorting data by cumulative_position
        self.data.sort_values(by=["cumulative_position"], ascending=True, inplace=True)

        tick_mininums = (
            self.data.groupby("chromosome", sort=False)["cumulative_position"]
            .apply(lambda x: x.iloc[0])
            .to_dict()
        )
        tick_maxinums = (
            self.data.groupby("chromosome", sort=False)["cumulative_position"]
            .apply(lambda x: x.iloc[-1])
            .to_dict()
        )

        self.xlabel = "Chromosomes"

        for tick_min, tick_max in zip(tick_mininums.values(), tick_maxinums.values()):
            self.ticks.append(int((tick_min + tick_max) / 2) + 1)

        self.ticksLabels = self.chromosomes

    def get_hover_text(self, df: pd.DataFrame):
        return (
            "SNP: "
            + df[self.snp_colname].astype(str)
            + "<br>POSITION: "
            + df[self.position_colname].astype(str)
            + "<br>GENE: "
            + df[self.gene_colname].astype(str)
            + "<br>P-VALUE (QUANTILE 5%): "
            + df[self.pvalue_colname].astype(str)
            + "<br>MEAN QUALITY: "
            + df[self.quality_colname].astype(str)
            + "<br>MEAN TOTAL DEPTH: "
            + df[self.total_depth_colname].astype(str)
            + "<br><br>ALLELE COUNTS (REF / ALT):<br>"
            + df[self.allele_counts_colname].astype(str)
        )

    def get_yvalues(self, data: pd.DataFrame):
        return (
            -np.log10(data[self.pvalue_colname].values)
            if self.logp
            else data[self.pvalue_colname].values
        )

    def figure(
        self,
        title="Manhattan Plot",
        showgrid=True,
        xlabel=None,
        ylabel="-log10(p)",
        point_size=5,
        showlegend=True,
        col=None,
        suggestiveline_value=-np.log10(1e-8),
        suggestiveline_color="blue",
        suggestiveline_width=1,
        genomewideline_value=-np.log10(5e-8),
        genomewideline_color="red",
        genomewideline_width=1,
        highlight=True,
        highlight_color="black",
    ):
        """Keyword arguments:
        - title (string; default 'Manhattan Plot'): The title of the
            graph.
        - showgrid (bool; default True): Boolean indicating whether
            gridlines should be shown.
        - xlabel (string; optional): Label of the x axis.
        - ylabel (string; default '-log10(p)'): Label of the y axis.
        - point_size (number; default 5): Size of the points of the
            scatter plot.
        - showlegend (bool; default True): Boolean indicating whether
            legends should be shown.
        - col (string; optional): A string representing the color of the
            points of the Scatter plot. Can be in any color format
            accepted by plotly.graph_objects.
        - suggestiveline_value (bool | float; default 8): A value which
            must be either False to deactivate the option, or a numerical value
            corresponding to the p-value at which the line should be
            drawn. The line has no influence on the data points.
        - suggestiveline_color (string; default 'grey'): Color of the
            suggestive line.
        - suggestiveline_width (number; default 2): Width of the
            suggestive line.
        - genomewideline_value (bool | float; default -log10(5e-8)): A
            boolean which must be either False to deactivate the option, or a
            numerical value corresponding to the p-value above which the
            data points are considered significant.
        - genomewideline_color (string; default 'red'): Color of the
            genome-wide line. Can be in any color format accepted by
            plotly.graph_objects.
        - genomewideline_width (number; default 1): Width of the genome
          wide line.
        - highlight (bool; default True): Whether to turn on or off the
            highlighting of data points considered significant.

        Returns:
        - A figure formatted for plotly.graph_objects.

        """

        xmin = 0
        # getting the max x position
        last_chromosome = self.chromosomes[-1]
        xmax = self.cumulative_sizes[last_chromosome]

        fig = go.Figure()

        #######################################################################
        # HORIZONTAL LINES
        #######################################################################

        if suggestiveline_value:
            fig.add_hline(
                y=suggestiveline_value,
                line_width=suggestiveline_width,
                line_dash="dash",
                line_color=suggestiveline_color,
            )

        if genomewideline_value:
            fig.add_hline(
                y=genomewideline_value,
                line_width=genomewideline_width,
                line_dash="dash",
                line_color=genomewideline_color,
            )

        #######################################################################
        # VERTICAL LINES
        #######################################################################

        logger.info("Adding vertical lines")
        for chrom in tqdm(self.chromosomes[1:]):
            fig.add_vline(
                x=self.cumulative_previous_sizes[chrom], line_width=1, line_color="grey"
            )

        #######################################################################
        # DATA TO HIGHLIGHT
        #######################################################################

        tmp = pd.DataFrame()  # Empty DataFrame to contain the highlighted data

        if highlight:
            if not isinstance(highlight, bool):
                if self.snp_colname not in self.data.columns.values:
                    raise KeyError(
                        "snp argument specified for highlight as %s but "
                        "column not found in the data.frame" % self.snp_colname
                    )
            else:
                if not genomewideline_value:
                    raise Warning(
                        "The genomewideline_value you entered is not a "
                        "positive value, or False, you cannot set highlight "
                        "to True in that case."
                    )
                tmp = self.data

                # Sort the p-values (or -log10(p-values) above the line
                if genomewideline_value:
                    if self.logp:
                        tmp = tmp.loc[
                            -np.log10(tmp[self.pvalue_colname]) > genomewideline_value
                        ]
                    else:
                        tmp = tmp.loc[tmp[self.pvalue_colname] > genomewideline_value]

                highlight_hover_text = self.get_hover_text(tmp)

                x_values = tmp["cumulative_position"].values
                y_values = self.get_yvalues(tmp)

                if not tmp.empty:
                    fig.add_trace(
                        go.Scattergl(
                            x=x_values,
                            y=y_values,
                            mode="markers",
                            text=highlight_hover_text,
                            marker=dict(color=highlight_color, size=point_size),
                            name="Point(s) of interest",
                        )
                    )

        #######################################################################
        # PLOT VALUES
        #######################################################################

        # Remove the highlighted data from the DataFrame if not empty
        if tmp.empty:
            data = self.data
        else:
            data = self.data.loc[~self.data.index.isin(tmp.index)]

        if len(self.chromosomes) == 1:
            # If single chromosome, ticks and labels automatic.
            layout = go.Layout(
                title=title,
                xaxis={
                    "title": self.xlabel if xlabel is None else xlabel,
                    "showgrid": showgrid,
                    "range": [xmin, xmax],
                },
                yaxis={"title": ylabel},
                hovermode="closest",
            )

            hover_text = self.get_hover_text(data)
            y_values = self.get_yvalues(data)

            fig.add_trace(
                go.Scattergl(
                    x=data[self.position_colname].values,
                    y=y_values,
                    mode="markers",
                    showlegend=showlegend,
                    name=self.chromosomes[0],
                    marker={"size": point_size},
                    text=hover_text,
                )
            )
        else:
            # if multiple chrms, use the ticks and labels you created above.
            layout = go.Layout(
                title=title,
                xaxis={
                    "title": self.xlabel if xlabel is None else xlabel,
                    "showgrid": showgrid,
                    "range": [xmin, xmax],
                    "tickmode": "array",
                    "tickvals": self.ticks,
                    "ticktext": self.ticksLabels,
                    "ticks": "outside",
                    "autotickangles": [0, 30, 45, 60, 90],
                    # "tickangle": 45,
                },
                yaxis={"title": ylabel},
                hovermode="closest",
            )

            for i, chrom in tqdm(enumerate(self.chromosomes)):
                chrom_data = data[data[self.chrom_colname] == chrom]
                hover_text = self.get_hover_text(chrom_data)

                x_values = chrom_data["cumulative_position"].values
                y_values = self.get_yvalues(chrom_data)

                fig.add_trace(
                    go.Scattergl(
                        x=x_values,
                        y=y_values,
                        mode="markers",
                        showlegend=showlegend,
                        name=chrom,
                        marker={"size": point_size},
                        text=hover_text,
                    )
                )

        fig.update_layout(layout)
        return fig
