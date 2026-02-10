"""
This original script was directly copied from https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/_manhattan.py
The _get_hover_text function was copied https://github.com/plotly/dash-bio/blob/master/dash_bio/component_factory/utils.py

As of February 2026, Dash Bio has not been updated since 3 years, and there are bugs.
This code is a modified version of the original code.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go

SUGGESTIVE_LINE_LABEL = "suggestive line"
GENOMEWIDE_LINE_LABEL = "genomewide line"


def _get_hover_text(df: pd.DataFrame):
    return (
        "SNP: "
        + df["snp"].astype(str)
        + "<br>GENE: "
        + df["gene"].astype(str)
        + "<br>P-VALUE (QUANTILE 5%): "
        + df["pvalue"].astype(str)
        + "<br>MEAN QUALITY: "
        + df["quality"].astype(str)
        + "<br>MEAN TOTAL DEPTH: "
        + df["total_depth"].astype(str)
        + "<br><br>ALLELE COUNTS (REF / ALT):<br>"
        + df["allele_counts"].astype(str)
    )


def ManhattanPlot(
    dataframe: pd.DataFrame,
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
    highlight_color="red",
):
    """Returns a figure for a manhattan plot.

    Keyword arguments:
    - dataframe (dataframe; required): A pandas dataframe which must contain at
        least the following three columns:
                - the chromosome number
                - genomic base-pair corrected_position
                - a numeric quantity to plot such as a p-value or zscore
    - logp (bool; optional): If True, the -log10 of the p-value is
        plotted. It isn't very useful to plot raw p-values; however,
        plotting the raw value could be useful for other genome-wide plots
        (e.g., peak heights, Bayes factors, test statistics, other
        "scores", etc.)
    - title (string; default 'Manhattan Plot'): The title of the graph.
    - showgrid (bool; default true): Boolean indicating whether gridlines
        should be shown.
    - xlabel (string; optional): Label of the x axis.
    - ylabel (string; default '-log10(p)'): Label of the y axis.
    - point_size (number; default 5): Size of the points of the Scatter
        plot.
    - showlegend (bool; default true): Boolean indicating whether legends
        should be shown.
    - col (string; optional): A string representing the color of the
        points of the scatter plot. Can be in any color format accepted by
        plotly.graph_objects.
    - suggestiveline_value (bool | float; default 8): A value which must
        be either False to deactivate the option, or a numerical value
        corresponding to the p-value at which the line should be drawn.
        The line has no influence on the data points.
    - suggestiveline_color (string; default 'grey'): Color of the suggestive
      line.
    - suggestiveline_width (number; default 2): Width of the suggestive
        line.
    - genomewideline_value (bool | float; default -log10(5e-8)): A boolean
        which must be either False to deactivate the option, or a numerical value
        corresponding to the p-value above which the data points are
        considered significant.
    - genomewideline_color (string; default 'red'): Color of the genome-wide
        line. Can be in any color format accepted by plotly.graph_objects.
    - genomewideline_width (number; default 1): Width of the genome-wide
      line.
    - highlight (bool; default True): turning on/off the highlighting of
        data points considered significant.
    - highlight_color (string; default 'red'): Color of the data points
        highlighted because they are significant. Can be in any color
        format accepted by plotly.graph_objects.

        # ...
        Example 1: Random Manhattan Plot
        '''
        dataframe = pd.DataFrame(
            np.random.randint(0,100,size=(100, 3)),
            columns=['P', 'CHR', '"corrected_position"'])
        fig = create_manhattan(dataframe, title='XYZ Manhattan plot')

        plotly.offline.plot(fig, image='png')
        '''

    """

    mh = _ManhattanPlot(dataframe, logp=logp)

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
        df: pd.DataFrame,
        logp=True,
    ):
        """
        Keyword arguments:
        - dataframe (dataframe; required): A pandas dataframe which
        must contain at least the following three columns:
            - the chromosome number
            - genomic base-pair corrected_position
            - a numeric quantity to plot such as a p-value or zscore
        - logp (bool; default True): If True, the -log10 of the p-value is
        plotted.  It isn't very useful to plot raw p-values; however,
        plotting the raw value could be useful for other genome-wide plots
        (e.g., peak heights, Bayes factors, test statistics, other
        "scores", etc.).

        Returns:
        - A ManhattanPlot object."""

        self.data = df
        self.xlabel = ""
        self.ticks = []
        self.ticksLabels = []
        self.nChr = len(df["chromosome"].unique())
        self.logp = logp

        # This section sets up positions and ticks. Ticks should be placed in
        # the middle of a chromosome. The new pos column is added that keeps
        # a running sum of the corrected_positions of each successive chromosome.
        # For example:
        # chrm "corrected_position" pos
        # 1   1  1
        # 1   2  2
        # 2   1  3
        # 2   2  4
        # 3   1  5

        if self.nChr == 1:
            # For a single chromosome
            self.data["corrected_position"] = self.data["position"]
            self.ticks.append(int(len(self.data["corrected_position"]) / 2.0) + 1)
            self.xlabel = f"Position in {self.data['chromosome'].unique()[0]}"
            self.ticksLabels = self.ticks
        else:
            # For multiple chromosomes
            lastbase = 0
            previous_chrom = None
            for i, chrom in enumerate(self.data["chromosome"].unique()):
                if i == 0:
                    self.data.loc[
                        self.data["chromosome"] == chrom, "corrected_position"
                    ] = self.data.loc[
                        self.data["chromosome"] == chrom, "position"
                    ].values

                    previous_chrom = chrom

                else:
                    print("previous_chrom", previous_chrom)
                    print("data", self.data)
                    previous_corrected_position = self.data.loc[
                        self.data["chromosome"] == previous_chrom, "position"
                    ]
                    # Shift the basepair corrected_position by the largest "corrected_position" of the
                    # current chromosome
                    print("previous_corrected_position", previous_corrected_position)
                    print("lastbase", lastbase)
                    lastbase = lastbase + previous_corrected_position.iat[-1]

                    self.data.loc[
                        self.data["chromosome"] == chrom, "corrected_position"
                    ] = (
                        self.data.loc[
                            self.data["chromosome"] == chrom, "position"
                        ].values
                        + lastbase
                    )

                tmin = min(
                    self.data.loc[
                        self.data["chromosome"] == chrom, "corrected_position"
                    ]
                )
                tmax = max(
                    self.data.loc[
                        self.data["chromosome"] == chrom, "corrected_position"
                    ]
                )
                self.ticks.append(int((tmin + tmax) / 2.0) + 1)

                previous_chrom = chrom

            self.xlabel = "Chromosome"
            self.data["corrected_position"] = self.data["corrected_position"].astype(
                self.data["corrected_position"].dtype
            )

            if self.nChr > 10:  # To avoid crowded labels
                self.ticksLabels = [
                    t
                    if np.mod(int(t), 2)  # Only every two ticks
                    else ""
                    for t in self.data["chromosome"].unique()
                ]
            else:
                self.ticksLabels = self.data["chromosome"].unique()  # All the ticks

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
        highlight_color="red",
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
        - highlight_color (string; default 'red'): Color of the data
            points highlighted because they are significant. Can be in any
            color format accepted by plotly.graph_objects.

        Returns:
        - A figure formatted for plotly.graph_objects.

        """

        xmin = min(self.data["corrected_position"].values)
        xmax = max(self.data["corrected_position"].values)

        horizontallines = []

        if suggestiveline_value:
            suggestiveline = go.layout.Shape(
                name=SUGGESTIVE_LINE_LABEL,
                type="line",
                fillcolor=suggestiveline_color,
                line=dict(color=suggestiveline_color, width=suggestiveline_width),
                x0=xmin,
                x1=xmax,
                xref="x",
                y0=suggestiveline_value,
                y1=suggestiveline_value,
                yref="y",
            )
            horizontallines.append(suggestiveline)

        if genomewideline_value:
            genomewideline = go.layout.Shape(
                name=GENOMEWIDE_LINE_LABEL,
                type="line",
                fillcolor=genomewideline_color,
                line=dict(color=genomewideline_color, width=genomewideline_width),
                x0=xmin,
                x1=xmax,
                xref="x",
                y0=genomewideline_value,
                y1=genomewideline_value,
                yref="y",
            )
            horizontallines.append(genomewideline)

        data_to_plot = []  # To contain the data traces
        tmp = pd.DataFrame()  # Empty DataFrame to contain the highlighted data

        if highlight:
            if not isinstance(highlight, bool):
                if "snp" not in self.data.columns.values:
                    raise KeyError(
                        "snp argument specified for highlight as %s but "
                        "column not found in the data.frame" % "snp"
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
                        tmp = tmp.loc[-np.log10(tmp["pvalue"]) > genomewideline_value]
                    else:
                        tmp = tmp.loc[tmp["pvalue"] > genomewideline_value]

                highlight_hover_text = _get_hover_text(tmp)

                y_values = (
                    -np.log10(tmp["pvalue"].values)
                    if self.logp
                    else tmp["pvalue"].values
                )

                if not tmp.empty:
                    data_to_plot.append(
                        go.Scattergl(
                            x=tmp["corrected_position"].values,
                            y=y_values,
                            mode="markers",
                            text=highlight_hover_text,
                            marker=dict(color=highlight_color, size=point_size),
                            name="Point(s) of interest",
                        )
                    )

        # Remove the highlighted data from the DataFrame if not empty
        if tmp.empty:
            data = self.data
        else:
            data = self.data.drop(self.data.index[tmp.index])

        if self.nChr == 1:
            if col is None:
                col = ["black"]

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

            hover_text = _get_hover_text(data)

            chromosome_name = str(data["chromosome"].unique()[0])

            y_values = (
                -np.log10(data["pvalue"].values) if self.logp else data["pvalue"].values
            )

            data_to_plot.append(
                go.Scattergl(
                    x=data["corrected_position"].values,
                    y=y_values,
                    mode="markers",
                    showlegend=showlegend,
                    name=chromosome_name,
                    marker={"color": col[0], "size": point_size},
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
                },
                yaxis={"title": ylabel},
                hovermode="closest",
            )

            icol = 0
            if col is None:
                col = ["black" if np.mod(i, 2) else "grey" for i in range(self.nChr)]

            for chrom in data["chromosome"].unique():
                tmp = data[data["chromosome"] == chrom]

                chromosome_name = str(data["chromosome"].unique()[0])

                hover_text = _get_hover_text(tmp)

                y_values = (
                    np.log10(tmp["pvalue"].values)
                    if self.logp
                    else tmp["pvalue"].values
                )

                data_to_plot.append(
                    go.Scattergl(
                        x=tmp["corrected_position"].values,
                        y=-y_values,
                        mode="markers",
                        showlegend=showlegend,
                        name=chromosome_name,
                        marker={"color": col[icol], "size": point_size},
                        text=hover_text,
                    )
                )

                icol = icol + 1

        layout.shapes = horizontallines

        return go.Figure(data=data_to_plot, layout=layout)
