import gzip
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import logging
from collections import OrderedDict
from bokeh.plotting import figure, output_file, show, gridplot, ColumnDataSource, save
from bokeh.models import HoverTool, PrintfTickFormatter


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


TOOLS = "pan,box_zoom,resize,wheel_zoom,"
TOOLS += "reset,previewsave,crosshair,hover"


def eliminate_duplicates(df):

    colors = df['color'].tolist()
    x = df['x'].tolist()[0]
    y = df['y'].tolist()[0]
    if 'r' in colors:
        return x, y, 'red'
    elif 'g' in colors:
        return x, y, 'green'
    else:
        return x, y, 'grey'


def add_symbol(df):
    ensemble_file = os.path.join(__location__, "ensemble_genes_75.txt.gz")
    gene_conversion = {line.split("\t")[0]: line.strip().split("\t")[-1]
                       for line in gzip.open(ensemble_file, 'rt').readlines()}
    df.loc[:, 'symbol'] = df[df.columns[0]].apply(lambda e: gene_conversion.get(e, e))
    return df


def qqplot_png(input_file, output_file, showit=False):
    pvalue='pvalue'
    qvalue='qvalue'
    min_samples=2
    draw_greys=True
    annotate=True
    cut=True
    #############################

    MIN_PVALUE = 10000
    MAX_ANNOTATED_GENES = 50

    df = pd.read_csv(input_file, header=0, sep="\t")
    df.dropna(subset=[pvalue], inplace=True)

    # Define the shape of the figure
    NROW = 1
    NCOL = 1

    fig = plt.figure(figsize=(6, 6))
    axs = [plt.subplot2grid((NROW, NCOL), (item // NCOL, item % NCOL)) for item in range(NROW * NCOL)]

    # Plot is on the right or on the left?
    dx_side = True
    ax = axs[0]

    colors = ['royalblue', 'blue']
    obs_pvalues = df[pvalue].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / MIN_PVALUE))
    obs_color = df['samples_mut'].map(lambda x: colors[1] if x >= min_samples else colors[0])
    obs_alpha = df['samples_mut'].map(lambda x: 0.7 if x >= min_samples else 0.3)

    data = pd.DataFrame({'HugoID': df['symbol'],
                         'observed': obs_pvalues,
                         'color': obs_color,
                         'alpha': obs_alpha,
                         'fdr': df[qvalue] if qvalue is not None else 1
                         })

    data.sort(columns=['observed'], inplace=True)
    exp_pvalues = -np.log10(np.arange(1, len(data) + 1) / float(len(data)))
    exp_pvalues.sort()
    data['expected'] = exp_pvalues

    # Get the maximum pvalues (observed and expected)
    max_x = float(data[['expected']].apply(np.max))
    max_y = float(data[['observed']].apply(np.max))

    # Give some extra space (+-5%)
    min_x = max_x * -0.05
    min_y = max_y * -0.05
    max_x *= 1.1
    max_y *= 1.1

    grey = data[data['color'] == colors[0]]
    blue = data[data['color'] == colors[1]]

    # Plot the data
    if draw_greys and len(grey['expected']) > 0:
        ax.scatter(grey['expected'].tolist(),
                   grey['observed'].tolist(),
                   color=grey['color'].tolist(),
                   alpha=grey['alpha'].tolist()[0],
                   s=30)

    if len(blue['expected']) > 0:
        ax.scatter(blue['expected'].tolist(),
                   blue['observed'].tolist(),
                   color=blue['color'].tolist(),
                   alpha=blue['alpha'].tolist()[0],
                   s=30)

    # Get the data that are significant with a FDR < 0.1 and FDR 0.25
    genes_to_annotate = []
    for fdr_cutoff, fdr_color in zip((0.25, 0.1), ('g-', 'r-')):
        fdr = data[data['fdr'] < fdr_cutoff]['observed']
        if len(fdr) > 0:
            fdr_y = np.min(fdr)
            fdr_x = np.min(data[data['observed'] == fdr_y]['expected'])
            ax.plot((fdr_x - max_x * 0.025, fdr_x + max_x * 0.025), (fdr_y, fdr_y), fdr_color)
            # Add the name of the significant genes
            genes = data[(data['observed'] >= fdr_y) & (data['expected'] >= fdr_x)]
            for count, line in genes.iterrows():
                if line['color'] == colors[0]:
                    continue
                genes_to_annotate.append({'x': line['expected'],
                                          'y': line['observed'],
                                          'HugoID': line['HugoID'],
                                          'color': colors[0] if line['color'] == colors[0] else fdr_color[0]})

    # Annotate the genes
    genes_annotated = 0
    if annotate and len(genes_to_annotate) > 0:
        genes_to_annotate = pd.DataFrame(genes_to_annotate)

        # Get rid of redundant genes
        grouped = genes_to_annotate.groupby(['HugoID'])
        grouped = grouped.apply(eliminate_duplicates)
        grouped = pd.DataFrame({'HugoID': grouped.index.tolist(),
                                'x': [x[0] for x in grouped],
                                'y': [y[1] for y in grouped],
                                'color': [c[2] for c in grouped]})
        grouped.sort(columns=['y', 'x'], inplace=True, ascending=[False, False])

        x_text = max_x * 1.1 if dx_side is True else min_x - (max_x * 0.2)
        y_text = np.floor(max_y)
        distance_between_genes = max_y * 0.05
        for count, line in grouped.iterrows():
            x, y = line['x'], line['y']
            ax.annotate(line['HugoID'], xy=(x, y), xytext=(x_text, y_text),
                        arrowprops=dict(facecolor="black", shrink=0.03,
                                        width=1, headwidth=6, alpha=0.3),
                        horizontalalignment="left", verticalalignment="center",
                        color=line['color'],
                        weight = 'normal'
                        )
            y_text -= distance_between_genes

            # This avoid the crash for ValueError: width and height must each be below 32768
            genes_annotated += 1
            if genes_annotated >= MAX_ANNOTATED_GENES:
                logging.warning("Annotations cut to {} genes".format(MAX_ANNOTATED_GENES))
                break

    # Add labels
    ax.set_xlabel("expected pvalues")
    ax.set_ylabel("observed pvalues")

    # Add the dashed diagonal
    ax.plot(np.linspace(0, np.floor(max(max_x, max_y))),
            np.linspace(0, np.floor(max(max_x, max_y))),
            'r--')
    ax.grid(True)

    # Redefine the limits of the plot
    if cut:
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)

    # Set the title: project, cancer_type
    #ax.set_title(title)

    # Adjust the plot
    try:
        plt.tight_layout()
    except ValueError as e:
        logging.warning('Ignoring tight_layout()')

    # Save the plot
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')

    # Show the plot
    if showit:
        plt.show()

    # Close the figure
    plt.close()


def qqplot_html(input_file, output_path, showit=False):

    pvalue = 'pvalue'
    qvalue = 'qvalue'
    min_samples = 2
    draw_greys = True
    MIN_PVALUE = 10000

    df = pd.read_csv(input_file, header=0, sep="\t")
    df.dropna(subset=[pvalue], inplace=True)

    colors = ['royalblue', 'blue']
    obs_pvalues = df[pvalue].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / MIN_PVALUE))
    obs_color = df['samples_mut'].map(lambda x: colors[1] if x >= min_samples else colors[0])
    obs_alpha = df['samples_mut'].map(lambda x: 0.7 if x >= min_samples else 0.3)

    data = pd.DataFrame({'HugoID': df['symbol'],
                         'GeneID': df['index'],
                         'observed': obs_pvalues,
                         'color': obs_color,
                         'alpha': obs_alpha,
                         'pvalue': df[pvalue].tolist(),
                         'fdr': df[qvalue] if qvalue is not None else 1
                         })

    data.sort(columns=['observed'], inplace=True)
    exp_pvalues = -np.log10(np.arange(1, len(data) + 1) / float(len(data)))
    exp_pvalues.sort()
    data['expected'] = exp_pvalues

    # Open the figure
    ax = figure(width=600, plot_height=600, tools=TOOLS)

    # Set the labels
    ax.xaxis.axis_label = 'Expected p-values'
    ax.yaxis.axis_label = 'Observed p-values'
    ax.axis.major_label_text_font_size = '14pt'
    ax.xaxis[0].formatter = PrintfTickFormatter(format="1e-%1.0f")
    ax.yaxis[0].formatter = PrintfTickFormatter(format="1e-%1.0f")

    # Create two set of data
    grey = data[data['color'] == colors[0]]
    blue = data[data['color'] == colors[1]]

    source_grey = ColumnDataSource(
        data=dict(
            x=grey['expected'].tolist(),
            y=grey['observed'].tolist(),
            symbol=grey['HugoID'].tolist(),
            geneid=grey['GeneID'].tolist(),
            pvalue=[str(x) for x in grey["pvalue"]],
            qvalue=[str(x) for x in grey["fdr"]]
        )
    )

    source_blue = ColumnDataSource(
        data=dict(
            x=blue['expected'].tolist(),
            y=blue['observed'].tolist(),
            symbol=blue['HugoID'].tolist(),
            geneid=blue['GeneID'].tolist(),
            pvalue=[str(x) for x in blue["pvalue"]],
            qvalue=[str(x) for x in blue["fdr"]]
        )
    )

    # Plot the first set of data
    if draw_greys and len(grey['expected']) > 0:
        ax.circle(x=grey['expected'],
                  y=grey['observed'],
                  color=grey['color'],
                  source=source_grey,
                  fill_alpha=grey['alpha'].tolist()[0],
                  size=10,
                  name='greys')

    if len(blue['expected']) > 0:
        ax.circle(x=blue['expected'],
                  y=blue['observed'],
                  color=blue['color'],
                  source=source_blue,
                  fill_alpha=blue['alpha'].tolist()[0],
                  size=10,
                  name='blues')

    # Get the maximum pvalues (observed and expected)
    max_x = float(data[['expected']].apply(np.max))
    max_y = float(data[['observed']].apply(np.max))

    # Give some extra space (+-5%)
    max_x *= 1.1
    max_y *= 1.1

    # Add a dashed diagonal from (min_x, max_x) to (min_y, max_y)
    ax.line(np.linspace(0, min(max_x, max_y)),
            np.linspace(0, min(max_x, max_y)),
            color='red', line_width=2, line_dash=[5, 5])

    # Set the grid
    #ax.xgrid.grid_line_color = None
    ax.grid.grid_line_alpha = 0.8
    ax.grid.grid_line_dash = [6, 4]

    # FDR
    genes_to_annotate = []
    for fdr_cutoff, fdr_color in zip((0.25, 0.1), ('green', 'red')):
        fdr = blue[blue['fdr'] < fdr_cutoff]['observed']
        if len(fdr) > 0:
            fdr_y = np.min(fdr)
            fdr_x = np.min(blue[blue['observed'] == fdr_y]['expected'])
            ax.line((fdr_x - max_x * 0.025, fdr_x + max_x * 0.025), (fdr_y, fdr_y),
                    color=fdr_color, line_width=2)

            # Add the name of the significant genes
            genes = blue[(blue['observed'] >= fdr_y) & (blue['expected'] >= fdr_x)]
            for count, line in genes.iterrows():
                genes_to_annotate.append({'x': line['expected'],
                                          'y': line['observed'],
                                          'HugoID': line['HugoID'],
                                          'color': colors[0] if line['color'] == colors[0] else fdr_color[0]})

    # Hover
    hover = ax.select(dict(type=HoverTool))
    hover.point_policy = "snap_to_data"
    hover.names = ['greys', 'blues']
    hover.tooltips = OrderedDict([
        ("symbol", "@symbol"),
        ("geneid", "@geneid"),
        ("p-value", "@pvalue"),
        ("q-value", "@qvalue"),
    ])

    # Put the subplots in a gridplot
    fig = gridplot([[ax]], toolbar_location='above')

    # Save the plot
    if output_path is not None:
        output_file(output_path)

    if showit:
        show(fig)
    else:
        save(fig)
