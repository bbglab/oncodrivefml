import gzip
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import pandas as pd
import numpy as np
from math import pi
from bokeh import __version__ as bokeh_version
from bokeh.io import vform
from bokeh.plotting import output_notebook
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.models import HoverTool, PrintfTickFormatter, CustomJS, TextInput, Circle


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


TOOLS = "pan,box_zoom,resize,wheel_zoom,"
TOOLS += "reset,previewsave,crosshair,hover"


class QQPlot(object):

    def __init__(self):

        # Open the figure
        self._fig = figure(width=600, plot_height=600, tools=TOOLS)

        # Set the labels
        self._fig.xaxis.axis_label = 'Expected p-values'
        self._fig.yaxis.axis_label = 'Observed p-values'
        self._fig.axis.major_label_text_font_size = '14pt'
        self._fig.xaxis[0].formatter = PrintfTickFormatter(format="1e-%1.4f")
        self._fig.yaxis[0].formatter = PrintfTickFormatter(format="1e-%1.4f")
        self._fig.xaxis[0].major_label_orientation = pi/4
        hover = self._fig.select(dict(type=HoverTool))
        hover.mode = 'mouse'
        hover.tooltips = None
        self._widgets = None
        self._source = None
        self._cr = None
        self._df = None
        self._dict_for_df = None
        self._extra_fields = None
        self._data = None
        self._max_x = None
        self._max_y = None

    def add_tooltips(self, tooltips):

        # Hover
        code = """if (cb_data.index['1d'].indices.length > 0)\n"""
        code += """{\n"""
        code += """     var s = source.get('data');\n"""
        code += """     if ( $( ".bk-tooltip.bk-tooltip-custom.bk-left" ).length == 0 )\n"""
        code += """     {\n"""
        code += """         $( ".bk-canvas-overlays" ).append( '<div class="bk-tooltip bk-tooltip-custom bk-left" style="z-index: 1010; top: 0px; left: 0px; display: block;"></div>' );\n"""
        code += """     }\n"""
        code += """     var inner = "<div><div><div>";\n"""
        code += """     for(i in cb_data.index['1d'].indices)\n"""
        code += """     {\n"""
        code += """         index = cb_data.index['1d'].indices[i];\n"""
        code += """         if (i > 2) break;\n"""
        code += """         inner = inner + """+tooltips+"""\n"""
        code += """     };\n"""
        code += """     inner = inner + "<div><span style='font-size: 15px; color: #009;'>TOTAL OF: " + cb_data.index['1d'].indices.length + "</span></div>";\n"""
        code += """     inner = inner + "</div></div></div>";\n"""
        code += """     $('.bk-tooltip.bk-tooltip-custom.bk-left')[0].innerHTML = inner;\n"""
        code += """     $('.bk-tooltip.bk-tooltip-custom.bk-left').attr('style', 'left:' + (cb_data.geometry.sx+10) + 'px; top:' + (cb_data.geometry.sy-5-$('.bk-tooltip.bk-tooltip-custom.bk-left').height()/2) + 'px; z-index: 1010; display: block;');\n"""
        code += """}else\n"""
        code += """{\n"""
        code += """     $( "div" ).remove( ".bk-tooltip.bk-tooltip-custom.bk-left" );\n"""
        code += """}\n"""

        callback = CustomJS(args={'source': self._source}, code=code)
        self._fig.add_tools(HoverTool(tooltips=None, callback=callback, renderers=[self._cr], mode='mouse'))

    def add_search_fields(self, search_fields, position=0):
        for key, value in search_fields.items():
            code_text_box = """
                origSearch = cb_obj.get('value');
                search = origSearch.toUpperCase();
                var selected = source.get('selected')['1d'].indices;
                searcher = source.get('data')."""+value+"""
                selected.length = 0
                for(index in searcher)
                {
                    if ( searcher[index].toUpperCase().indexOf(search) > -1)
                        selected.push(index);
                }
                if (selected.length == 0)
                    /*alert("Value not found: '"+origSearch+"'")*/
                    swal("Value not found: '"+origSearch+"'")
                source.trigger('change');"""
            #self._widgets.insert( position, TextInput(value="", title=key, name=key, callback=CustomJS(args=dict(source=self._source),  code=code_text_box)))

    def load(self, input_file, fields=None, basic_fields=None, extra_fields=None):

        self._df = pd.read_csv(input_file, header=0, sep="\t")
        inv_fields = {v: k for k, v in basic_fields.items()}
        self._df = self._df.rename(columns = inv_fields)

        self._dict_for_df = { 'pvalue': self._df['pvalue'].tolist(), 'fdr': self._df['qvalue']}
        self._extra_fields = extra_fields

        for key, value in self._extra_fields.items():
            self._dict_for_df[key] = self._df[value].tolist()

    def add_basic_plot(self):

        # Default settings
        min_samples = 2
        min_pvalue = 10000

        colors = ['royalblue', 'blue']

        self._dict_for_df ['observed'] = self._df['pvalue'].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / min_pvalue))
        self._dict_for_df ['color'] = self._df['num_samples'].map(lambda x: colors[1] if x >= min_samples else colors[0])
        self._dict_for_df ['alpha'] = self._df['num_samples'].map(lambda x: 0.7 if x >= min_samples else 0.3)

        self._data = pd.DataFrame( self._dict_for_df)

        self._data.sort(columns=['observed'], inplace=True)
        exp_pvalues = -np.log10(np.arange(1, len(self._data) + 1) / float(len(self._data)))
        exp_pvalues.sort()
        self._data['expected'] = exp_pvalues

        dict_for_source=dict(
            x=self._data['expected'].tolist(),
            y=self._data['observed'].tolist(),
            color=self._data['color'].tolist(),
            alpha=self._data['alpha'].tolist(),
            pvalue=[str(x) for x in self._data["pvalue"]],
            qvalue=[str(x) for x in self._data["fdr"]]
        )
        for key, value in self._extra_fields.items():
            dict_for_source[key] = self._data[key].tolist()

        self._source = ColumnDataSource( data = dict_for_source )

        # Plot the first set of data
        if len(self._data['expected']) > 0:
            invisible_circle = Circle(x='x', y='y', fill_color='color', fill_alpha='alpha', line_color=None, size=10)
            visible_circle = Circle(x='x', y='y', fill_color='color', fill_alpha=0.9, line_color='red', size=10)
            self._cr = self._fig.add_glyph(self._source, invisible_circle, selection_glyph=visible_circle, nonselection_glyph=invisible_circle)

        # Get the maximum pvalues (observed and expected)
        self._max_x = float(self._data[['expected']].apply(np.max))
        self._max_y = float(self._data[['observed']].apply(np.max))

        # Give some extra space (+-5%)
        self._max_x *= 1.1
        self._max_y *= 1.1

        # Add a dashed diagonal from (min_x, max_x) to (min_y, max_y)
        self._fig.line(np.linspace(0, min(self._max_x, self._max_y)),
                np.linspace(0, min(self._max_x, self._max_y)),
                color='red', line_width=2, line_dash=[5, 5])

        # Set the grid
        self._fig.grid.grid_line_alpha = 0.8
        self._fig.grid.grid_line_dash = [6, 4]
        self._widgets = [self._fig]

    def add_cutoff(self):
        colors = ['royalblue', 'blue']

        # FDR
        genes_to_annotate = []
        for fdr_cutoff, fdr_color in zip((0.25, 0.1), ('green', 'red')):
            fdr = self._data[self._data['fdr'] < fdr_cutoff]['observed']
            if len(fdr) > 0:
                fdr_y = np.min(fdr)
                fdr_x = np.min(self._data[self._data['observed'] == fdr_y]['expected'])
                self._fig.line((fdr_x - self._max_x * 0.025, fdr_x + self._max_x * 0.025), (fdr_y, fdr_y),
                        color=fdr_color, line_width=2)

                # Add the name of the significant genes
                genes = self._data[(self._data['observed'] >= fdr_y) & (self._data['expected'] >= fdr_x)]
                for count, line in genes.iterrows():
                    genes_to_annotate.append({'x': line['expected'],
                                              'y': line['observed'],
                                              'HugoID': line['HugoID'],
                                              'color': colors[0] if line['color'] == colors[0] else fdr_color[0]})

    def show(self, output_path, showit=True, notebook=False):
        """Create an interactive qqplot"""

        # Import modules
        if notebook:
            output_notebook()
        layout = vform( *self._widgets)

        # Save the plot
        if output_path is not None:
            output_file(output_path)

        if showit:
            show(layout)
        else:
            script, div = components(layout)
            html = """  <!DOCTYPE html>
                        <html>
                        <script src="https://code.jquery.com/jquery-1.11.3.min.js">
                        </script><script src="https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.1.3/sweetalert.min.js"></script>
                        <link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/sweetalert/1.1.3/sweetalert.min.css">
                        <link href="http://cdn.pydata.org/bokeh/release/bokeh-"""+bokeh_version+""".min.css" rel="stylesheet" type="text/css">
                        <script src="http://cdn.pydata.org/bokeh/release/bokeh-"""+bokeh_version+""".min.js"></script>\n"""+\
                        script+\
                        """<body>\n"""+\
                        div+"""
                        </body>
                        </html>"""
            text_file = open(output_path, "w")
            text_file.write(html)
            text_file.close()


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

    qqp = QQPlot()
    qqp.load(input_file = input_file, basic_fields = {'num_samples': 'samples_mut', 'pvalue': 'pvalue', 'qvalue': 'qvalue' },
            extra_fields = {'HugoID': 'symbol', 'EnsblID': 'index'})
    qqp.add_basic_plot()
    qqp.add_cutoff()
    qqp.add_search_fields( {'Hugo ID': 'HugoID', 'Ensembl ID': 'EnsblID'}, position = 0)
    qqp.add_tooltips(""" "<div>\\
                             <span style='font-size: 17px; font-weight: bold;'>\" + s.HugoID[index] + \"</span> \\
                             <span style='font-size: 15px; color: #966;'>[\" + s.EnsblID[index] + \"]</span> \\
                          </div> \\
                          <div> \\
                             <span style='font-size: 15px;'>p/q-value</span> \\
                             <span style='font-size: 10px; color: #696;'>(\" + s.pvalue[index] + \", \" + s.qvalue[index] + \")</span> \\
                          </div> \\
                          </br>" """)
    qqp.show(output_path = output_path, notebook = False, showit=showit)