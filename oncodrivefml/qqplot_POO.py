# Import modules
import pandas as pd
import numpy as np
from bokeh.plotting import figure, show, gridplot, ColumnDataSource, save
from bokeh.models import HoverTool, PrintfTickFormatter, CustomJS, TextInput, Circle
from math import pi
from bokeh.io import vform
from bokeh.plotting import output_notebook, output_file
from bokeh.embed import components

# Global variables
TOOLS = "pan,box_zoom,resize,wheel_zoom,reset,previewsave,crosshair,hover,tap"


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



    def add_tooltips( self, tooltips):

        # Hover
        #code =  """debugger;\n"""
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

    def add_search_fields(self, searchFields, position = 0):
        for key, value in searchFields.items():
            code_text_box = """
                search = cb_obj.get('value').toUpperCase();
                var selected = source.get('selected')['1d'].indices;
                searcher = source.get('data')."""+value+"""
                selected.length = 0
                for(index in searcher)
                {
                    if ( searcher[index].toUpperCase().indexOf(search) > -1)
                        selected.push(index);
                }
                source.trigger('change');"""
            self._widgets.insert( position, TextInput(value="", title=key, name=key, callback=CustomJS(args=dict(source=self._source),  code=code_text_box)))


    def load(self, input_file, fields=None, basic_fields = None, extra_fields = None):

        self._df = pd.read_csv(input_file, header=0, sep="\t")
        inv_fields = {v: k for k, v in basic_fields.items()}
        self._df = self._df.rename(columns = inv_fields)

        self._dict_for_df = { 'pvalue': self._df['pvalue'].tolist(),
                        'fdr': self._df['qvalue']
                        }
        self._extra_fields = extra_fields

        for key, value in self._extra_fields.items():
            self._dict_for_df[key] = self._df[value].tolist()


    def add_basic_plot(self):

        # Default settings
        min_samples = 2
        MIN_PVALUE = 10000

        colors = ['royalblue', 'blue']

        self._dict_for_df ['observed'] = self._df['pvalue'].map(lambda x: -np.log10(x) if x > 0 else -np.log10(1 / MIN_PVALUE))
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
        #ax.xgrid.grid_line_color = None
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
        '''Create an interactive qqplot'''
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
                        <script src="https://code.jquery.com/jquery-1.11.3.min.js"></script>
                        <link href="http://cdn.pydata.org/bokeh/release/bokeh-0.10.0.min.css" rel="stylesheet" type="text/css">
                        <script src="http://cdn.pydata.org/bokeh/release/bokeh-0.10.0.min.js"></script>\n"""+\
                        script+\
                        """<body>\n"""+\
                        div+"""
                        </body>
                        </html>"""
            text_file = open( output_path, "w")
            text_file.write(html)
            text_file.close()

            #save(layout)
