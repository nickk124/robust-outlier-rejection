#----- Author: Nick Konz -----#

from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button, BoxZoomTool, Div
from bokeh.plotting import Figure, output_file, show, save
from bokeh.models.tools import HoverTool

import numpy as np

output_file("RCRhistogram_value.html")

arr_hist = []
left = []
right = []
xMin = 0
xMax = 1

src = ColumnDataSource(data=dict(arr_hist = arr_hist, left = left, right = right))

TOOLS="xpan,xwheel_zoom, xzoom_in, xzoom_out, reset"

p = Figure(plot_height = 400, plot_width = 600, x_range = (xMin, xMax),
                    x_axis_label = 'measured value',
                    y_axis_label = 'weighted number of measurements', 
                    tools = TOOLS, 
                    active_scroll='xwheel_zoom', 
                    active_drag = "xpan")

p.add_tools(HoverTool(tooltips=[("count", "@arr_hist")]))

p.quad(source = src, bottom = 0, top = 'arr_hist', left = 'left', right = 'right', fill_color = 'cornflowerblue', line_color = 'black')

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""
    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData(); //getData also transforms
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    var hist_result = getHistBins(y_data, w_data);
    axis.axis_label = xaxislabel;

    makeDefaultBinCount(y_data); //if this callback was activated by a bin changer button, this doesn't reset bincount

    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = [];
    var left = [];
    var right = [];

    arr_hist = data['arr_hist'];
    left = data['left'];
    right = data['right'];

    for (var i = 0; i < bincount; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    if (canChangeRange){
        x_range.start = left[0];
        x_range.end = right[bincount-1];

        xMin = x_range.start;
        xMax = x_range.end;
    } else {
        p.x_range.start = xMin 
        p.x_range.end = xMax;
    }

    p.change.emit();
    src.change.emit();
""")

callback_addbins = CustomJS(code="""
    bincount += 1;
    canChangeBins = true;
    canChangeRange = false;
""")

callback_subtractbins = CustomJS(code="""
    bincount -= 1;
    canChangeBins = true;
    canChangeRange = false;
""")

callback_changetolinear = CustomJS(code="""
    basis = "linear"
    xaxislabel = 'measured value'
    canChangeBins = true;
    canChangeRange = false;
""")

callback_changetolog= CustomJS(code="""
    basis = "log"
    xaxislabel = 'log10(measured value)'
    canChangeBins = true;
""")

callback_changetoexp = CustomJS(code="""
    basis = "exp"
    xaxislabel = '10^(measured value)'
    canChangeBins = true;
""")

callback_resetrange = CustomJS(code="""
    canChangeRange = true;
""")

#interactivity

button_plot = Button(label = "Plot Histogram", button_type = "primary")
button_addbins = Button(label = "Add Bins", button_type = "primary")
button_subtractbins = Button(label = "Subtract Bins", button_type = "primary")

button_linear = Button(label = "Linear Basis (Default)", button_type = "primary")
button_log = Button(label = "Logarithmic Basis", button_type = "primary")
button_exp = Button(label = "Exponential Basis", button_type = "primary")

button_plot.js_on_click(callback_plot)
button_plot.js_on_click(callback_resetrange)

button_subtractbins.js_on_click(callback_subtractbins)
button_addbins.js_on_click(callback_addbins)
button_subtractbins.js_on_click(callback_plot)
button_addbins.js_on_click(callback_plot)

button_linear.js_on_click(callback_changetolinear)
button_log.js_on_click(callback_changetolog)
button_exp.js_on_click(callback_changetoexp)
button_linear.js_on_click(callback_plot)
button_log.js_on_click(callback_plot)
button_exp.js_on_click(callback_plot)

#p.js_on_event(Pan, callbackPlot)

# layout
buttons_0 = widgetbox(button_plot)
buttons_1 = widgetbox(button_addbins, button_subtractbins)#, sizing_mode = 'stretch_both')
buttons = column(buttons_0, buttons_1, sizing_mode = 'scale_width')
buttons_2 = row(button_linear, button_log, button_exp, sizing_mode = 'fixed')#, sizing_mode = 'stretch_both')

basisLabel = Div(text="""
<br>
<h2 style="position:absolute; left:-25px;>Select Basis</h2>
""",
width=200, height=70)

layoutplot = column(row(p, buttons), column(basisLabel, buttons_2))

# show(layoutplot)
save(layoutplot)