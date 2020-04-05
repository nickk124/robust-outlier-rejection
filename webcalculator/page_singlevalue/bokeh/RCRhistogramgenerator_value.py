#----- Author: Nick Konz -----#

from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button, BoxZoomTool, Div
from bokeh.plotting import Figure, output_file, show, save
from bokeh.models.tools import HoverTool
from bokeh.events import LODEnd
from copy import copy

import numpy as np

output_file("RCRhistogram_value.html")

arr_hist = []
left = []
right = []
#xMin = 0
#xMax = 1

src = ColumnDataSource(data=dict(arr_hist = arr_hist, left = left, right = right))

TOOLS="xpan,xwheel_zoom, xzoom_in, xzoom_out, reset"

p = Figure(plot_height = 400, plot_width = 600, x_range = (0, 1),# x_range must be manually initialized to avoid Bokeh's auto-ranging
                    x_axis_label = 'measured value',
                    y_axis_label = 'weighted number of measurements', 
                    tools = TOOLS, 
                    active_scroll='xwheel_zoom', 
                    active_drag = "xpan")

p.add_tools(HoverTool(tooltips=[("count", "@arr_hist")]))

p.quad(source = src, bottom = 0, top = 'arr_hist', left = 'left', right = 'right', fill_color = 'cornflowerblue', line_color = 'black')

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range, y_range=p.y_range), code="""
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
        if (!currentlyChangingBasis){
            x_range.start = left[0];
            x_range.end = right[bincount-1];
        }
        // xMin = x_range.start;
        // xMax = x_range.end;

        //console.log("Changed range based off of hist, stored:", xMin, xMax);
    } else {
        p.x_range.start = xMin 
        p.x_range.end = xMax;
        x_range.start = xMin 
        x_range.end = xMax;

        //console.log("Fixed range to previous stored:", xMin, xMax)
    }

    y_range.start = 0;

    src.change.emit();
    p.change.emit();
""")

callback_plot2 = copy(callback_plot)

callback_store_range = CustomJS(args=dict(x_range=p.x_range), code="""
    xMin = x_range.start;
    xMax = x_range.end;
    console.log("storing range as:", xMin, xMax);
""")

callback_maintain_range = CustomJS(args=dict(x_range=p.x_range), code="""
    let xMinNew = xMin;
    let xMaxNew = xMax;
    x_range.start = xMinNew;
    x_range.end = xMaxNew;
    x_range.change.emit();
""")

callback_maintain_range_basis = CustomJS(args=dict(x_range=p.x_range), code="""
    console.log("prev. range:", xMin, xMax);
    let xMinNew = xMin;
    let xMaxNew = xMax;
    if (xMinNew == 0){
        xMinNew = Number.MIN_VALUE;
    }
    if (xMaxNew == 0){
        xMaxNew = Number.MIN_VALUE;
    }

    console.log(prevbasis);
    console.log(basis);
    switch (prevbasis){
        case "linear":
            switch (basis){
                case "linear":
                    break;
                case "log":
                    xMinNew = Math.log10(xMinNew);
                    xMaxNew = Math.log10(xMaxNew);
                    break;
                case "exp":
                    xMinNew = Math.pow(10, xMinNew);
                    xMaxNew = Math.pow(10, xMaxNew);
                    break;
                default:
                    break;
            }
            break;
        case "exp":
            switch (basis){
                case "linear":
                    xMinNew = Math.log10(xMinNew);
                    xMaxNew = Math.log10(xMaxNew);
                    break;
                case "log":
                    xMinNew = Math.log10(Math.log10(xMinNew));
                    xMaxNew = Math.log10(Math.log10(xMaxNew));
                    break;
                case "exp":
                    break;
                default:
                    break;
            }
            break;
        case "log":
            switch (basis){
                case "linear":
                    xMinNew = Math.pow(10, xMinNew);
                    xMaxNew = Math.pow(10, xMaxNew);
                    break;
                case "log":
                    break;
                case "exp":
                    xMinNew = Math.pow(10, Math.pow(10, xMinNew));
                    xMaxNew = Math.pow(10, Math.pow(10, xMaxNew));
                    break;
                default:
                    break;
            }
            break;
        default:
            break;
	}

    if (isNaN(xMinNew)){
        xMinNew = 0;
    }

    if (isNaN(xMaxNew)){
        xMaxNew = 1;
    }

    if (!isFinite(xMinNew) || !isFinite(xMaxNew)){
        alert("Error: Datapoint too large for exponential basis!");
        xMinNew = 0;
        xMaxNew = 1;
    }

    x_range.start = xMinNew;
    x_range.end = xMaxNew;
    console.log("new range:", xMinNew, xMaxNew);
    x_range.change.emit();
    currentlyChangingBasis = false;
""")

# callback_addbins = CustomJS(args=dict(x_range=p.x_range), code="""
#     bincount += 1;
#     canChangeBins = true;
#     canChangeRange = false;
# """)

# callback_subtractbins = CustomJS(code="""
#     bincount -= 1;
#     canChangeBins = true;
#     canChangeRange = false;
# """)

callback_changetolinear = CustomJS(code="""
    currentlyChangingBasis = true;
    prevbasis = basis
    basis = "linear"
    xaxislabel = 'measured value'
    canChangeBins = true;
""")

callback_changetolog= CustomJS(code="""
    currentlyChangingBasis = true;
    prevbasis = basis
    basis = "log"
    xaxislabel = 'log10(measured value)'
    canChangeBins = true;
""")

callback_changetoexp = CustomJS(code="""
    currentlyChangingBasis = true;
    prevbasis = basis
    basis = "exp"
    xaxislabel = '10^(measured value)'
    canChangeBins = true;
""")

callback_resetrange = CustomJS(code="""
    canChangeRange = true;
""")

callback_add_visible_bins = CustomJS(args=dict(x_range=p.x_range), code="""
    let bincount_prev = bincount;
    let xrange_displayed = x_range.end - x_range.start;
    let xrange_data = xMaxData - xMinData;
    let bincount_displayed = Math.round(bincount_default * xrange_displayed / xrange_data)
    console.log("current displayed bincount:", bincount_displayed); // this is correct, I checked
    
    let binwidth = xrange_displayed / bincount_displayed; // also correct
    console.log("current (displayed) binwidth:", binwidth); 

    console.log("default bincount:", bincount_default);
    bincount = Math.round((bincount_default + 1) * xrange_data / xrange_displayed);
    canChangeBins = true;
    canChangeRange = false;
    console.log("Autoadjusting bincount from:", bincount_prev, "to:", bincount);
""")

callback_subtract_visible_bins = CustomJS(args=dict(x_range=p.x_range), code="""
    let bincount_prev = bincount;
    let xrange_displayed = x_range.end - x_range.start;
    let xrange_data = xMaxData - xMinData;
    // let bincount_displayed = Math.round(bincount_default * xrange_displayed / xrange_data)
    // console.log("current displayed bincount:", bincount_displayed); // this is correct, I checked
    
    // let binwidth = xrange_displayed / bincount_displayed; // also correct
    // console.log("current (displayed) binwidth:", binwidth); 

    // console.log("default bincount:", bincount_default);
    bincount = Math.round((bincount_default - 1) * xrange_data / xrange_displayed);
    canChangeBins = true;
    canChangeRange = false;
    console.log("Autoadjusting bincount from:", bincount_prev, "to:", bincount);
""")
#interactivity

button_plot = Button(label = "Plot Histogram", button_type = "primary")
button_addbins = Button(label = "Add Bins", button_type = "primary")
button_subtractbins = Button(label = "Subtract Bins", button_type = "primary")

button_linear = Button(label = "Linear Basis (Default)", button_type = "primary")
button_log = Button(label = "Logarithmic Basis", button_type = "primary")
button_exp = Button(label = "Exponential Basis", button_type = "primary")

button_plot.js_on_click(callback_resetrange)
button_plot.js_on_click(callback_plot)
button_plot.js_on_click(callback_plot2)

# button_subtractbins.js_on_click(callback_store_range)
# button_addbins.js_on_click(callback_store_range)

# button_subtractbins.js_on_click(callback_subtractbins)
# button_addbins.js_on_click(callback_addbins)

# button_subtractbins.js_on_click(callback_plot)
# button_addbins.js_on_click(callback_plot)

# button_subtractbins.js_on_click(callback_maintain_range)
# button_addbins.js_on_click(callback_maintain_range)

button_linear.js_on_click(callback_store_range)
button_log.js_on_click(callback_store_range)
button_exp.js_on_click(callback_store_range)

button_linear.js_on_click(callback_changetolinear)
button_log.js_on_click(callback_changetolog)
button_exp.js_on_click(callback_changetoexp)

button_linear.js_on_click(callback_plot)
button_log.js_on_click(callback_plot)
button_exp.js_on_click(callback_plot)

button_linear.js_on_click(callback_maintain_range_basis)
button_log.js_on_click(callback_maintain_range_basis)
button_exp.js_on_click(callback_maintain_range_basis)

# autoadjust bins upon adding/subtracting

button_subtractbins.js_on_click(callback_store_range)
button_subtractbins.js_on_click(callback_subtract_visible_bins)
button_subtractbins.js_on_click(callback_plot)
button_subtractbins.js_on_click(callback_maintain_range)

button_addbins.js_on_click(callback_store_range)
button_addbins.js_on_click(callback_add_visible_bins)
button_addbins.js_on_click(callback_plot)
button_addbins.js_on_click(callback_maintain_range)

# layout
buttons_0 = widgetbox(button_plot)
buttons_1 = widgetbox(button_addbins, button_subtractbins)#, sizing_mode = 'stretch_both')
buttons = column(buttons_0, buttons_1, sizing_mode = 'scale_width')
buttons_2 = row(button_linear, button_log, button_exp, sizing_mode = 'fixed')#, sizing_mode = 'stretch_both')

basisLabel = Div(text="""
<br>
<h2 style="position:absolute; left:-25px; font-size:1.66666666667em;">Select Basis:</h2>
""",
width=200, height=70)

layoutplot = column(row(p, buttons), column(basisLabel, buttons_2))

#show(layoutplot)
save(layoutplot)