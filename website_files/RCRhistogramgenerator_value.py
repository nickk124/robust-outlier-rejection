from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button, BoxZoomTool
from bokeh.plotting import Figure, output_file, show, save

import numpy as np

output_file("RCRhistogram_value.html")

"""
y = [1,1,1,2,3,4,9,10,11,15,16,17,18,55]

N = len(y)

bincount = np.int(np.sqrt(float(N)))
range = [np.min(y), np.max(y)]

arr_hist, edges = np.histogram(y, bins = bincount, range = range ) #bins is the number of bins
left = edges[:-1]
right = edges[1:]
print(edges)
print(bincount)
print(range[1] - range[0])
"""

arr_hist = []
left = []
right = []
xMin = 0
xMax = 1

src = ColumnDataSource(data=dict(arr_hist = arr_hist, left = left, right = right))

p = Figure(plot_height = 400, plot_width = 600, x_range = (xMin, xMax),
                    x_axis_label = 'measured value',
                    y_axis_label = 'weighted number of measurements', tools = "xpan, xwheel_zoom", active_scroll='xwheel_zoom', active_drag = "xpan")

p.quad(source = src, bottom = 0, top = 'arr_hist', left = 'left', right = 'right', fill_color = 'cornflowerblue', line_color = 'black')

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""
    testCallBegin();

    testLog('Initial Plot');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    makeDefaultBinCount(y_data);

    var hist_result;

    if (basis === "exp"){
        var y_trans = transformExp(y_data);
        y_data = y_trans;
        hist_result = getHistBins(y_trans, w_data);
        axis.axis_label = "10^(measured value)";

    } else if (basis === "log"){
        var y_trans = transformLog(y_data);
        y_data = y_trans;
        hist_result = getHistBins(y_trans, w_data);
        axis.axis_label = "log(measured value) (base 10)";

    } else if (basis === "linear"){
        hist_result = getHistBins(y_data, w_data);
        axis.axis_label = "measured value";
    }

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

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    x_range.start = left[0];
    x_range.end = right[bincount-1];

    testLog(right[bincount-1]);

    p.change.emit();
    src.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_test = CustomJS(args=dict(src=src), code="""
    testCallBegin();

    testLog('Hello');

    var y_data = [1.0,1.1,1.2,1.3,5.0,6.0,7.0,9.0,44.2,100.1];

    testLog(y_data);

    var inpData = getData()[0];

    testLog(inpData);

    var hist_result = getHistBins(y_data);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    testLog(arr_hist_result);
    testLog(left_result);
    testLog(right_result);

    testCallEnd();
""")

callback_addbins = CustomJS(args=dict(src=src, p=p, x_range=p.x_range), code="""
    testCallBegin();

    var oldxMin = xMin;
    var oldxMax = xMax;

    testLog(oldxMin);
    testLog(oldxMax);

    //x_range.setv({"start": oldxMin, "end": oldxMax});

    testLog('Plus a bin');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    if (basis === 'exp'){
        y_data = transformExp(y_data);
    } else if (basis === 'log'){
        y_data = transformLog(y_data);
    }


    bincount += 1;

    var hist_result = getHistBins(y_data, w_data);
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

    //x_range.setv({"start": oldxMin, "end": oldxMax});

    testLog(oldxMin);
    testLog(oldxMax);

    testLog(xMin);
    testLog(xMax);

    x_range.start = oldxMin;
    x_range.end = oldxMax;

    src.change.emit();

    //x_range.setv({"start": oldxMin, "end": oldxMax});
    p.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_subtractbins = CustomJS(args=dict(src=src, p=p, x_range=p.x_range), code="""
    testCallBegin();

    var oldxMin = xMin;
    var oldxMax = xMax;

    testLog('Minus a bin');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    if (basis === 'exp'){
        y_data = transformExp(y_data);
    } else if (basis === 'log'){
        y_data = transformLog(y_data);
    }

    bincount -= 1;

    var hist_result = getHistBins(y_data, w_data);
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

    arr_hist.pop();
    left.pop();
    right.pop();

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    x_range.start = oldxMin;
    x_range.end = oldxMax;

    src.change.emit();
    p.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_linear = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""
    testCallBegin();

    var oldxMin = xMin;
    var oldxMax = xMax;

    basis = 'linear';

    testLog('Changing to linear');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    var hist_result = getHistBins(y_data, w_data);
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

    axis.axis_label = "measured value";

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    x_range.start = oldxMin;
    x_range.end = oldxMax;

    src.change.emit();
    p.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_log = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""
    testCallBegin();

    var oldxMin = xMin;
    var oldxMax = xMax;

    basis = 'log'

    testLog('Changing to logarithmic');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    var y_trans = transformLog(y_data);

    var hist_result = getHistBins(y_trans, w_data);
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

    axis.axis_label = "log(measured value) (base 10)";

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    if (oldxMin > 0){
        x_range.start = Math.log10(oldxMin);
    } else {
        x_range.start = 0
    }
    x_range.end = Math.log10(oldxMax);

    src.change.emit();
    p.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_exp = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""
    testCallBegin();

    var oldxMin = xMin;
    var oldxMax = xMax;

    basis = 'exp';

    testLog('Changing to exponential');

    p.reset.emit();

    hasWeights = $('#Weighted').data('clicked');
	hasErrorBars = $('#ErrorBars').data('clicked');

    var data_res = getData();
    var y_data = data_res[0];
    var w_data = [];

    if (hasWeights || hasErrorBars){
        w_data = data_res[1];
    }

    var y_trans = transformExp(y_data);

    var hist_result = getHistBins(y_trans, w_data);
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

    axis.axis_label = "10^(measured value)";

    arr_hist = arr_hist.slice(0,bincount);
    left = left.slice(0,bincount);
    right = right.slice(0,bincount);

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    x_range.start = Math.pow(10,oldxMin);
    x_range.end = Math.pow(10,oldxMax);

    src.change.emit();
    p.change.emit();

    testLog(arr_hist);
    testLog(left);
    testLog(right);

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

# make it so that if the user has zoomed in or panned around, that when they change basis or change bin count, the plot boundaries are unchanged
#the global values for these boundaries will be in linear basis
p.x_range.callback = CustomJS(args=dict(src=src), code="""
    var start = cb_obj.start;
    var end = cb_obj.end;

    switch (basis){
        case 'exp':
            xMin = Math.log10(start);
            xMax = Math.log10(end);
            break;
        case 'log':
            xMin = Math.pow(10, start);
            xMax = Math.pow(10, end);
            break;
        default:
            xMin = start;
            xMax = end;
            break;
    }
""")


button_plot = Button(label = "Plot Histogram", button_type = "primary")
#button_test= Button(label = "Test", button_type = "success")

button_addbins = Button(label = "Add Bins", button_type = "primary")
button_subtractbins = Button(label = "Subtract Bins", button_type = "primary")

button_linear = Button(label = "Linear Basis (Default)", button_type = "primary")
button_log = Button(label = "Logarithmic Basis", button_type = "primary")
button_exp = Button(label = "Exponential Basis", button_type = "primary")

button_plot.js_on_click(callback_plot)

button_subtractbins.js_on_click(callback_subtractbins)
button_addbins.js_on_click(callback_addbins)

button_linear.js_on_click(callback_linear)
button_log.js_on_click(callback_log)
button_exp.js_on_click(callback_exp)

#button_test.js_on_click(callback_test)
buttons_0 = widgetbox(button_plot)
buttons_1 = widgetbox(button_addbins, button_subtractbins, sizing_mode = 'stretch_both')
buttons_2 = widgetbox(button_linear, button_log, button_exp, sizing_mode = 'stretch_both')
buttons = column(buttons_0, buttons_1, buttons_2, sizing_mode = 'scale_width')
layout = row(p, buttons)

#show(layout)
save(layout)
