from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button
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

src = ColumnDataSource(data=dict(arr_hist = arr_hist, left = left, right = right))

p = Figure(plot_height = 400, plot_width = 600,
                    title = 'Histogram of Inputted Data',
                    x_axis_label = 'y (data)',
                    y_axis_label = 'number of datapoints')

p.quad(source = src, bottom = 0, top = 'arr_hist', left = 'left', right = 'right', fill_color = 'cornflowerblue', line_color = 'black')

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p), code="""
    testCallBegin();

    testLog('Initial Plot');

    p.reset.emit();

    var y_data = getData()[0];

    makeDefaultBinCount(y_data);

    var hist_result = getHistBins(y_data);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    src.change.emit();

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

callback_addbins = CustomJS(args=dict(src=src, p=p), code="""
    testCallBegin();

    testLog('Plus a bin');

    p.reset.emit();

    var y_data = getData()[0];

    if (basis === 'exp'){
        y_data = transformExp(y_data);
    } else if (basis === 'log'){
        y_data = transformLog(y_data);
    }


    bincount += 1;

    var hist_result = getHistBins(y_data);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    src.change.emit();

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_subtractbins = CustomJS(args=dict(src=src, p=p), code="""
    testCallBegin();

    testLog('Minus a bin');

    p.reset.emit();

    var y_data = getData()[0];

    if (basis === 'exp'){
        y_data = transformExp(y_data);
    } else if (basis === 'log'){
        y_data = transformLog(y_data);
    }

    bincount -= 1;

    var hist_result = getHistBins(y_data);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    arr_hist.pop();
    left.pop();
    right.pop();


    src.change.emit();

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_linear = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    basis = 'linear';

    testLog('Changing to linear');

    p.reset.emit();

    var y_data = getData()[0];

    var hist_result = getHistBins(y_data);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    axis.axis_label = "y (data)";

    p.change.emit()

    src.change.emit();

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_log = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    basis = 'log'

    testLog('Changing to logarithmic');

    p.reset.emit();

    var y_data = getData()[0];

    var y_trans = transformLog(y_data);

    var hist_result = getHistBins(y_trans);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    axis.axis_label = "log(y) (base 10)";

    p.change.emit();

    src.change.emit();

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")

callback_exp = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    basis = 'exp';

    testLog('Changing to exponential');

    p.reset.emit();

    var y_data = getData()[0];

    var y_trans = transformExp(y_data);

    var hist_result = getHistBins(y_trans);
    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1];
    var right_result = hist_result[2];

    var data = src.data;

    var arr_hist = data['arr_hist'];
    var left = data['left'];
    var right = data['right'];

    for (var i = 0; i < arr_hist_result.length; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    axis.axis_label = "10^y";

    p.change.emit();

    src.change.emit();

    //testLog([arr_hist.length, left.length, right.length]);

    testCallEnd();
""")


button_plot = Button(label = "Plot Histogram", button_type = "primary")
#button_test= Button(label = "Test", button_type = "success")

button_addbins = Button(label = "Add Bins", button_type = "primary")
button_subtractbins = Button(label = "Subtract Bins", button_type = "primary")

button_linear = Button(label = "Linear Basis", button_type = "primary")
button_log = Button(label = "Logarithmic Basis", button_type = "primary")
button_exp = Button(label = "Exponential Basis", button_type = "primary")

button_plot.js_on_click(callback_plot)

button_subtractbins.js_on_click(callback_subtractbins)
button_addbins.js_on_click(callback_addbins)

button_linear.js_on_click(callback_linear)
button_log.js_on_click(callback_log)
button_exp.js_on_click(callback_exp)

#button_test.js_on_click(callback_test)

buttons_1 = widgetbox(button_plot, button_addbins, button_subtractbins, sizing_mode = 'stretch_both')
buttons_2 = widgetbox(button_linear, button_log, button_exp, sizing_mode = 'stretch_both')
buttons = column(buttons_1, buttons_2, sizing_mode = 'scale_width')
layout = row(p, buttons)

#show(layout)
save(layout)
