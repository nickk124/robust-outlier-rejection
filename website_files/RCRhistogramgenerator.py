from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button
from bokeh.plotting import Figure, output_file, show, save

import numpy as np

output_file("RCRhistogram.html")

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
callback_plot = CustomJS(args=dict(src=src), code="""
    testCallBegin();

    var y_data = getData();

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

    testCallEnd();
""")

callback_test = CustomJS(args=dict(src=src), code="""
    testCallBegin();

    testLog('Hello');

    var y_data = [1.0,1.1,1.2,1.3,5.0,6.0,7.0,9.0,44.2,100.1];

    testLog(y_data);

    var inpData = getData();

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

callback_addbins = CustomJS(args=dict(src=src), code="""
    testCallBegin();

    var y_data = getData();

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

    testCallEnd();
""")

callback_subtractbins = CustomJS(args=dict(src=src), code="""
    testCallBegin();

    var y_data = getData();

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

    src.change.emit();

    testCallEnd();
""")

button_plot = Button(label = "Plot Histogram", button_type = "primary")
#button_test= Button(label = "Test", button_type = "success")

button_addbins = Button(label = "Add Bins", button_type = "primary")
button_subtractbins = Button(label = "Subtract Bins", button_type = "primary")
button_linear = Button(label = "Linear", button_type = "primary")
button_log = Button(label = "Logarithmic", button_type = "primary")
button_exp = Button(label = "Exponential", button_type = "primary")

button_plot.js_on_click(callback_plot)
button_subtractbins.js_on_click(callback_subtractbins)
button_addbins.js_on_click(callback_addbins)

#button_test.js_on_click(callback_test)

buttons_1 = row(button_plot, button_addbins, button_subtractbins)
buttons_2 = row(button_linear, button_log, button_exp)
buttons = column(buttons_1, buttons_2, sizing_mode = 'scale_width')
layout = row(p, buttons)

show(layout)
save(layout)
