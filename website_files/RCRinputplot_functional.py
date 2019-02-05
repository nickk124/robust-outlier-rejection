from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button
from bokeh.plotting import Figure, output_file, show, save

import numpy as np

output_file("RCRinputplot_functional.html")

x = []
y = []

src = ColumnDataSource(data=dict(x=x, y=y))

p = Figure(plot_height = 400, plot_width = 600,
                    x_axis_label = 'x',
                    y_axis_label = 'y(x)')

p.circle('x', 'y', source=src, size=5, color="navy", alpha=0.5)

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p), code="""
    testCallBegin();

    testLog('Plotting Input');

    p.reset.emit();

    var input_data = getData();
    var x_data = input_data[0];
    var y_data = input_data[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    src.change.emit();

    testCallEnd();
""")
callback_linear_x = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    testLog('x -> linear');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformXLinear(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "x";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")
callback_linear_y = CustomJS(args=dict(src=src, p=p, axis=p.yaxis[0]), code="""
    testCallBegin();

    testLog('y -> linear');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformYLinear(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "y";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")
callback_exp_x = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    testLog('x -> exponential');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformXExp(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "10^x";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")
callback_exp_y = CustomJS(args=dict(src=src, p=p, axis=p.yaxis[0]), code="""
    testCallBegin();

    testLog('y -> exponential');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformYExp(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "10^y";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")
callback_log_x = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0]), code="""
    testCallBegin();

    testLog('x -> logarithmic');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformXLog(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "log(x) (base 10)";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")
callback_log_y = CustomJS(args=dict(src=src, p=p, axis=p.yaxis[0]), code="""
    testCallBegin();

    testLog('y -> logarithmic');

    p.reset.emit();

    var input_data = getData();
    var oldx_data = input_data[0];
    var oldy_data = input_data[1];

    var trans_result = transformYLog(oldx_data, oldy_data);

    var x_data = trans_result[0];
    var y_data = trans_result[1];

    var data = src.data;

    var x = [];
    var y = [];

    x = data['x'];
    y = data['y'];

    if (x_data.length != y_data.length){
        testLog('Alert: inputted data must be same size');
    }

    for (var i = 0; i < y_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
    }

    axis.axis_label = "log(y) (base 10)";

    p.change.emit()

    src.change.emit();

    testCallEnd();
""")


button_plot = Button(label = "Plot Data", button_type = "primary")
#button_test= Button(label = "Test", button_type = "success")

button_linear_x = Button(label = "Linear x Basis", button_type = "primary")
button_log_x = Button(label = "Logarithmic x Basis", button_type = "primary")
button_exp_x = Button(label = "Exponential x Basis", button_type = "primary")

button_linear_y = Button(label = "Linear y Basis", button_type = "primary")
button_log_y = Button(label = "Logarithmic y Basis", button_type = "primary")
button_exp_y = Button(label = "Exponential y Basis", button_type = "primary")

button_plot.js_on_click(callback_plot)

button_linear_x.js_on_click(callback_linear_x)
button_log_x.js_on_click(callback_log_x)
button_exp_x.js_on_click(callback_exp_x)

button_linear_y.js_on_click(callback_linear_y)
button_log_y.js_on_click(callback_log_y)
button_exp_y.js_on_click(callback_exp_y)

#button_test.js_on_click(callback_test)

buttons_1 = widgetbox(button_plot, sizing_mode = 'stretch_both')
buttons_2 = widgetbox(button_linear_x, button_log_x, button_exp_x, sizing_mode = 'stretch_both')
buttons_3 = widgetbox(button_linear_y, button_log_y, button_exp_y, sizing_mode = 'stretch_both')
buttons = column(buttons_1, buttons_2, buttons_3, sizing_mode = 'scale_width')
layout = row(p, buttons)

#show(layout)
save(layout)
