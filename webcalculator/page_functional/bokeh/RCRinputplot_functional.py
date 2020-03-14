from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button, Div
from bokeh.plotting import Figure, output_file, show, save
from bokeh.events import MouseEnter

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
callback_plot = CustomJS(args=dict(src=src, p=p, xaxis=p.xaxis[0], yaxis=p.yaxis[0]), code="""
    p.reset.emit();

    var input_data = getData(); //getData also transforms
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

    x = x.slice(0,y_data.length);
    y = y.slice(0,y_data.length);

    data['x'] = x;
    data['y'] = y;

    //xaxis.axis_label = xaxislabel;
    //yaxis.axis_label = yaxislabel;
    xaxislabel = 'x'
    yaxislabel = 'y'

    p.change.emit();
    src.change.emit();
""")
callback_linear_x = CustomJS(code="""
    x_basis = "linear"
    xaxislabel = 'x'
""")
callback_linear_y = CustomJS(code="""
    y_basis = "linear"
    yaxislabel = 'y'
""")
callback_exp_x = CustomJS(code="""
    x_basis = "exp"
    xaxislabel = '10^x'
""")
callback_exp_y = CustomJS(code="""
    y_basis = "exp"
    yaxislabel = '10^y'
""")
callback_log_x = CustomJS(code="""
    x_basis = "log"
    xaxislabel = 'log10(x)'
""")
callback_log_y = CustomJS(code="""
    y_basis = "log"
    yaxislabel = 'log10(y)'
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

button_linear_x.js_on_click(callback_plot)
button_log_x.js_on_click(callback_plot)
button_exp_x.js_on_click(callback_plot)

button_linear_y.js_on_click(callback_plot)
button_log_y.js_on_click(callback_plot)
button_exp_y.js_on_click(callback_plot)

buttons_1 = widgetbox(button_plot)#, sizing_mode = 'stretch_both')
buttons_2 = row(button_linear_x, button_log_x, button_exp_x)#, sizing_mode = 'stretch_both')
buttons_3 = row(button_linear_y, button_log_y, button_exp_y)#, sizing_mode = 'stretch_both')
# buttons = column(buttons_1, buttons_2, buttons_3, sizing_mode = 'scale_width')

basisLabel = Div(text="""
<br>
<h1 style="font-size:2.5em;">Select the Basis of the Dependent Variable</h1>
""",
width=800, height=70)

layout = column(row(p, buttons_1), column(basisLabel, buttons_2, buttons_3))

# show(layout)
save(layout)