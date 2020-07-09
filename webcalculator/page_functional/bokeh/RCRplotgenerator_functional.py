from bokeh.layouts import column
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button
from bokeh.plotting import Figure, output_file, show, save

output_file("RCRplot_functional.html")

x_rejected = [] #initial data (list)
y_rejected = []
x_nonrejected = [] #initial data (list)
y_nonrejected = []
x_original = []
y_fitted = []

sourcenon = ColumnDataSource(data=dict(x_nonrejected=x_nonrejected, y_nonrejected=y_nonrejected)) #ColumnDataSource is the object where the data of the graph is stored.
sourcerej = ColumnDataSource(data=dict(x_rejected=x_rejected, y_rejected=y_rejected)) #ColumnDataSource is the object where the data of the graph is stored.
sourceall = ColumnDataSource(data=dict(x_original=x_original, y_fitted=y_fitted)) #ColumnDataSource is the object where the data of the graph is stored.

plot = Figure(plot_width=600, plot_height=600)
#plotting code
plot.circle('x_nonrejected', 'y_nonrejected', source=sourcenon, size=5, color="navy", alpha=0.5, legend="nonrejected data")# nonrejected
plot.circle('x_rejected', 'y_rejected', source=sourcerej, size=5, color="red", alpha=0.5, legend="rejected data")# rejected
plot.line('x_original', 'y_fitted', source=sourceall, color="green", line_width=2, legend="fitted model")

#plot.circle('x', 'y', source=source, line_width=3, line_alpha=0.6) #plotting by name

#defines the callback to be used:
callback = CustomJS(args=dict(sourcenon=sourcenon, sourcerej=sourcerej, sourceall=sourceall, p=plot), code="""
    console.log('Plot Generating...');

    p.reset.emit();

    var datanon = sourcenon.data;
    var datarej = sourcerej.data;
    var dataall = sourceall.data;

    var x_rejected = datarej['x_rejected']
    var y_rejected = datarej['y_rejected']
    var x_nonrejected = datanon['x_nonrejected']
    var y_nonrejected = datanon['y_nonrejected']
    var x_original = dataall['x_original']
    var y_fitted = dataall['y_fitted']

    for (var i = 0; i < x_rejected_result.length; i++) {
        x_rejected[i] = x_rejected_result[i]
        y_rejected[i] = y_rejected_result[i]
    }
    for (var i = 0; i < x_nonrejected_result.length; i++) {
        x_nonrejected[i] = x_nonrejected_result[i]
        y_nonrejected[i] = y_nonrejected_result[i]
    }

    x_original = [];
    y_fitted = [];

    // console.log(x_original_result.length);

    for (var i = 0; i < x_original_result.length; i++) {
        x_original.push(x_original_result[i]);
        y_fitted.push(y_fitted_result[i]);
    }

    datarej['x_rejected'] = x_rejected;
    datarej['y_rejected'] = y_rejected;
    datanon['x_nonrejected'] = x_nonrejected;
    datanon['y_nonrejected'] = y_nonrejected;
    dataall['x_original'] = x_original;
    dataall['y_fitted'] = y_fitted;

    // console.log(dataall['x_original']);
    // console.log(dataall['y_fitted']);

    sourcenon.change.emit();
    sourcerej.change.emit();
    sourceall.change.emit();
    p.change.emit();
""")

#slider = Slider(start=0.1, end=4, value=1, step=.1, title="power")
#slider.js_on_change('value', callback)

plot.xaxis.axis_label = 'x'
plot.yaxis.axis_label = 'y'

button = Button(label = "Plot Results", button_type = "primary")
button.js_on_click(callback)
layout = column(button, plot)

#show(layout)
save(layout)
