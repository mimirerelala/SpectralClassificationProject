import os
import numpy as np
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.embed import components
from starspectre.settings import BASE_DIR as root
colors = ['green', 'red', 'blue', 'yellow']
temp = os.path.join(root, 'website/temp')


def handle_uploaded_file(file, filename):
    with open(os.path.join(temp, filename), 'wb+') as dataout:
        for chunk in file.chunks():
            dataout.write(chunk)


def plot_handle(coordinates):
    PLOT_OPTIONS = dict(plot_width=800, plot_height=400)

    plot = figure(responsive=True, tools='reset,box_zoom,wheel_zoom,pan,resize,save', **PLOT_OPTIONS)
    for index in range(len(coordinates)):
        plot.line(coordinates[index][0], coordinates[index][1], line_color='{}'.format(colors[index]))
    
    resources = INLINE
    js_resources = resources.render_js()
    css_resources = resources.render_css()

    script, div = components({'plot': plot})

    return (js_resources, css_resources, script , div)



def file_parse(folder, file_name):
    for_plot = []
    coords = []
    for file in os.listdir(folder):
        if file.startswith(file_name):
            for_plot.append(os.path.join(folder, file_name))
    for item in for_plot:
        data = np.genfromtxt(item, dtype=None)
        data_x = data[:, 0]
        data_y = data[:, 1]
        coords.append([data_x, data_y])
    return coords
