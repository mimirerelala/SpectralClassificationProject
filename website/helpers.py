import os
import numpy as np
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.embed import components
from starspectre.settings import BASE_DIR as root
colors = ['green', 'red', 'orange', 'blue', 'black', 'yellow', 'brown', 'grey']
temp = os.path.join(root, 'website/temp')


def handle_uploaded_file(foldername, file, filename):
    os.chdir(temp)
    try:
        os.mkdir(foldername)
    except:
        pass
    session_folder = os.path.join(temp, foldername)
    with open(os.path.join(session_folder, filename), 'wb+') as dataout:
        for chunk in file.chunks():
            dataout.write(chunk)


def file_parse(folder):
    for_plot = []
    for file in os.listdir(folder):
            for_plot.append(os.path.join(folder, file))
    return for_plot


def add_coords(filename):
    '''
    data_x, data_y = np.loadtxt(filename, usecols=(0,1), unpack=True)
    '''
    coords = []
    data = np.genfromtxt(filename, dtype=None)
    data_x = data[:, 0]
    data_y = data[:, 1]
    coords.append(data_x)
    coords.append(data_y)
    return coords


def gen_coordinates(folder):
    coordinates = []
    for i in file_parse(folder):
        a = add_coords(i)
        coordinates.append(a)
    return [coordinates, file_parse(folder)]


def plot_handle(data_list):
    coordinates = data_list[0]
    files = [os.path.basename(i) for i in data_list[1]]
    PLOT_OPTIONS = dict(plot_width=800, plot_height=400, toolbar_location="left")
    plot = figure(responsive=True, tools='reset,box_zoom,wheel_zoom,pan,resize,save,hover', **PLOT_OPTIONS)
    for index in range(len(coordinates)):
        plot.line(coordinates[index][0], coordinates[index][1], line_color='{}'.format(colors[index]), legend='{}'.format(files[index]))

    resources = INLINE
    js_resources = resources.render_js()
    css_resources = resources.render_css()
    script, div = components({'plot': plot})

    return (js_resources, css_resources, script , div)