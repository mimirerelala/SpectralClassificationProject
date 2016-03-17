import random
from django.shortcuts import render
from django.http import HttpResponse
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.embed import components
from bokeh.util.browser import view
from . sshkey import ssh_key
from example_input_2 import data_x as data_1_x
from example_input_2 import data_y as data_1_y
from example_input_3 import data_x as data_2_x
from example_input_3 import data_y as data_2_y
# Create your views here.


def index(request):
    return HttpResponse("Hi there")


def smart_chart(request):
    PLOT_OPTIONS = dict(plot_width=800, plot_height=400)
    # SCATTER_OPTIONS = dict(size=12, alpha=0.5)

    # red = figure(responsive=True, tools='pan', **PLOT_OPTIONS)
    # red.scatter(data(), data(), color="red", **SCATTER_OPTIONS)
    # blue = figure(responsive=False, tools='pan', **PLOT_OPTIONS)
    # blue.scatter(data(), data(), color="blue", **SCATTER_OPTIONS)
    green = figure(responsive=True, tools='reset,box_zoom,wheel_zoom,pan,resize,save', **PLOT_OPTIONS)
    green.line(data_1_x, data_1_y, line_color='green')
    green.line(data_2_x, data_2_y, line_color='red')

    resources = INLINE

    js_resources = resources.render_js()
    css_resources = resources.render_css()

    script, div = components({'green': green})

    return render(request, 'smart_chart.html', {'js_resources':js_resources, 'css_resources':css_resources, 'plot_script':script, 'plot_div':div})


def ssh_key_show(request):
    return HttpResponse(ssh_key)


