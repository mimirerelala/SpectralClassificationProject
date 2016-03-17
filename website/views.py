import os
import numpy as np
from django.shortcuts import render
from django.http import HttpResponse
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.embed import components
from bokeh.util.browser import view
from . sshkey import ssh_key
from starspectre.settings import BASE_DIR as root
from website import forms as f
from website import helpers as h


# from example_input_2 import data_x as data_1_x
# from example_input_2 import data_y as data_1_y
# from example_input_3 import data_x as data_2_x
# from example_input_3 import data_y as data_2_y
data_1_x, data_1_y = np.loadtxt(os.path.join(root,'example_input2.txt'), delimiter='  ', usecols=(0, 1), unpack=True)
data_2_x, data_2_y = np.loadtxt(os.path.join(root,'example_input3.txt'), delimiter='  ', usecols=(0, 1), unpack=True)
# Create your views here.


def index(request):
    return HttpResponse("Hi there")


def smart_chart(request):
    # file upload
    if request.method == 'POST':
        form = f.UploadFileForm(request.POST, request.FILES)
        filename = request.POST.get('session_key') or 'testfilename'
        if form.is_valid():
            h.handle_uploaded_file(request.FILES['file'], filename)
            feedback = 'Success'
        else:
            feedback = 'Error'
    else:
        form = f.UploadFileForm()


    # Plot creation
    PLOT_OPTIONS = dict(plot_width=800, plot_height=400)
    
    green = figure(responsive=True, tools='reset,box_zoom,wheel_zoom,pan,resize,save', **PLOT_OPTIONS)
    green.line(data_1_x, data_1_y, line_color='green')
    green.line(data_2_x, data_2_y, line_color='red')

    resources = INLINE

    js_resources = resources.render_js()
    css_resources = resources.render_css()

    script, div = components({'green': green})

    return render(request, 'smart_chart.html', locals())


def ssh_key_show(request):
    return HttpResponse(ssh_key)


