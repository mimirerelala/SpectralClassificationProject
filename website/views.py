import os
from random import choice
from string import ascii_letters
import numpy as np
from django.shortcuts import render
from django.http import HttpResponse
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.embed import components
from . sshkey import ssh_key
from starspectre.settings import BASE_DIR as root
from website import forms as f
from website import helpers as h
filename = 0

# from example_input_2 import data_x as data_1_x
# from example_input_2 import data_y as data_1_y
# from example_input_3 import data_x as data_2_x
# from example_input_3 import data_y as data_2_y
data_1_x, data_1_y = np.loadtxt(os.path.join(root,'example_input2.txt'), delimiter='  ', usecols=(0, 1), unpack=True)
data_2_x, data_2_y = np.loadtxt(os.path.join(root,'example_input3.txt'), delimiter='  ', usecols=(0, 1), unpack=True)

data1 = np.genfromtxt(os.path.join(root,'example_input2.txt'))
data_1_x = data1[:,0]
data_1_y = data1[:,1]
# Create your views here.


def index(request):
    return HttpResponse("Hi there<br><a href=/smart_chart>Smart chart</a>")


def smart_chart(request):
    global filename
    if request.method == 'POST':
        if request.POST.get('postname') == 'FILEPOST':
            form = f.UploadFileForm(request.POST, request.FILES)
            # comment right before first + sign if we are OK! with file overwrite for every session
            filename = request.POST.get('csrfmiddlewaretoken') + '-' + ''.join([i for j in range(5) for i in choice(ascii_letters) ])
            if form.is_valid():
                h.handle_uploaded_file(request.FILES['file'], filename)
                feedback = 'Success upload'
            else:
                feedback = 'Error uploading'
                filename = 0
        elif request.POST.get('postname') == 'CLASSIFY':
            if filename != 0:
                feedback = 'Succsessfully plotting'
                form = f.UploadFileForm()
                filename = 0
            else:
                feedback = 'No data submitet for analisys'
                form = f.UploadFileForm()
                filename = 0
    else:
        form = f.UploadFileForm()
        filename = 0


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


