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
                feedback = 'Successfully uploaded'
                js_resources, css_resources, script, div = h.plot_handle(h.gen_coordinates(h.temp, filename))
            else:
                feedback = 'Error uploading'
                filename = 0
        elif request.POST.get('postname') == 'CLASSIFY':
            if filename != 0:
                feedback = 'Succsessfully plotting'
                js_resources, css_resources, script, div = h.plot_handle(h.gen_coordinates(h.temp, filename))
                #form = f.UploadFileForm()
                #filename = 0
            else:
                feedback = 'No data submitet for analisys'
                form = f.UploadFileForm()
                filename = 0
    else:
        form = f.UploadFileForm()
        filename = 0

    return render(request, 'smart_chart.html', locals())
    


def ssh_key_show(request):
    return HttpResponse(ssh_key)


