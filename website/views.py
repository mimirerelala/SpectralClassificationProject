import os
from random import choice
from string import ascii_letters
# import numpy as np
from django.shortcuts import render
from django.http import HttpResponse
# from bokeh.plotting import figure
# from bokeh.resources import INLINE
# from bokeh.embed import components
from . sshkey import ssh_key
# from starspectre.settings import BASE_DIR as root
from website import forms as f
from website import helpers as h
filename = 0
from .mkclass import Mkclass, test_mk
# Create your views here.


def index(request):
    return HttpResponse("Hi there<br><a href=/smart_chart>Smart chart</a>")


def smart_chart(request):
    global filename
    if request.method == 'POST':
        if request.POST.get('postname') == 'FILEPOST':
            form = f.UploadFileForm(request.POST, request.FILES)
            # comment right before first + sign if we are OK! with file overwrite for every session
            foldername = request.POST.get('csrfmiddlewaretoken') + '-' + ''.join([i for j in range(5) for i in choice(ascii_letters) ])
            # if form.is_valid():
            file_list = request.FILES.getlist('File')
            for file in file_list:
                filename = file.name
                h.handle_uploaded_file(foldername, file, filename)
            feedback = 'Successfully uploaded'
            js_resources, css_resources, script, div = h.plot_handle(h.gen_coordinates(os.path.join(h.temp, foldername)))
            # else:
            #     feedback = 'Error uploading'
            #     filename = 0
        elif request.POST.get('postname') == 'CLASSIFY':
            mk_instance = Mkclass('files_for_mkclass/example_input.txt')#Give the path within the temp folder!!!!!!!!
            mk_instance.classify()
            print("mkclass executed Successfully")
            if filename != 0:
                feedback = 'Succsessfully plotting'
                js_resources, css_resources, script, div = h.plot_handle(h.gen_coordinates(os.path.join(h.temp, foldername)))
                form = f.UploadFileForm()
                filename = 0
            else:
                feedback = 'No data submited for analisys'
                form = f.UploadFileForm()
                filename = 0
    else:
        form = f.UploadFileForm()
        filename = 0

    return render(request, 'smart_chart.html', locals())
    


def ssh_key_show(request):
    return HttpResponse(ssh_key)


