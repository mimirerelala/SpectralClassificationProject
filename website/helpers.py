import os
from starspectre.settings import BASE_DIR as root
temp = os.path.join(root, 'website/temp')


def handle_uploaded_file(file, filename):
    with open(os.path.join(temp, filename), 'wb+') as dataout:
        for chunk in file.chunks():
            dataout.write(chunk)