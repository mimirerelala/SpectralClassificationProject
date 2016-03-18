from django.db import models
from .template_model_data import template_model_data
# Create your models here.

from django.db import models


class Plot_data_point(models.Model):
    filename = models.CharField(max_length=500)
    x_data = models.DecimalField(max_digits = 10, decimal_places=2)
    y_data = models.DecimalField(max_digits = 10, decimal_places=6)

class Results_from_mkclass(models.Model):
	filename = models.CharField(max_length=500)
	resulting_class = models.CharField(max_length=100)
	resulting_quality = models.CharField(max_length=30)
	calc_date = models.DateTimeField()

class Template_choice(models.Model):
	template_types_file_names = models.CharField(max_length=2,
                                      choices=template_model_data,
                                      default='t160l50p00.rbn')

