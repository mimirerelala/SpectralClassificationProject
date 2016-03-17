import numpy as n
import json

list_templates = n.genfromtxt('libnor36_list', dtype = str)


temp_to_classes = {320: 'G5', 130: 'B7', 260: 'F5', 70: 'B0',  77:'B1',
455: 'M5', 445: 'M4', 200: 'A5', 330: 'G8', 140: 'B8', 270: 'F6', 80: 'B1', 85:'B1.5',
440: 'M3.5', 210: 'A7', 340: 'K0', 150: 'B9', 407: 'M0', 410:'M0', 280: 'F8', 
90: 'B2', 30: 'O6', 160: 'A0', 400: 'K7', 290: 'F9', 100: 'B3', 385:'K4.5', 390: 'K5', 
230: 'F0',235:'F1', 40: 'O7', 425: 'M2', 427:'M2', 170: 'A1', 300: 'G0', 370: 'K3', 240: 'F2', 
350: 'K1', 360: 'K2', 50: 'O8',307:'G1', 310: 'G2', 120: 'B5', 250: 'F3', 60: 'O9', 190: 'A3', 
474: 'M6', 495:'M8'}
#IMA IZMISLENI!!!!



lumins_formatting = {'00':'Ia','10':'Ib','15':'Ib','20':'II','25':'III','30':'III','35':'III','40':'IV','50':'V'}

classes = set()
lumins = set()

for template in list_templates:
	class_temp = template[1:4]
	lumin_class = template[5:7]
	classes.add(class_temp)
	lumins.add(lumin_class)
	print(template, temp_to_classes[int(class_temp)], lumins_formatting[lumin_class])





#print(classes)
#print(lumins)
#print(temp_to_classes)