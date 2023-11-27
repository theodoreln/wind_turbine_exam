# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 10:26:47 2023

@author: theodoreln
"""


""" This file just read the geometry of the blade and store it in memory """

import csv

path = 'C:/Users/theod/Documents/Documents importants/DTU/Courses/46300 - Wind Turbine and Aerodynamics/Assignements/'
file = [path + 'bladedat.txt']
n = 18
r_list = [0]*n
c_list = [0]*n
twist_list = [0]*n
thick_list = [0]*n

with open(file[0]) as f:
    reader = csv.reader(f, delimiter='\t')
    i = 0
    for row in reader:
        r_list[i] = float(row[0])
        twist_list[i] = float(row[1])
        c_list[i] = float(row[2])
        thick_list[i] = float(row[3])
        i = i + 1

