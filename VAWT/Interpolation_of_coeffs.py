# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:56:42 2022

@author: cgrinde
"""


""" This file is used to search for the right value of Cl and Cd by doing a double interpolation on angle of attach and thickness """

import csv

# The path were the different blades profils are
path = 'C:/Users/theod/Documents/Documents importants/DTU/Courses/46300 - Wind Turbine and Aerodynamics/For exam/VAWT/'

# We charge the different document
file=[path + 'airfoil.txt']

#Initializing tables were we are going to store our values from blades profils
cl_list=[0]*95
cd_list=[0]*95
aoa_list=[0]*95



with open(file[0]) as f:
    next(f)
    i=0
    for line in f:
        # Split each line into three values
        values = line.split()   
        # Convert values to float and append to respective lists
        aoa_list[i]=float(values[0])
        cl_list[i]=float(values[1])
        cd_list[i]=float(values[2])
        i+=1
# Print the lists
#print("Angle of Attack (AOA):", aoa_list)
#print("Coefficient of Lift (CL):", cl_list)
#print("Coefficient of Drag (CD):", cd_list)