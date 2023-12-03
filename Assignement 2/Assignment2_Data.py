# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 19:45:25 2023

@author: mathe
"""

import cmath as cm
from math import*
import numpy as np
import matplotlib.pyplot as plt

Pmech = 4#MW
Pphase = Pmech/3*1000000 #W
Vwind = 12#m/s
N = 22.5#rpm
Vcutin = 3#m/s

"""Generator"""
poles = 52;
Omega=N*2*pi/60*poles/2
#Stator
Ls = 5.573/1000#H
Rs = 14.821/1000#Ohm
Lflux = 15.826#Wb

"""Wind turbine Transformer"""
R1 = 10/1000 #Ohm
R2p = 10/1000 #m#Ohm
Rc = 1240 #Ohm
L1 = 20/1000000 #microH
L2p = 20/1000000 #microH
Lm = 25/1000 #mH

a = 690/33000 #transformateur

"""Cable connection"""
r_cable = 0.5 #Ohm/km
l_cable = 0.5/1000 #mH/km
c_cable = 0.1/1000000 #microF/km

"""Grid"""
Nb_Turbines = 15
L_cable = 60 #km
Omega_grid = 2*pi*50

"""Impedence"""
Z1 = R1+1j*L1*Omega_grid
Z2p = R2p + 1j*L2p*Omega_grid
Zm = (1j*Rc*Lm*Omega_grid)/(Rc+1j*Lm*Omega_grid)
Zcable1 = (r_cable*L_cable + 1j*l_cable*L_cable*Omega_grid)*a*a
Zcable2 = (a*a)/(1j*c_cable*L_cable*Omega_grid)

