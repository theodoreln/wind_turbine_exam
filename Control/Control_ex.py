import matplotlib.pyplot as plt
import numpy as np
from math import *

I = 15000

""" Read information about Cp(lambda) and V0(t) """

path = 'C:/Users/theod/Documents/Documents importants/DTU/Courses/46300 - Wind Turbine and Aerodynamics/For exam/Control/'
file = [path + 'Cp_Lambda.txt']
lambda_list = []
Cp_list = []

with open(file[0]) as f:
    for line in f :
        values = line.split()
        if len(values) == 2:
            col1, col2 = map(float, values)
            lambda_list.append(col1)
            Cp_list.append(col2)
            

file = [path + 'V0.txt']
t_list = []
V0_list = []

with open(file[0]) as f:
    for line in f :
        values = line.split()
        if len(values) == 2:
            col1, col2 = map(float, values)
            t_list.append(col1)
            V0_list.append(col2)
            
""" Plotting of Cp and V0 """

plt.plot(lambda_list, Cp_list, label='Cp')
plt.legend()
plt.grid()
plt.show()

plt.plot(t_list, V0_list, label='V0')
plt.legend()
plt.grid()
plt.show()

""" Function to calculate w, Mr, Mg """

# Interpolation of Cp on lambda
def interpolation_Cp(lamb, lambda_list, Cp_list) :
    Cp = np.interp(lamb, lambda_list, Cp_list)
    return (Cp)

def solve(lambda_list, Cp_list, t_list, V0_list, w_initial, V_initial) :
    n=len(t_list)
    w_list = []
    Mr_list = []
    Mg_list = []
    for i in range(n) :
        if i == 0 :
            w = w_initial
            V = V_initial
            lamb = w*9.5/V
        else :
            w = w_list[i-1]
            V = V0_list[i-1]
            lamb = w*9.5/V
        Cp = np.interp(lamb, lambda_list, Cp_list)
        Mr = (0.5*1.225*V*V*V*pi*9.5*9.5*Cp)/w
        if w < 5 :
            Mg = 0
        else :
            Mg = 255000*(w-5)
        w = w+(Mr-Mg)*(0.1)/I
        w_list.append(w)
        if i != 0 :
            Mr_list.append(Mr)
            Mg_list.append(Mg)
        if i == n-1 :
            w = w_list[i]
            V = V0_list[i]
            lamb = w*9.5/V
            Cp = np.interp(lamb, lambda_list, Cp_list)
            Mr = (0.5*1.225*V*V*V*pi*9.5*9.5*Cp)/w
            if w < 5 :
                Mg = 0
            else :
                Mg = 255000*(w-5)
            Mr_list.append(Mr)
            Mg_list.append(Mg)
    return(w_list, Mr_list, Mg_list)
        

w_list, Mr_list, Mg_list = solve(lambda_list, Cp_list, t_list, [10]*len(t_list), 1, 10)
plt.plot(t_list, w_list, label='w')
plt.legend()
plt.grid()
plt.show()

plt.plot(t_list, Mr_list, label='Mr')
plt.plot(t_list, Mg_list, label='Mg')
plt.legend()
plt.grid()
plt.show()

w_list, Mr_list, Mg_list = solve(lambda_list, Cp_list, t_list, V0_list, 1, 10)
plt.plot(t_list, w_list, label='w')
plt.legend()
plt.grid()
plt.show()

plt.plot(t_list, Mr_list, label='Mr')
plt.plot(t_list, Mg_list, label='Mg')
plt.legend()
plt.grid()
plt.show()


""" Free turbine """

def free(lambda_list, Cp_list, t_list, V0_list, w_initial, V_initial) :
    n=len(t_list)
    w_list = []
    Mr_list = []
    Mg_list = []
    for i in range(n) :
        if i == 0 :
            w = w_initial
            V = V_initial
            lamb = w*9.5/V
        else :
            w = w_list[i-1]
            V = V0_list[i-1]
            lamb = w*9.5/V
        Cp = np.interp(lamb, lambda_list, Cp_list)
        Mr = (0.5*1.225*V*V*V*pi*9.5*9.5*Cp)/w
        w = w+(Mr)*(0.1)/I
        w_list.append(w)
        if i != 0 :
            Mr_list.append(Mr)
        if i == n-1 :
            w = w_list[i]
            V = V0_list[i]
            lamb = w*9.5/V
            Cp = np.interp(lamb, lambda_list, Cp_list)
            Mr = (0.5*1.225*V*V*V*pi*9.5*9.5*Cp)/w
            Mr_list.append(Mr)
    return(w_list, Mr_list)


w_list, Mr_list = free(lambda_list, Cp_list, t_list, [10]*len(t_list), 1, 10)
plt.plot(t_list, w_list, label='w')
plt.legend()
plt.grid()
plt.show()

plt.plot(t_list, Mr_list, label='Mr')
plt.legend()
plt.grid()
plt.show()








