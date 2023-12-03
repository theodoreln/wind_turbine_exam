import matplotlib.pyplot as plt
from math import *

""" Importing electrical parameters -> Modify values """

# Factor of the transformer : a = V1/V2
a = 9

# Frequency 
omega = 50*2*pi

# Primary side leakage impedance
R1 = 4.52
L1 = 0.0169
Z1 = R1 + 1j*omega*L1

# Secondary side leakage impedance 
# !!! Watch out, transform the value from the secondary to the primary by multiplying by a^2 !!!
R2 = 4.05
L2 = 0.0134
Z2 = R2 + 1j*omega*L2

# Core resistance and magnetizing inductance
# !!! Watch out, be certain that the values are taken to the right side, here transformation to primary side !!!
Rc = 10.5
Lm = 0.05
Zm = Zm = a*a*(1j*omega*Rc*Lm)/(Rc+1j*omega*Lm)

# Power factor of the primary (ind)
PF1 = 0.8

# Voltage of the primary, phase(V1)=0°
V1 = 3600

# Amplitude of the primary current (can also give the power)
I1 = 50

""" Function for the computation of every value, taken to the primary !!!!! """

def configuration(V1, I1, PF1, Z1, Z2, Zm) :
    # Complexe current, + instead of - si PF capacitif
    I1 = I1*(PF1-1j*sin(acos(PF1)))
    Vm = V1-Z1*I1
    Im = Vm/Zm
    I2 = I1 - Im
    V2 = Vm-Z2*I2
    S1 = 3*V1*I1.conjugate()
    S2 = 3*V2*I2.conjugate()
    Eff = S2.real/S1.real
    return(I1, Vm, Im, I2, V2, S1, S2, Eff)

I1, Vm, Im, I2, V2, S1, S2, Eff = configuration(V1, I1, PF1, Z1, Z2, Zm)
