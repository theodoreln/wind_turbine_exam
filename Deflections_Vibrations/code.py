# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 10:20:04 2023

@author: theodoreln
"""

from math import pi, sin, cos
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import eig
from math import sqrt

""" Importation of all datas """

path = 'C:/Users/theod/Documents/Documents importants/DTU/Courses/46300 - Wind Turbine and Aerodynamics/For exam/Deflections/'

# Structural properties    
r, stru_pitch, dist_mass, EI_one, EI_two, twist = [], [], [], [], [], []

file = 'bladestruc.txt'
with open(path+file, 'r') as file:
    lines = file.readlines()

for line in lines:
    columns = line.split()
    r.append(float(columns[0]))
    stru_pitch.append(float(columns[1])*pi/180)
    dist_mass.append(float(columns[2]))
    EI_one.append(float(columns[3]))
    EI_two.append(float(columns[4]))
    twist.append(float(columns[5])*pi/180)
    
# Loads
pn_6, pt_6 = [], []
V0_6, pitch_6, omega_6 = 6, 0.896*pi/180, 0.6283

file = 'loads6.txt'
with open(path+file, 'r') as file:
    lines = file.readlines()

for line in lines[2:] :
    columns = line.split()
    pn_6.append(float(columns[1]))
    pt_6.append(float(columns[2]))

pn_11, pt_11 = [], []
V0_11, pitch_11, omega_11 = 11, 0*pi/180, 0.9253

file = 'loads11.txt'
with open(path+file, 'r') as file:
    lines = file.readlines()

for line in lines[2:] :
    columns = line.split()
    pn_11.append(float(columns[1]))
    pt_11.append(float(columns[2]))

pn_20, pt_20 = [], []
V0_20, pitch_20, omega_20 = 20, 17.35*pi/180, 1.0053

file = 'loads20.txt'
with open(path+file, 'r') as file:
    lines = file.readlines()

for line in lines[2:] :
    columns = line.split()
    pn_20.append(float(columns[1]))
    pt_20.append(float(columns[2]))


""" Internal Forces """

# Watch out, in this exercice, we want to have a look at tangeantial and normal deformation
# So we are using pt as py and pn as pz !!!!! 
# The same for every y and z in the rest of the exercice !!!

# Shear force and Moment - Watch out forces are in kN
def T_M(r, py, pz) :
    N = len(r)
    Ty, Tz, My, Mz = [0]*N, [0]*N, [0]*N, [0]*N
    for i in range(N-1) :
        Ty[N-(i+2)] = Ty[N-(i+1)] + 0.5*(py[N-(i+1)]+py[N-(i+2)])*(r[N-(i+1)]-r[N-(i+2)])
        Tz[N-(i+2)] = Tz[N-(i+1)] + 0.5*(pz[N-(i+1)]+pz[N-(i+2)])*(r[N-(i+1)]-r[N-(i+2)])
        My[N-(i+2)] = My[N-(i+1)] - Tz[N-(i+1)]*(r[N-(i+1)]-r[N-(i+2)]) - ((1/6)*pz[N-(i+2)]+(1/3)*pz[N-(i+1)])*(r[N-(i+1)]-r[N-(i+2)])*(r[N-(i+1)]-r[N-(i+2)])      
        Mz[N-(i+2)] = Mz[N-(i+1)] + Ty[N-(i+1)]*(r[N-(i+1)]-r[N-(i+2)]) + ((1/6)*py[N-(i+2)]+(1/3)*py[N-(i+1)])*(r[N-(i+1)]-r[N-(i+2)])*(r[N-(i+1)]-r[N-(i+2)])
    return (Ty, Tz, My, Mz)

Ty_6, Tz_6, My_6, Mz_6 = T_M(r, pt_6, pn_6)
Ty_11, Tz_11, My_11, Mz_11 = T_M(r, pt_11, pn_11)
Ty_20, Tz_20, My_20, Mz_20 = T_M(r, pt_20, pn_20)


""" Curvature """

# Curvature computation
def curvature(My, Mz, stru_pitch, twist, pitch, EI_one, EI_two) :
    N = len(My)
    M1, M2, k1, k2, ky, kz = [0]*N, [0]*N, [0]*N, [0]*N, [0]*N, [0]*N
    for i in range(N) :
        M1[i] = My[i]*cos(stru_pitch[i] + twist[i] + pitch) - Mz[i]*sin(stru_pitch[i] + twist[i] + pitch)
        M2[i] = Mz[i]*cos(stru_pitch[i] + twist[i] + pitch) + My[i]*sin(stru_pitch[i] + twist[i] + pitch)
        # Here we us N and not kN for the computation of the curvature
        k1[i] = M1[i]*1000/EI_one[i]
        k2[i] = M2[i]*1000/EI_two[i]
        ky[i] = (k1[i]*cos(stru_pitch[i] + twist[i] + pitch) + k2[i]*sin(stru_pitch[i] + twist[i] + pitch))
        kz[i] = (k2[i]*cos(stru_pitch[i] + twist[i] + pitch) - k1[i]*sin(stru_pitch[i] + twist[i] + pitch))
    return(ky, kz)

ky_6, kz_6 = curvature(My_6, Mz_6, stru_pitch, twist, pitch_6, EI_one, EI_two)
ky_11, kz_11 = curvature(My_11, Mz_11, stru_pitch, twist, pitch_11, EI_one, EI_two)
ky_20, kz_20 = curvature(My_20, Mz_20, stru_pitch, twist, pitch_20, EI_one, EI_two)


""" Angle and deflection """

def angle_defl(r, ky, kz) :
    N = len(r)
    angley, anglez, defly, deflz = [0]*N, [0]*N, [0]*N, [0]*N
    for i in range(N-1) :
        angley[i+1] = angley[i] + 0.5*(ky[i+1]+ky[i])*(r[i+1]-r[i])
        anglez[i+1] = anglez[i] + 0.5*(kz[i+1]+kz[i])*(r[i+1]-r[i])
        defly[i+1] = defly[i] + anglez[i]*(r[i+1]-r[i]) + ((1/6)*kz[i+1]+(1/3)*kz[i])*(r[i+1]-r[i])*(r[i+1]-r[i])
        deflz[i+1] = deflz[i] - angley[i]*(r[i+1]-r[i]) - ((1/6)*ky[i+1]+(1/3)*ky[i])*(r[i+1]-r[i])*(r[i+1]-r[i])
    return(angley, anglez, defly, deflz)

angley_6, anglez_6, defly_6, deflz_6 = angle_defl(r, ky_6, kz_6)
angley_11, anglez_11, defly_11, deflz_11 = angle_defl(r, ky_11, kz_11)
angley_20, anglez_20, defly_20, deflz_20 = angle_defl(r, ky_20, kz_20)
        

""" Plotting """

def plot(r, V0, py, pz, Ty, Tz, My, Mz, ky, kz, angley, anglez, defly, deflz) :
     
    plt.plot(r, py, 'b.-')
    plt.plot(r, pz, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('p(x) [kN/m]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["py", "pz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
    plt.plot(r, Ty, 'b.-')
    plt.plot(r, Tz, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('T(x) [kN]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["Ty", "Tz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
    plt.plot(r, My, 'b.-')
    plt.plot(r, Mz, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('M(x) [kNm]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["My", "Mz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
    plt.plot(r, ky, 'b.-')
    plt.plot(r, kz, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('k(x) [rad/m]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["ky", "kz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
    plt.plot(r, angley, 'b.-')
    plt.plot(r, anglez, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('theta(x) [rad]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["thetay", "thetaz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
    plt.plot(r, defly, 'b.-')
    plt.plot(r, deflz, 'r.-')
    plt.xlabel('x [m]')
    plt.ylabel('u(x) [m]')
    plt.title('V0 = ' + str(V0) + ' m/s')
    plt.legend(["uy", "uz"]) 
    plt.axhline(0, color='gray', linewidth=0.5)
    plt.show()
    
plot(r, V0_6, pt_6, pn_6, Ty_6, Tz_6, My_6, Mz_6, ky_6, kz_6, angley_6, anglez_6, defly_6, deflz_6)
plot(r, V0_11, pt_11, pn_11, Ty_11, Tz_11, My_11, Mz_11, ky_11, kz_11, angley_11, anglez_11, defly_11, deflz_11)
plot(r, V0_20, pt_20, pn_20, Ty_20, Tz_20, My_20, Mz_20, ky_20, kz_20, angley_20, anglez_20, defly_20, deflz_20)

plt.plot(r, pn_6, 'b.-')
plt.plot(r, pn_11, 'r.-')
plt.plot(r, pn_20, 'y.-')
plt.xlabel('x [m]')
plt.ylabel('Pz(x) [kN/m]')
plt.legend(["V0 = 6 m/s", "V0 = 11 m/s", "V0 = 20 m/s"]) 
plt.show()

plt.plot(r, pt_6, 'b.-')
plt.plot(r, pt_11, 'r.-')
plt.plot(r, pt_20, 'y.-')
plt.xlabel('x [m]')
plt.ylabel('Py(x) [kN/m]')
plt.legend(["V0 = 6 m/s", "V0 = 11 m/s", "V0 = 20 m/s"]) 
plt.show()


""" Deflection calculation from the beginning """

def defl(py, pz, r, EI_one, EI_two, stru_pitch, twist, pitch) :
    Ty, Tz, My, Mz = T_M(r, py, pz)
    ky, kz = curvature(My, Mz, stru_pitch, twist, pitch, EI_one, EI_two)
    angley, anglez, uy, uz = angle_defl(r, ky, kz)
    return (uy, uz)


""" Flexibility matrix """

# Computation of the flexibility matrix
def Flexibility(r, EI_one, EI_two, stru_pitch, twist, pitch) :
    N = len(r)
    F = np.empty(((N-1)*2,(N-1)*2))
    for i in range((N-1)*2) :
        P = [0]*(N-1)*2
        # Here 0.001 and not 1 because loads supposed to be in kN !!!
        P[i] = 0.001
        py = [0]+P[:N-1]
        pz = [0]+P[N-1:]
        uy, uz = defl(py, pz, r, EI_one, EI_two, stru_pitch, twist, pitch)
        for j in range(N-1) :
            F[j][i] = uy[j+1]
            F[j+N-1][i] = uz[j+1]
    return (F)

F = Flexibility(r, EI_one, EI_two, stru_pitch, twist, pitch_11)

""" Mass matrix """

def Mass(r, dist_mass) :
    N = len(r)
    M = np.zeros(((N-1)*2,(N-1)*2))
    for i in range(N-1) :
        M[i][i] = dist_mass[i+1]
        M[i+N-1][i+N-1] = dist_mass[i+1]
    return(M)

M = Mass(r, dist_mass)

""" Solve the eigen values problem """

def Solve(F, M) :
    Matrix = np.matmul(F,M)
    V,D = eig(Matrix)
    return(V,D)

V, D = Solve(F,M)
omega1 = sqrt(1/V[0])
omega2 = sqrt(1/V[1])
omega3 = sqrt(1/V[2])

f1 = omega1/(2*pi)
f2 = omega2/(2*pi)
f3 = omega3/(2*pi)

""" Matrix with mode shapes """

M1y, M1z, M2y, M2z, M3y, M3z = [0], [0], [0], [0], [0], [0],
n = len(r)
for i in range(n-1):
    M1y.append(D[i,0]/np.max(np.abs(D[:,0])))
    M1z.append(D[n-1+i,0]/np.max(np.abs(D[:,0])))
    M2y.append(D[i,1]/np.max(np.abs(D[:,1])))
    M2z.append(D[n-1+i,1]/np.max(np.abs(D[:,1])))
    M3y.append(D[i,2]/np.max(np.abs(D[:,2])))
    M3z.append(D[n-1+i,2]/np.max(np.abs(D[:,2])))
    
plt.plot(r, M1y, '-b', label='M1y')
plt.plot(r, M1z, '-r', label='M1z')
plt.legend()
plt.grid()
plt.show()

plt.plot(r, M2y, '-b', label='M2y')
plt.plot(r, M2z, '-r', label='M2z')
plt.legend()
plt.grid()
plt.show()

plt.plot(r, M3y, '-b', label='M3y')
plt.plot(r, M3z, '-r', label='M3z')
plt.legend()
plt.grid()
plt.show()



