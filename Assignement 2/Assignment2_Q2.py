# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 13:33:51 2023

@author: theodoreln
"""

from Assignment2_Data import *
from Assignment2_Q1_a import *
from Assignment2_Q1_b import *

""" Needed value """

P_Vsc = 3865310
V_Poc_goal = 33000*a/sqrt(3)


""" Functions for votages and currents computation """

# I_1 argument computation
def theta_1_compute(P_Vsc) :
    theta_1 = atan(-(0.2*P_Vsc)/P_Vsc)
    return(theta_1)

# theta_1 = theta_1_compute(P_Vsc)

# I_1 module 
def I_1_compute(P_Vsc, V_1, theta_1) :
    I_1 = P_Vsc/(3*V_1*cos(theta_1))
    return(I_1)

# I_1 = I_1_compute(P_Vsc, 492.17, theta_1)

# I_m computation
def I_m_computation(V_1, I_1, theta_1, Z1, Zm) : 
    I_m = (V_1-I_1*(cos(theta_1)+1j*sin(theta_1))*Z1)/Zm
    I_m, theta_m = cm.polar(I_m)
    return(I_m, theta_m)

# I_m, theta_m = I_m_computation(492.17, I_1, theta_1, Z1, Zm)

# I_2p computation
def I_2p_computation(I_1, theta_1, I_m, theta_m) :
    I_2p = I_1*(cos(theta_1)+1j*sin(theta_1))-I_m*(cos(theta_m)+1j*sin(theta_m))
    I_2p, theta_2pI = cm.polar(I_2p)
    return(I_2p, theta_2pI)

# I_2p, theta_2pI = I_2p_computation(I_1, theta_1, I_m, theta_m)

# V_2p computation
def V_2p_computation(I_m, theta_m, I_2p, theta_2pI, Zm, Z2p) :
    V_2p = I_m*(cos(theta_m)+1j*sin(theta_m))*Zm-I_2p*(cos(theta_2pI)+1j*sin(theta_2pI))*Z2p
    V_2p, theta_2pV = cm.polar(V_2p)
    return(V_2p, theta_2pV)

# V_2p, theta_2pV = V_2p_computation(I_m, theta_m, I_2p, theta_2pI, Zm, Z2p)

# V_cab computation
# def V_cab_computation(I_2p, theta_2pI, Zcable1) :
#     V_cab = I_2p*(cos(theta_2pI)+1j*sin(theta_2pI))*Zcable1
#     V_cab, theta_cab = cm.polar(V_cab)
#     return(V_cab, theta_cab)

# V_cab, theta_cab = V_cab_computation(I_2p, theta_2pI, Zcable1)

# V_Poc computation
def V_Poc_computation(V_2p, theta_2pV, I_2p, theta_2pI, Zcable1) : 
    V_Poc = V_2p*(cos(theta_2pV)+1j*sin(theta_2pV))-I_2p*(cos(theta_2pI)+1j*sin(theta_2pI))*Zcable1
    V_Poc, theta_PocV = cm.polar(V_Poc)
    return(V_Poc, theta_PocV)

# V_Poc, theta_PocV = V_Poc_computation(V_2p, theta_2pV, I_2p, theta_2pI, Zcable1)

# I_cap computation
def I_cap_computation(V_Poc, theta_PocV, Zcable2) :
    I_cap = V_Poc*(cos(theta_PocV)+1j*sin(theta_PocV))/Zcable2
    I_cap, theta_cap = cm.polar(I_cap)
    return(I_cap, theta_cap)

# I_cap, theta_capI = I_cap_computation(V_Poc, theta_PocV, Zcable2)

# I_Poc_computation
def I_Poc_computation(I_2p, theta_2pI, I_cap, theta_cap) :
    I_Poc = I_2p*(cos(theta_2pI)+1j*sin(theta_2pI))-I_cap*(cos(theta_cap)+1j*sin(theta_cap))
    I_Poc, theta_PocI = cm.polar(I_Poc)
    return(I_Poc, theta_PocI)

# I_Poc, theta_PocI = I_Poc_computation(I_2p, theta_2pI, I_cap, theta_capI)

# S_2p computation
def S_computation(V, theta_V, I, theta_I) :
    S = 3*V*(cos(theta_V)+1j*sin(theta_V))*I*(cos(theta_I)-1j*sin(theta_I))
    S, theta_S = cm.polar(S)
    return(S, theta_S)

# S_2p, theta_2pS = S_computation(V_2p, theta_2pV, I_2p, theta_2pI)
# S_cap, theta_capS = S_computation(V_Poc, theta_PocV, I_cap, theta_capI)
# S_Poc, theta_PocS = S_computation(V_Poc, theta_PocV, I_Poc, theta_PocI)

# Eta computation
def Eta_computation(P_Vsc, S_Poc, theta_PocS) :
    Eta = (S_Poc*cos(theta_PocS))/P_Vsc
    return(Eta)

# Eta = Eta_computation(P_Vsc, S_Poc, theta_PocS)

"""  Iteration on V_1 to find the right V_Poc """

def Q2_solve(P_Vsc, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2) :
    Solution = {}
    stop_criteria = 0
    error_goal = 0.001  
    V_1 = 0.001    
    while stop_criteria != 1 :
        # Status 
        # print('V_1 = ', V_1)
        # Computation of the all electrical system
        theta_1 = theta_1_compute(P_Vsc)
        I_1 = I_1_compute(P_Vsc, V_1, theta_1)
        I_m, theta_m = I_m_computation(V_1, I_1, theta_1, Z1, Zm)
        I_2p, theta_2pI = I_2p_computation(I_1, theta_1, I_m, theta_m)
        V_2p, theta_2pV = V_2p_computation(I_m, theta_m, I_2p, theta_2pI, Zm, Z2p)
        V_Poc, theta_PocV = V_Poc_computation(V_2p, theta_2pV, I_2p, theta_2pI, Zcable1)
        I_cap, theta_capI = I_cap_computation(V_Poc, theta_PocV, Zcable2)
        I_Poc, theta_PocI = I_Poc_computation(I_2p, theta_2pI, I_cap, theta_capI)
        S_2p, theta_2pS = S_computation(V_2p, theta_2pV, I_2p, theta_2pI)
        S_cap, theta_capS = S_computation(V_Poc, theta_PocV, I_cap, theta_capI)
        S_Poc, theta_PocS = S_computation(V_Poc, theta_PocV, I_Poc, theta_PocI)
        Eta = Eta_computation(P_Vsc, S_Poc, theta_PocS)
        
        # Condition for the error goal
        if abs(V_Poc-V_Poc_goal) < error_goal :
            Solution['V_1'] = V_1
            Solution['I_1'] = I_1
            Solution['theta_1'] = theta_1
            Solution['I_m'] = I_m
            Solution['theta_m'] = theta_m
            Solution['I_2p'] = I_2p
            Solution['theta_2pI'] = theta_2pI
            Solution['V_2p'] = V_2p
            Solution['theta_2pV'] = theta_2pV
            Solution['V_Poc'] = V_Poc
            Solution['theta_PocV'] = theta_PocV
            Solution['I_cap'] = I_cap
            Solution['theta_capI'] = theta_capI
            Solution['I_Poc'] = I_Poc
            Solution['theta_PocI'] = theta_PocI   
            Solution['S_2p'] = S_2p
            Solution['theta_2pS'] = theta_2pS
            Solution['S_cap'] = S_cap
            Solution['theta_capS'] = theta_capS
            Solution['S_Poc'] = S_Poc
            Solution['theta_PocS'] = theta_PocS
            Solution['Eta'] = Eta
            
            # To be 100% sure the solution is relevant
            if abs(theta_PocV*180/pi)<90:
                # Verification
                # print('V_1 = ', cm.rect(V_1,0))
                # print('I_1 = ', cm.rect(I_1,theta_1))
                # print('I_m = ', cm.rect(I_m,theta_m))
                # print('I_2p = ', cm.rect(I_2p,theta_2pI))
                # print('V_2p = ', cm.rect(V_2p,theta_2pV))
                # print('I_cap = ', cm.rect(I_cap,theta_capI))
                # print('I_Poc = ', cm.rect(I_Poc,theta_PocI))
                # print('V_Poc = ', cm.rect(V_Poc,theta_PocV))
                # print('S_2p = ', cm.rect(S_2p,theta_2pS))
                # print('S_cap = ', cm.rect(S_cap,theta_capS))
                # print('S_Poc = ', cm.rect(S_Poc,theta_PocS))
                # print('Eta = ', Eta)
                # Stop criteria
                stop_criteria = 1
        
        # Iteration 
        V_1 += 0.001
        
    return(Solution)

# Solution_a = Q2_solve(3865314, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2)
# Solution_b = Q2_solve(3915910, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2)


""" Algorithm for the iteration on the rotational speed """

def Q2_solve_list (Sg_list, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2) :
    N_min=0
    N_max=22.5
    N_step=0.5
    
    # Initialisation of N_list and the number of element (here 50)
    i=int((N_max-N_min)/N_step)+1
    N_list=[0]*i
    V_1_list = [0]*i
    I_1_list = [0]*i
    theta_1_list = [0]*i
    I_m_list = [0]*i
    theta_m_list = [0]*i
    I_2p_list = [0]*i
    theta_2pI_list = [0]*i
    V_2p_list = [0]*i
    theta_2pV_list = [0]*i
    V_Poc_list = [0]*i
    theta_PocV_list = [0]*i
    I_cap_list = [0]*i
    theta_capI_list = [0]*i
    I_Poc_list = [0]*i
    theta_Poc_list = [0]*i
    S_2p_list = [0]*i
    theta_2pS_list = [0]*i
    S_cap_list = [0]*i
    theta_capS_list = [0]*i
    S_Poc_list = [0]*i
    theta_PocS_list = [0]*i
    Eta_list = [0]*i
    
    # Definition of maximum, minimum and step for pitch
    N = N_min
    for k in range(i):
        print(k/i*100,'%')
        if N==0:
            N_list[k],V_1_list[k],I_1_list[k],theta_1_list[k],I_m_list[k],theta_m_list[k],I_2p_list[k],theta_2pI_list[k],V_2p_list[k],theta_2pV_list[k],V_Poc_list[k],theta_PocV_list[k],I_cap_list[k],theta_capI_list[k],I_Poc_list[k],theta_Poc_list[k],S_2p_list[k],theta_2pS_list[k],S_cap_list[k],theta_capS_list[k],S_Poc_list[k],theta_PocS_list[k],Eta_list[k] = 0,397.8359999975043,0,0,0,0,0,0,398.1569690477785,0,398.37069640903763,0,0,0,0,0,60495.073850544235,-1.5776727220290727,0,0,1992120.0235204976,1.5710558058837207,0
            
        else:
            Solution = Q2_solve(Sg_list[k], V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2)
            N_list[k]=N
            V_1_list[k],I_1_list[k],theta_1_list[k],I_m_list[k],theta_m_list[k],I_2p_list[k],theta_2pI_list[k],V_2p_list[k],theta_2pV_list[k],V_Poc_list[k],theta_PocV_list[k],I_cap_list[k],theta_capI_list[k],I_Poc_list[k],theta_Poc_list[k],S_2p_list[k],theta_2pS_list[k],S_cap_list[k],theta_capS_list[k],S_Poc_list[k],theta_PocS_list[k],Eta_list[k] = Solution['V_1'],Solution['I_1'],Solution['theta_1'],Solution['I_m'],Solution['theta_m'],Solution['I_2p'],Solution['theta_2pI'],Solution['V_2p'],Solution['theta_2pV'],Solution['V_Poc'],Solution['theta_PocV'],Solution['I_cap'],Solution['theta_capI'],Solution['I_Poc'],Solution['theta_PocI'],Solution['S_2p'],Solution['theta_2pS'],Solution['S_cap'],Solution['theta_capS'],Solution['S_Poc'],Solution['theta_PocS'],Solution['Eta']
            
            
        N+=N_step
        
    return(N_list,V_1_list,I_1_list,theta_1_list,I_m_list,theta_m_list,I_2p_list,theta_2pI_list,V_2p_list,theta_2pV_list,V_Poc_list,theta_PocV_list,I_cap_list,theta_capI_list,I_Poc_list,theta_Poc_list,S_2p_list,theta_2pS_list,S_cap_list,theta_capS_list,S_Poc_list,theta_PocS_list,Eta_list)

N_list_a,V_1_list_a,I_1_list_a,theta_1_list_a,I_m_list_a,theta_m_list_a,I_2p_list_a,theta_2pI_list_a,V_2p_list_a,theta_2pV_list_a,V_Poc_list_a,theta_PocV_list_a,I_cap_list_a,theta_capI_list_a,I_Poc_list_a,theta_Poc_list_a,S_2p_list_a,theta_2pS_list_a,S_cap_list_a,theta_capS_list_a,S_Poc_list_a,theta_PocS_list_a,Eta_list_a = Q2_solve_list(Sg_list_a, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2)

P_1_list_b = []
for i in range(len(Sg_list_b)) :
    P_1_list_b.append(Sg_list_b[i].real)
    
N_list_b,V_1_list_b,I_1_list_b,theta_1_list_b,I_m_list_b,theta_m_list_b,I_2p_list_b,theta_2pI_list_b,V_2p_list_b,theta_2pV_list_b,V_Poc_list_b,theta_PocV_list_b,I_cap_list_b,theta_capI_list_b,I_Poc_list_b,theta_Poc_list_b,S_2p_list_b,theta_2pS_list_b,S_cap_list_b,theta_capS_list_b,S_Poc_list_b,theta_PocS_list_b,Eta_list_b = Q2_solve_list(P_1_list_b, V_Poc_goal, Z1, Z2p, Zm, Zcable1, Zcable2)

""" Plotting for case A """ 

# Voltage amplitude
plt.plot(N_list_a, V_1_list_a, '-b', label = 'Transformer primary')
plt.plot(N_list_a, V_2p_list_a, '-r', label = 'Transformer secondary')
plt.plot(N_list_a, V_Poc_list_a, '-g', label = 'POC')
plt.ylabel('Voltage (V)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_V_a.pdf')
plt.show()

# Active Power
P_1_list_a = []
P_2p_list_a = []
P_Poc_list_a = []
for i in range(len(S_2p_list_a)) :
    P_1_list_a.append(Sg_list_a[i]/1000000)
    P_2p_list_a.append(S_2p_list_a[i]*cos(theta_2pS_list_a[i])/1000000)
    P_Poc_list_a.append(S_Poc_list_a[i]*cos(theta_PocS_list_a[i])/1000000)

plt.plot(N_list_a, P_1_list_a, '-b', label = 'Transformer primary')
plt.plot(N_list_a, P_2p_list_a, '-r', label = 'Transformer secondary')
plt.plot(N_list_a, P_Poc_list_a, '-g', label = 'POC')
plt.ylabel('Active Power (MW)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_PowerP_a.pdf')
plt.show()

# Reactive Power
Q_1_list_a = []
Q_2p_list_a = []
Q_Poc_list_a = []
for i in range(len(S_2p_list_a)) :
    Q_1_list_a.append(Sg_list_a[i]*0.2/1000000)
    Q_2p_list_a.append(S_2p_list_a[i]*sin(theta_2pS_list_a[i])/1000000)
    Q_Poc_list_a.append(S_Poc_list_a[i]*sin(theta_PocS_list_a[i])/1000000)

plt.plot(N_list_a, Q_1_list_a, '-b', label = 'Transformer primary')
plt.plot(N_list_a, Q_2p_list_a, '-r', label = 'Transformer secondary')
plt.plot(N_list_a, Q_Poc_list_a, '-g', label = 'POC')
plt.ylabel('Reactive Power (MVA)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_PowerQ_a.pdf')
plt.show()

# Efficiency 
eta_transfo_a = [0]
eta_cable_a = [0]
for i in range(len(Eta_list_a)) :
    # Eta_list_a[i] = Eta_list_a[i]*100
    if i != 0 :
        eta_transfo_a.append(P_2p_list_a[i]/P_1_list_a[i]*100)
        eta_cable_a.append(P_Poc_list_a[i]/P_2p_list_a[i]*100)
    
plt.plot(N_list_a, eta_transfo_a, '-b', label = 'Transformer')
plt.plot(N_list_a, eta_cable_a, '-r', label = 'Cable')
plt.plot(N_list_a, Eta_list_a, '-g',  label = 'Transmission System')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.ylim(0,100)
plt.legend()
plt.grid()
plt.savefig('plot_2_Eff_a.pdf')
plt.show()


""" Plotting for case B """ 

# Voltage amplitude
plt.plot(N_list_b, V_1_list_b, '-b', label = 'Transformer primary')
plt.plot(N_list_b, V_2p_list_b, '-r', label = 'Transformer secondary')
plt.plot(N_list_b, V_Poc_list_b, '-g', label = 'POC')
plt.ylabel('Voltage (V)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_V_b.pdf')
plt.show()

# Active Power
P_1_list_b = []
P_2p_list_b = []
P_Poc_list_b = []
for i in range(len(S_2p_list_b)) :
    P_1_list_b.append(Sg_list_b[i].real/1000000)
    P_2p_list_b.append(S_2p_list_b[i]*cos(theta_2pS_list_b[i])/1000000)
    P_Poc_list_b.append(S_Poc_list_b[i]*cos(theta_PocS_list_b[i])/1000000)

plt.plot(N_list_b, P_1_list_b, '-b', label = 'Transformer primary')
plt.plot(N_list_b, P_2p_list_b, '-r', label = 'Transformer secondary')
plt.plot(N_list_b, P_Poc_list_b, '-g', label = 'POC')
plt.ylabel('Active Power (MW)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_PowerP_b.pdf')
plt.show()

# Reactive Power
Q_1_list_b = []
Q_2p_list_b = []
Q_Poc_list_b = []
for i in range(len(S_2p_list_b)) :
    Q_1_list_b.append(Sg_list_b[i].real*0.2/1000000)
    Q_2p_list_b.append(S_2p_list_b[i]*sin(theta_2pS_list_b[i])/1000000)
    Q_Poc_list_b.append(S_Poc_list_b[i]*sin(theta_PocS_list_b[i])/1000000)

plt.plot(N_list_b, Q_1_list_b, '-b', label = 'Transformer primary')
plt.plot(N_list_b, Q_2p_list_b, '-r', label = 'Transformer secondary')
plt.plot(N_list_b, Q_Poc_list_b, '-g', label = 'POC')
plt.ylabel('Reactive Power (MVA)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_2_PowerQ_b.pdf')
plt.show()

# Efficiency 
eta_transfo_b = [0]
eta_cable_b = [0]
for i in range(len(Eta_list_b)) :
    Eta_list_b[i] = Eta_list_b[i]*100
    if i != 0 :
        eta_transfo_b.append(P_2p_list_b[i]/P_1_list_b[i]*100)
        eta_cable_b.append(P_Poc_list_b[i]/P_2p_list_b[i]*100)
    
plt.plot(N_list_b, eta_transfo_b, '-b', label = 'Transformer')
plt.plot(N_list_b, eta_cable_b, '-r', label = 'Cable')
plt.plot(N_list_b, Eta_list_b, '-g', label = 'Transmission System')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.ylim(0,100)
plt.legend()
plt.grid()
plt.savefig('plot_2_Eff_b.pdf')
plt.show()


""" Q3 """

def Q3_solve(S_Poc_list, theta_PocS_list, N_list, kp) :
    Eff_total = [0]
    for i in range(len(N_list)) :
        if i != 0 :
            Eff_total.append(100*S_Poc_list[i]*cos(theta_PocS_list[i])/(kp*N_list[i]**3))
    return(Eff_total)

Eff_total_a = Q3_solve(S_Poc_list_a, theta_PocS_list_a, N_list_a, kp)
Eff_total_b = Q3_solve(S_Poc_list_b, theta_PocS_list_b, N_list_b, kp)


plt.plot(N_list_a, Eff_total_a, label = 'Efficiency WT')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.ylim(0,100)
plt.legend()
plt.grid()
plt.savefig('plot_3_Eff_a.pdf')
plt.show()

plt.plot(N_list_b, Eff_total_b, label = 'Efficiency WT')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.ylim(0,100)
plt.legend()
plt.grid()

plt.savefig('plot_3_Eff_b.pdf')
plt.show()


S1_list_a = []
S1_list_b = []
for i in range(len(Sg_list_a)) :
    S1_list_a.append(Sg_list_a[i]+1j*0.2*Sg_list_a[i])
    S1_list_b.append(Sg_list_b[i].real+1j*0.2*Sg_list_a[i].real)
    S1_list_a[i] = cm.polar(S1_list_a[i])
    S1_list_b[i] = cm.polar(S1_list_b[i])




