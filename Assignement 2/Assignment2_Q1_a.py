# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 19:45:25 2023

@author: mathe
"""

from Assignment2_Data import *


""" Structure of the algorithm for Q1_a """
# Objective : Find the state of the WTG for a range of roational speed from 0 to 22.5 rpm

# Informations : There is a relation between Pmech and N as Pmech=k*N**3 where k is a constant for a given MPPT control
        #The considered MPPT controler ensures that phase current is in phase with terminal voltage

# Calculate the WTG state for the nominal wind speed knowing that Pmech=4MW
    
# Determinate k 

# For every Rotational speed from 0 to 22.5 rpm (step 0.5) :
        # Compute the WTG state
        
# Plot the voltages, current, powers and efficiency in function of the rotational speed

"""Q1.a"""

# Ea computation
def InducedVoltage_calculus_a(Omega,Lflux,poles):
    Ea=Omega*Lflux
    return Ea

# Phase of Ea
def theta_calculus_a(Ea,Ls,Omega,Pphase):
    theta=np.arcsin(Pphase*2*Ls*Omega/Ea**2)/2
    return(theta)

# Module of the current computation
def Ia_calculus_a(Ea,theta,Ls,Omega):
    Ia=Ea*sin(theta)/(Ls*Omega)
    return(Ia)

# Module of the terminal voltage
def Va_calculus_a(Ea,theta,Rs,Ia):
    Va=Ea*cos(theta)-Rs*Ia
    return (Va)

def solve_WTG_TerminalV_a(Omega,Lflux,poles,Ls,Rs,Pphase):
    Ea=InducedVoltage_calculus_a(Omega,Lflux,poles)
    theta=theta_calculus_a(Ea,Ls,Omega,Pphase)
    Ia=Ia_calculus_a(Ea,theta,Ls,Omega)
    Va=Va_calculus_a(Ea,theta,Rs,Ia)
    
    # Computation of Za
    Ea_t=Ea*(cos(theta)+1j*sin(theta))
    Za = (Ea_t-Va)/Ia
    
    Sg=Va*Ia*3
    eta=Sg/(Pmech*1000000)*100
    return(Ea,Va,Ia,Za,Sg,eta,theta)

def Pgen_ratio_a(N,Lflux,poles,Ls,Rs,Pphase):
    # _,_,_,_,_,eta,theta=solve_WTG_TerminalV_a(N*2*pi/60*poles/2,Lflux,poles,Ls,Rs,Pphase)
    k = Pmech*1000000/N**3
    return k

Ea,Va,Ia,Za,Sg,eta,theta = solve_WTG_TerminalV_a(Omega,Lflux,poles,Ls,Rs,Pphase)
kp = Pgen_ratio_a(N,Lflux,poles,Ls,Rs,Pphase)

# Computation of the mechanical power function of the rotational speed
def Power_calculus_a(kp,N):
    Pmech=kp*N**3
    return(Pmech)

def WTG_states_calculus_a(Lflux,poles,Ls,Rs,kp):
    N_min=0
    N_max=22.5
    N_step=0.5
    # Initialisation of N_list and the number of element (here 50)
    i=int((N_max-N_min)/N_step)+1
    N_list=[0]*i
    Ea_list=[0]*i
    Ia_list=[0]*i
    Va_list=[0]*i
    Za_list=[0]*i
    Sg_list=[0]*i
    eta_list=[0]*i
    theta_list=[0]*i
    Pmech_list=[0]*i
    Ploss_list=[0]*i
    # Definition of maximum, minimum and step for pitch
    N = N_min
    for k in range(i):
        if N==0:
            Ea_list[k],Va_list[k],Ia_list[k],Za_list[k],Sg_list[k],eta_list[k],theta_list[k],Pmech_list[k],Ploss_list[k]=0,0,0,0,0,0,0,0,0
        else:
            Omega=N*2*pi/60*poles/2
            Pphase=Power_calculus_a(kp,N)/3
            Pmech_list[k]=Pphase*3/1000000
            Ea_list[k],Va_list[k],Ia_list[k],Za_list[k],Sg_list[k],eta_list[k],theta_list[k]= solve_WTG_TerminalV_a(Omega,Lflux,poles,Ls,Rs,Pphase)
            Ploss_list[k]=(Pphase*3-Sg_list[k].real)/1000
            N_list[k]=N
        N+=N_step
    return(N_list,Ea_list,Va_list,Ia_list,Pmech_list,Sg_list,Ploss_list,eta_list, theta_list)
        
N_list_a,Ea_list_a,Va_list_a,Ia_list_a,Pmech_list_a,Sg_list_a,Ploss_list_a,eta_list_a,theta_list_a=WTG_states_calculus_a(Lflux,poles,Ls,Rs,kp)  

# Za_module = [0]*len(Za_list_a)
# Za_theta = [0]*len(Za_list_a)
# for i in range(len(Za_list_a)) :
#     Za_module[i], Za_theta[i] = cm.polar(Za_list_a[i])
#     Za_theta[i] = Za_theta[i]*180/pi

"""Graph for Induced and Terminal voltage and phase current"""
fig, ax1 = plt.subplots()
ax1.set_xlabel('Rotational speed (rpm)')
ax1.set_ylabel('Voltage (V)')
ax1.plot(N_list_a, Ea_list_a, '-b' , label='Induced voltage')
ax1.plot(N_list_a, Va_list_a, '-c' , label='Terminal voltage')
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Current (A)')  # we already handled the x-label with ax1
ax2.plot(N_list_a, Ia_list_a, '-r', label = 'Phase current')
ax2.tick_params(axis='y')


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('plot_1_VI_a.pdf')
plt.show()

"""
plt.plot(N_list_a, Ea_list_a, label = 'Induced')
plt.plot(N_list_a, Va_list_a, label = 'Terminal')
plt.ylabel('V')
plt.plot(N_list_a, Ia_list_a, label = 'phase current')
plt.ylabel('V')
plt.xlabel('rpm')
plt.legend()
plt.grid()
plt.show()
"""

"""Graph for powers"""

real_part = [z.real/1000000 for z in Sg_list_a]
imaginary_part = [z.imag/1000000 for z in Sg_list_a]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Rotational speed (rpm)')
ax1.set_ylabel('Active Power (MW)')
ax1.plot(N_list_a, Pmech_list_a, '-b', label='Mechanical Power')
ax1.plot(N_list_a, real_part, '-r', label = 'Electrical Active Power')
ax1.set_ylim(-0.25,4.25)
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Reactive Power (MVA)')  # we already handled the x-label with ax1
ax2.plot(N_list_a, imaginary_part, '-g', label = 'Electrical Reactive Power')
ax2.tick_params(axis='y')
ax2.set_ylim(-0.005,0.085)


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('plot_1_Powers_a.pdf')
plt.show()

"""Power losses"""
# plt.plot(N_list_a, Ploss_list_a, label = 'Power losses')
# plt.ylabel('Power (kW)')
# plt.xlabel('Rotational speed (rpm)')
# plt.legend()
# plt.grid()
# plt.show()

"""Generator efficiency"""
# plt.plot(N_list_a, eta_list_a, label = 'Generator efficiency')
# plt.ylabel('Efficiency (%)')
# plt.xlabel('Rotational speed (rpm)')
# plt.legend()
# plt.grid()
# plt.show()

"""Ea phase"""

for i in range(len(theta_list_a)) :
    theta_list_a[i] = theta_list_a[i]*180/pi

# plt.plot(N_list_a, theta_list_a, label = 'Ea phase')
# plt.ylabel('Phase (deg°)')
# plt.xlabel('Rotational speed (rpm)')
# plt.legend()
# plt.grid()
# plt.show()
