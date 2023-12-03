# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 19:45:25 2023

@author: mathe
"""

from Assignment2_Data import *
from Assignment2_Q1_a import *

""" Structure of the algorithm for Q1_b """
# Objective : Find the state of the WTG for a range of roational speed from 0 to 22.5 rpm

# Informations : There is a relation between Pmech and N as Pmech=k*N**3 where k is a constant for a given MPPT control
        #The considered MPPT controler ensures that the phase current is in phase with induce voltage 

# Calculate the WTG state for the nominal wind speed knowing that Pmech=4MW
    
# Determinate k

# For every Rotational speed from 0 to 22.5 rpm (step 0.5) :
        # Compute the WTG state
        
# Plot the voltages, current, powers and efficiency in function of the rotational speed


""" Computation of the Induced Voltage, the current, their phase, the terminal voltage """

# Calcul the Induced Voltage of one phase
def InducedVoltage_calculus_b(Omega,Lflux,poles):
    Ea=Omega*Lflux
    return Ea

# Calcul the Current of the phase
def Ia_calculus_b(Pphase,Ea):
    Ia=Pphase/Ea
    return(Ia)

# Calcul the phase of the terminal voltage
def thetacalculus_b(Ea,Ia,Rs,Ls,Omega):
    theta=np.arctan((Omega*Ls*Ia)/(Ea-Rs*Ia))
    return(theta)

# Calcul the module of the terminal voltage
def Va_calculus_b(Omega,Ia,Rs,Ea,Ls):
    Va=sqrt((Omega*Ls*Ia)**2+(Ea-Rs*Ia)**2)
    return (Va)

""" Computation of the state of the WTG for one rotational speed """
def solve_WTG_TerminalV_b(Omega,Lflux,poles,Ls,Rs,Pphase):
    Ea=InducedVoltage_calculus_b(Omega,Lflux,poles)
    Ia=Ia_calculus_b(Pphase,Ea)
    theta=thetacalculus_b(Ea,Ia,Rs,Ls,Omega)
    Va=Va_calculus_b(Omega,Ia,Rs,Ea,Ls)
    
    # Computation of Za
    Ea_t=Ea*(cos(theta)+1j*sin(theta))
    Ia_t=Ia*(cos(theta)+1j*sin(theta))
    Za=(Ea_t-Va)/Ia_t
    
    Sg=Va*Ia*(cos(theta)-1J*sin(theta))*3
    Ia=Ia*(cos(theta)+1j*sin(theta))
    eta=Sg.real/(Pmech*1000000)*100
    return(Ea,Va,Ia,Za,Sg,eta,theta)

#Calculus of the state of the WTG at the nominal Wind speed
Ea,Va,Ia,Za,Sg,eta,theta=solve_WTG_TerminalV_b(Omega,Lflux,poles,Ls,Rs,Pphase)


#Power ratio calculus
def Pgen_ratio_b (N,Lflux,poles,Ls,Rs,Pphase):
    k=Pmech*1000000/N**3
    return k

#Calculus of the power ratio of the system MPPT control of Q_b
kp=Pgen_ratio_b (N,Lflux,poles,Ls,Rs,Pphase)

#Power calculus
def Power_calculus_b (kp,N):
    Pmech=kp*N**3
    return(Pmech)

""" Algorithm for the iteration on the wind speed """

def WTG_states_calculus_b (Lflux,poles,Ls,Rs,kp):
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
    Pmech_list=[0]*i
    Ploss_list=[0]*i
    theta_list=[0]*i
    # Definition of maximum, minimum and step for pitch
    N = N_min
    for k in range(i):
        if N==0:
            Ea_list[k],Va_list[k],Ia_list[k],Za_list[k],Sg_list[k],eta_list[k],theta_list[k],Pmech_list[k],Ploss_list[k]=0,0,0,0,0,0,0,0,0
        else:
            Omega=N*2*pi/60*poles/2
            Pphase=Power_calculus_b(kp,N)/3
            Pmech_list[k]=Pphase*3/1000000
            Ea,Va,Ia,Za,Sg,eta,theta= solve_WTG_TerminalV_b(Omega,Lflux,poles,Ls,Rs,Pphase)
            Ea_list[k],Va_list[k],Ia_list[k],Za_list[k],Sg_list[k],eta_list[k],theta_list[k]=Ea,Va,abs(Ia),Za,Sg,eta,theta
            Ploss_list[k]=(Pphase*3-Sg_list[k].real)/1000
            N_list[k]=N
        N+=N_step
    return(N_list,Ea_list,Va_list,Ia_list,Pmech_list,Sg_list,Ploss_list,eta_list,theta_list)
        
N_list_b,Ea_list_b,Va_list_b,Ia_list_b,Pmech_list_b,Sg_list_b,Ploss_list_b,eta_list_b,theta_list_b=WTG_states_calculus_b (Lflux,poles,Ls,Rs,kp)  

# Za_module = [0]*len(Za_list_b)
# Za_theta = [0]*len(Za_list_b)
# Sg_module = [0]*len(Za_list_b)
# Sg_theta = [0]*len(Za_list_b)
# for i in range(len(Za_list_b)) :
#     Sg_module[i], Sg_theta[i] = cm.polar(Sg_list_b[i])
#     Sg_theta[i] = Sg_theta[i]*180/pi
#     Za_module[i], Za_theta[i] = cm.polar(Za_list_b[i])
#     Za_theta[i] = Za_theta[i]*180/pi

""" Graph for Induced and Terminal voltage and phase current """
#Method to plot 2 different scales in one graph
fig, ax1 = plt.subplots()
ax1.set_xlabel('Rotational speed (rpm)')
ax1.set_ylabel('Voltage (V)')
ax1.plot(N_list_b, Ea_list_b, '-b' , label='Induced voltage')
ax1.plot(N_list_b, Va_list_b, '-c' , label = 'Terminal voltage')
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Current (A)')  # we already handled the x-label with ax1
ax2.plot(N_list_b, Ia_list_b, '-r', label = 'Phase current')
ax2.tick_params(axis='y')


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('plot_1_VI_b.pdf')
plt.show()


""" Graph for powers """

real_part = [z.real/1000000 for z in Sg_list_b]
imaginary_part = [z.imag/1000000 for z in Sg_list_b]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Rotational speed (rpm)')
ax1.set_ylabel('Active Power (MW)')
ax1.plot(N_list_b, Pmech_list_b, '-b', label='Mechanical Power')
ax1.plot(N_list_b, real_part, '-r', label = 'Electrical Active Power')
ax1.set_ylim(-1,4.208)
ax1.grid()

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Reactive Power (MVA)')  # we already handled the x-label with ax1
ax2.plot(N_list_b, imaginary_part, '-g', label = 'Electrical Reactive Power')
ax2.tick_params(axis='y')
ax2.set_ylim(-2.5,10.5)


lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('plot_1_Powers_b.pdf')
plt.show()


""" Power losses """

plt.plot(N_list_a, Ploss_list_a, '-b', label = 'Power loss case a')
plt.plot(N_list_b, Ploss_list_b, '-r', label = 'Power loss case b')
plt.ylabel('Power (kW)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_1_PowerLoss.pdf')
plt.show()


""" Generator efficiency """

plt.plot(N_list_a, eta_list_a, '-b', label = 'Generator efficiency case a')
plt.plot(N_list_b, eta_list_b, '-r', label = 'Generator efficiency case b')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_1_Efficiency.pdf')
plt.show()


plt.plot(N_list_a, eta_list_a, '-b', label = 'Generator efficiency case a')
plt.plot(N_list_b, eta_list_b, '-r', label = 'Generator efficiency case b')
plt.ylabel('Efficiency (%)')
plt.xlabel('Rotational speed (rpm)')
plt.xlim(20,23)
plt.ylim(78,102)
plt.legend()
plt.grid()
plt.savefig('plot_1_Efficiency_zoom.pdf')
plt.show()


"""Ea phase"""

for i in range(len(theta_list_b)) :
    theta_list_b[i] = theta_list_b[i]*180/pi

plt.plot(N_list_a, theta_list_a, '-b', label='Ea phase case a')
plt.plot(N_list_b, theta_list_b, '-r', label = 'Ea/Ia phase case b')
plt.ylabel('Phase (deg°)')
plt.xlabel('Rotational speed (rpm)')
plt.legend()
plt.grid()
plt.savefig('plot_1_Phase.pdf')
plt.show()





