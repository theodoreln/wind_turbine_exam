# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 10:18:51 2023

@author: mathe
"""
from math import *
from Interpolation_of_coeffs import *
import matplotlib.pyplot as plt
import numpy as np

#Calculus of the state for indice i
def state_calculus (R,theta_i,V0,W_x_i,W_y_i,Omega,rho,c,aoa_list, cl_list, cd_list,pitch):
    x_i=-R*sin(theta_i)
    y_i=R*cos(theta_i)
    Vrel_x_i=Omega*y_i+V0-W_x_i
    Vrel_y_i=-Omega*x_i+W_y_i
    Vnorm_i=(V0-W_x_i)*sin(theta_i)-W_y_i*cos(theta_i)
    Vtan_i=(V0-W_x_i)*cos(theta_i)+W_y_i*sin(theta_i)+Omega*R
    Vrel_i=sqrt(Vnorm_i**2+Vtan_i**2)
    # Convert in degrees because the Cl Cd table in degrees 
    alpha_i=np.rad2deg(atan(Vnorm_i/Vtan_i))
    Cl_alpha,Cd_alpha=interpolate_airfoil_data(alpha_i, aoa_list, cl_list, cd_list)
    l_i=1/2*rho*Vrel_i**2*c*Cl_alpha
    d_i=1/2*rho*Vrel_i**2*c*Cd_alpha
    cos_Beta_i=Vrel_y_i/Vrel_i
    sin_Beta_i=Vrel_x_i/Vrel_i
    p_x_i=l_i*cos_Beta_i+d_i*sin_Beta_i
    p_y_i=-l_i*sin_Beta_i+d_i*cos_Beta_i
    
    #needed for the calculus of the coefficient CP
    phi=(alpha_i+pitch)*pi/180
    p_n_i = l_i*cos(phi)+d_i*sin(phi)
    p_t_i = l_i*sin(phi)-d_i*cos(phi)
    
    return(p_x_i,p_y_i,p_n_i,p_t_i)

#Calculus of the W_x and W_y for indice i
def W_i_n_calculus(W_0,theta_i,R,L):
    W_x_i_n=W_0*(1-R/L*sin(theta_i))
    W_y_i_n=R/L*W_x_i_n*cos(theta_i)
    return (W_x_i_n,W_y_i_n)

#Evaluation of the cl and cd coefficient by interpolation of the values for a given angle of attack (alpha)
def interpolate_airfoil_data(aoa, aoa_list,cl_list,cd_list):
    cl = np.interp(aoa, aoa_list, cl_list)
    cd = np.interp(aoa, aoa_list, cd_list)
    return cl, cd

#Unstaedy Bem algorythm using ntime which is the number of periods for which we calculate the state
def unsteady_Bem (ntime,Delta_t,Omega,V0,B,R,c,L,pitch):
    #Definition of the starting point
    a_n=0
    theta_1_n=0
    
    #creation of the list that we will return in order to plot the final results
    P_x_t_list=[0]*(ntime)
    P_y_t_list=[0]*(ntime)
    t_list=[0]*(ntime)
    CT_list=[0]*(ntime)
    CP_list=[0]*(ntime)
    
    # Creation of the list for the first blade (useful to see the difference of phase when moving the number of blades)
    P_x_t_table = np.zeros((ntime, B))
    P_y_t_table = np.zeros((ntime, B))
    
    #We iterate on thetime value
    for n in range(ntime):
        t_n=n*Delta_t
        theta_1_n=theta_1_n+Omega*Delta_t
        W_x_n=a_n*V0
        
        #Initialisation of the totals
        p_x_tot=0
        p_y_tot=0
        p_n_tot=0
        p_t_tot=0
        
        #We iterate on the number of blades
        for i in range(1,B+1):
            #new definition of the theta = angle of the blades
            theta_i_n=theta_1_n+2*pi*(i-1)/B
            
            #computation of all the state
            W_x_i_n,W_y_i_n=W_i_n_calculus(W_x_n,theta_i_n,R,L)
            p_x_i_n,p_y_i_n,p_n_i_n,p_t_i_n=state_calculus(R, theta_i_n, V0, W_x_i_n, W_y_i_n, Omega, rho, c, aoa_list, cl_list, cd_list,pitch)
            
            #Sum for the globales constraints at a time
            p_x_tot+=p_x_i_n
            p_y_tot+=p_y_i_n
            p_n_tot+=p_n_i_n
            p_t_tot+=p_t_i_n
            
            if B == 3 :
                P_x_t_table[n, i-1] = p_x_i_n
                P_y_t_table[n, i-1] = p_y_i_n
            
        #Calculus of the trust coefficient and power coefficient
        CT=p_x_tot/(rho*V0**2*R)
        CP=Omega*p_t_tot/(rho*V0**3)
        
        #Computation of the a according to the BEM code
        if a_n <= 1/3:
            fg = 1
        elif a_n> 1/3:
            0.25*(5-3*a_n)
        a_qs=CT/(4*(1-fg*a_n))
        tau=2*R/V0
        a_n=a_qs+(a_n-a_qs)*exp(-Delta_t/tau)
        
        #Saving of the values of p_x and p_y total, CT and CP to plot them 
        P_x_t_list[n]=p_x_tot
        P_y_t_list[n]=p_y_tot
        t_list[n]=t_n
        CT_list[n]=CT
        CP_list[n]=CP
    return(t_list,P_x_t_list,P_y_t_list,CT_list,CP_list, P_x_t_table, P_y_t_table)

    
# Values used in the exercice
R=3
L=2.5*R
pitch=0
V0=8
Omega=14
S=0.2
rho=1.225

# Plotting and calculating values for one number of blades
def Blade_calculus(ntime,Delta_t,Omega,V0,B,S,R,L,pitch) :
    c=S*R/B
    t_list,P_x_t_list,P_y_t_list,CT_list,CP_list, P_x_t_table, P_y_t_table=unsteady_Bem (ntime,Delta_t,Omega,V0,B,R,c,L,pitch)
    
    """Power losses"""
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    plt.plot(t_list, P_x_t_list, label = 'Total Thrust ($p_{xtot}$)')
    plt.plot(t_list, P_y_t_list, label = 'Total Side Force ($p_{ytot}$)')
    plt.ylabel('N/m')
    plt.xlabel('t')
    plt.title('Total Thrust and Side Force over Time')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 2, 2)
    plt.plot(t_list, CP_list, label='Power Coefficient ($C_p$)', color='green')
    plt.xlabel('Time')
    plt.ylabel('Cp')
    plt.title('Power Coefficient over Time')
    plt.legend()
    plt.grid(True)
    
    plt.subplot(2, 2, 3)
    plt.plot(t_list, CT_list, label='Thrust Coefficient ($C_T$)', color='red')
    plt.xlabel('Time')
    plt.ylabel('CT')
    plt.title('Thrust Coefficient over Time')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    
    plt.show()
    
    # Observe the difference of phases in the case of 3 blades
    if B == 3 :
        for k in range(B) :
            plt.plot(t_list, P_x_t_table[:, k], label='B = '+str(k+1))
        plt.xlim((6,7))
        plt.legend()
        plt.grid()
        plt.show()
        
        for k in range(B) :
            plt.plot(t_list, P_y_t_table[:, k], label='B = '+str(k+1))
        plt.xlim((6,7))
        plt.legend()
        plt.grid()
        plt.show()
    
    return(t_list, P_x_t_list, P_y_t_list, CT_list, CP_list)
    

# For all Blades
def All_Blades_calculus(Omega,V0,R,L,pitch) :
    # Time step and end of time
    t_f=10
    Delta_t=0.01
    # Number of time step 
    ntime=(int(t_f/Delta_t))
    # Number of blades max
    nblades = 3
    # List of everything
    P_x_t_table = np.zeros((ntime, nblades))
    P_y_t_table = np.zeros((ntime, nblades))
    CT_table = np.zeros((ntime, nblades))
    CP_table = np.zeros((ntime, nblades))
    # Initialization
    B=0
    #this loop for is just to plot the graph for every number of blades that we have
    for k in range(nblades):
        B+=1
        t_list, P_x_t_list, P_y_t_list, CT_list, CP_list= Blade_calculus(ntime, Delta_t, Omega, V0, B, S, R, L, pitch)
        for i in range(ntime) :
            P_x_t_table[i, k] = P_x_t_list[i]
            P_y_t_table[i, k] = P_y_t_list[i]
            CT_table[i, k] = CT_list[i]
            CP_table[i, k] = CP_list[i]
        
All_Blades_calculus(Omega,V0,R,L,pitch)













