# -*- coding: utf-8 -*-
""" Initialisation """
# import the classic functions (cos, sin, pi, etc...) and the interpolation function
from math import *
from Interpolation_of_coeffs import *
from Blade_geometry import *


# air density kg/m^3 for the all exercice
rho = 1.225 

""" Definition of functions useful for the BEM algorithm, just definition from the course """

# Computation of the flow angle from a, aprim, V0, omega, r
def flow_calculus(a, aprim, V0, omega, r) :
    flow = atan(((1-a)*V0)/((1+aprim)*omega*r))
    return flow

# Computation of the alpha (angle of attack) angle from the flow and the pitch
def alpha_calculus(flow, theta) :
    alpha = flow - theta
    return (alpha)

# Computation of the theta (total pitch) angle from the original pitch and the twist
def theta_calculus(pitch, twist) :
    theta = twist + pitch
    return theta

# Computation of the Cn coefficient from Cl, Cd coefficients and the flow angle
def Cn_calculus(Cl, Cd, flow) :
    Cn = Cl*cos(flow)+Cd*sin(flow)
    return Cn

# Computation of the Ct coefficient from Cl, Cd coefficients and the flow angle
def Ct_calculus(Cl, Cd, flow) :
    Ct = Cl*sin(flow)-Cd*cos(flow)
    return Ct

# Computation of the solidity from c, r and B
def solidity_calculus(c, r, B) :
    solidity = (c*B)/(2*pi*r)
    return solidity

# Computation of F coefficient from B, R, r and the flow angle
def F_calculus(B, R, r, flow) :
    F = (2/pi)*acos(exp((-B/2)*(R-r)/(r*sin(abs(flow)))))
    return F

# Computation of CT from a, Cn, solidity and the flow angle
def CT_calculu(a, Cn, solidity, flow) :
    CT = (1-a)*(1-a)*Cn*solidity/(sin(flow)*sin(flow))
    return CT

# Computation of Vrel from V0, a, aprim, omega, r
def Vrel_calculus(V0, a, aprim, omega, r) :
    Vrel = sqrt(((1-a)*V0)*((1-a)*V0)+((1+aprim)*omega*r)*((1+aprim)*omega*r))
    return Vrel

# Computation of pn from Vrel, c and Cn
def pn_calculus(Vrel, c, Cn) :
    pn = 0.5*rho*Vrel*Vrel*c*Cn
    return pn

# Computation of pt from Vrel, c and Ct
def pt_calculus(Vrel, c, Ct) :
    pt = 0.5*rho*Vrel*Vrel*c*Ct
    return pt

# Computation of the K variable for the WW algorithm
def K_calculus(F, flow, solidity, Cn) :
    K = (4*F*sin(flow)*sin(flow))/(solidity*Cn)
    return K

""" Definition of WW and Glauert algorithm to calculate a and aprim with thickness """

# Algorithm to calculate a and aprim, with glauert correction
def a_glauert_algo(V0, omega, c, R, r, B, pitch, twist, thick) :
    a0 = -1
    aprim0 = -1
    a1 = 0
    aprim1 = 0 
    # Degrees to radians
    pitch = pitch * pi / 180
    twist = twist * pi / 180 
    eps = 0.00001
    n=0
    while abs(a1-a0)>eps or abs(aprim1-aprim0)>eps :
        a0 = a1
        aprim0 = aprim1
        flow = flow_calculus(a0, aprim0, V0, omega, r)
        attack = alpha_calculus(flow, theta_calculus(pitch, twist))
        
        # If you want to use the interpolation to find Cl and Cd, outcomment below :
            
        attack_degree = attack * 180 / pi
        Cl, Cd, _ = force_coeffs_10MW(attack_degree,thick,aoa_tab,cl_tab,cd_tab,cm_tab)
        
        # If you want to use constant Cl and Cd parameters, outcomment below and change value :
            
        # Cl = 0.5
        # Cd = 0.01
        
        Cn = Cn_calculus(Cl, Cd, flow)
        Ct = Ct_calculus(Cl, Cd, flow)
        solidity = solidity_calculus(c, r, B)
        F = F_calculus(B, R, r, flow)
        if a0 < 1/3 :
            a1 = (solidity*Cn*(1-a0))/(4*F*sin(flow)*sin(flow))
        else :
            CT = CT_calculu(a0, Cn, solidity, flow)
            astar = CT/(4*F*(1-(1/4)*(5-3*a0)*a0))
            a1 = 0.1*astar+0.9*a0
        aprim1 = (solidity*Ct)*(1+aprim0)/(4*F*sin(flow)*cos(flow))
        n = n+1
    return (a1, aprim1, Cl, Cd)

# Algorithm to calculate a and aprim, with ww correction
def a_ww_algo(V0, omega, c, R, r, B, pitch, twist, thick, ac) :
    a0 = -1
    aprim0 = -1
    a1 = 0
    aprim1 = 0 
    pitch = pitch * pi / 180
    twist = twist * pi / 180 
    eps = 0.00001
    n=0
    while abs(a1-a0)>eps or abs(aprim1-aprim0)>eps :
        a0 = a1
        aprim0 = aprim1
        flow = flow_calculus(a0, aprim0, V0, omega, r)
        attack = alpha_calculus(flow, theta_calculus(pitch, twist))
        
        # If you want to use the interpolation to find Cl and Cd, outcomment below :
            
        attack_degree = attack * 180 / pi
        Cl, Cd, _ = force_coeffs_10MW(attack_degree,thick,aoa_tab,cl_tab,cd_tab,cm_tab)
        
        # If you want to use constant Cl and Cd parameters, outcomment below and change value :
            
        # Cl = 0.5
        # Cd = 0.01
        
        Cn = Cn_calculus(Cl, Cd, flow)
        Ct = Ct_calculus(Cl, Cd, flow)
        solidity = solidity_calculus(c, r, B)
        F = F_calculus(B, R, r, flow)
        if a0 < ac :
            a1 = (solidity*Cn*(1-a0))/(4*F*sin(flow)*sin(flow))
        else :
            K = K_calculus(F, flow, solidity, Cn)
            # Some problems with the ww correction at the tip
            if K < 0 :
                return (-1, -1, Cl, Cd)
            a1 = 1+K*(1-2*ac)/2-(1/2)*sqrt((K*(1-2*ac)+2)*(K*(1-2*ac)+2)+4*(K*ac*ac-1))
        aprim1 = (solidity*Ct)*(1+aprim0)/(4*F*sin(flow)*cos(flow))
        n = n+1
    return (a1, aprim1, Cl, Cd)


""" Full BEM algorithm to calculate pn and pt on one point of the blade """

# BEM algorithm, it is using Glauert correction, attention pitch and twist in degrees
def BEM_algo_glauert(V0, omega, c, R, r, B, pitch, twist, thick) :
    a, aprim, Cl, Cd = a_glauert_algo(V0, omega, c, R, r, B, pitch, twist, thick)
    flow = flow_calculus(a, aprim, V0, omega, r)
    pitch = pitch * pi / 180
    twist = twist * pi / 180 
    Cn = Cn_calculus(Cl, Cd, flow)
    Ct = Ct_calculus(Cl, Cd, flow)
    Vrel = Vrel_calculus(V0, a, aprim, omega, r)
    pn = pn_calculus(Vrel, c, Cn)
    pt = pt_calculus(Vrel, c, Ct)
    F = F_calculus(B, R, r, flow)
    return (pn, pt)

# BEM algorithm, it is using WW correction, attention pitch and twist in degrees
def BEM_algo_ww(V0, omega, c, R, r, B, pitch, twist, thick) :
    
    # 0.2 is the common value used for ac in the ww correction, so move it if it's specified to do it
    a, aprim, Cl, Cd = a_ww_algo(V0, omega, c, R, r, B, pitch, twist, thick, 0.2)
    
    # Problems with the WW correction at the tip, giving 0N loads at the end
    if a == -1 and aprim == -1 :
        return (0,0)
    
    flow = flow_calculus(a, aprim, V0, omega, r)
    pitch = pitch * pi / 180
    twist = twist * pi / 180 
    Cn = Cn_calculus(Cl, Cd, flow)
    Ct = Ct_calculus(Cl, Cd, flow)
    Vrel = Vrel_calculus(V0, a, aprim, omega, r)
    pn = pn_calculus(Vrel, c, Cn)
    pt = pt_calculus(Vrel, c, Ct)
    F = F_calculus(B, R, r, flow)
    return (pn, pt)


""" Test of the function """

# If you want to use the BEM function, it's up to you to change what you want to return in the function description
# But watch ou, the question of the assignment are meant to work with a and aprim in return !!!!!

# If you want to use the interpolation (and therefore the thickness) please see above what you need to uncomment in the a and aprim calculus functions
# But watch out again, the questions of the assignement are meant to work with the interpolation and not fix values !!!!!

# Examples with the course examples
pn, pt = BEM_algo_glauert(8, 2.61, 1.5, 31, 24.5, 3, -3, 2, 24.10)
pn, pt = BEM_algo_ww(8, 2.61, 1.5, 31, 24.5, 3, -3, 2, 24.10)

    
