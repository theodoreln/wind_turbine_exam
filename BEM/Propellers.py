""" Initialisation """
# import the classic functions (cos, sin, pi, etc...) and the interpolation function
from math import *
import matplotlib.pyplot as plt

# air density kg/m^3 for the all exercice
rho = 1.225 

""" Table of the exercice """

r_list = [0.2, 0.3, 0.4, 0.5, 0.6, 0.645, 0.690]
twist_list = [23, 20, 16, 14.5, 13, 12.3, 11.7]
c_list = [0.106, 0.117, 0.112, 0.103, 0.088, 0.082, 0]

""" Definition of functions useful for the BEM algorithm, just definition from the course """

# Computation of the flow angle from a, aprim, V0, omega, r
def flow_calculus(a, aprim, V0, omega, r) :
    flow = atan(((1+a)*V0)/((1-aprim)*omega*r))
    return flow

# Computation of the alpha (angle of attack) angle from the flow and the pitch
def alpha_calculus(flow, theta) :
    alpha = theta - flow
    return (alpha)

# Computation of the Cn coefficient from Cl, Cd coefficients and the flow angle
def Cn_calculus(Cl, Cd, flow) :
    Cn = Cl*cos(flow)-Cd*sin(flow)
    return Cn

# Computation of the Ct coefficient from Cl, Cd coefficients and the flow angle
def Ct_calculus(Cl, Cd, flow) :
    Ct = Cl*sin(flow)+Cd*cos(flow)
    return Ct

# Computation of the solidity from c, r and B
def solidity_calculus(c, r, B) :
    solidity = (c*B)/(2*pi*r)
    return solidity

# Computation of F coefficient from B, R, r and the flow angle
def F_calculus(B, R, r, flow) :
    F = (2/pi)*acos(exp((-B/2)*(R-r)/(r*sin(abs(flow)))))
    return F

# Computation of Vrel from V0, a, aprim, omega, r
def Vrel_calculus(V0, a, aprim, omega, r) :
    Vrel = sqrt(((1+a)*V0)*((1+a)*V0)+((1-aprim)*omega*r)*((1-aprim)*omega*r))
    return Vrel

# Computation of pn from Vrel, c and Cn
def pn_calculus(Vrel, c, Cn) :
    pn = 0.5*rho*Vrel*Vrel*c*Cn
    return pn

# Computation of pt from Vrel, c and Ct
def pt_calculus(Vrel, c, Ct) :
    pt = 0.5*rho*Vrel*Vrel*c*Ct
    return pt

""" Definition of BEM algorithm to calculate a, aprim, pn and pt """

# Algorithm to calculate a and aprim, with glauert correction
def BEM_propeller(V0, omega, c, R, r, B, twist) :
    a0 = -1
    aprim0 = -1
    a1 = 0
    aprim1 = 0 
    # Degrees to radians
    twist = twist * pi / 180 
    eps = 0.00001
    n=0
    while abs(a1-a0)>eps or abs(aprim1-aprim0)>eps :
        a0 = a1
        aprim0 = aprim1
        flow = flow_calculus(a0, aprim0, V0, omega, r)
        alpha = alpha_calculus(flow, twist)
        attack_degree = alpha * 180 / pi
        Cl, Cd = 0.1*attack_degree+0.4, 0.008
        Cn = Cn_calculus(Cl, Cd, flow)
        Ct = Ct_calculus(Cl, Cd, flow)
        solidity = solidity_calculus(c, r, B)
        F = F_calculus(B, R, r, flow)
        astar = (1+a0)*(solidity*Cn)/(4*F*sin(flow)*sin(flow))
        a1 = 0.1*astar+0.9*a0
        aprimstar = (1-aprim0)*(solidity*Ct)/(4*F*sin(flow)*cos(flow))
        aprim1 = 0.1*aprimstar+0.9*aprim0
        n = n+1
    a = a1
    aprim = aprim1
    Vrel = Vrel_calculus(V0, a, aprim, omega, r)
    pn = pn_calculus(Vrel, c, Cn)
    pt = pt_calculus(Vrel, c, Ct)
    return (a, aprim, pn, pt)


""" Iteration on the all turbine for one given wind speed """

# Calculation of power and thrust with the list of loads and distance from the rotor
def Power_Thrust_calculus(r_list, pn_list, pt_list, omega, B) :
    Power = 0
    Thrust = 0
    n = len(r_list)
    for i in range(n-1) :
        An = (pn_list[i+1]-pn_list[i])/(r_list[i+1]-r_list[i])
        Bn = (pn_list[i]*r_list[i+1]-pn_list[i+1]*r_list[i])/(r_list[i+1]-r_list[i])
        Mn = (An/2)*(r_list[i+1]**2-r_list[i]**2)+(Bn)*(r_list[i+1]-r_list[i])
        At = (pt_list[i+1]-pt_list[i])/(r_list[i+1]-r_list[i])
        Bt = (pt_list[i]*r_list[i+1]-pt_list[i+1]*r_list[i])/(r_list[i+1]-r_list[i])
        Mt = (At/3)*(r_list[i+1]**3-r_list[i]**3)+(Bt/2)*(r_list[i+1]**2-r_list[i]**2)
        Thrust = Thrust + Mn
        Power = Power + Mt
    Thrust = Thrust * B
    Power = Power * B * omega
    return (Power, Thrust)

# Calculation of CP with the Power and the area of the blades
def CP_calculus(Power, V0, A) :
    CP = Power/(0.5*rho*V0*V0*V0*A)
    return CP

# Calculation of CT with the Thrust and the area of the blades
def CT_calculus(Thrust, V0, A) :
    CT = Thrust/(0.5*rho*V0*V0*A)
    return CT

# Computation of Power, Thrust, CP and CT on the turbine for a given V0
def one_speed_config(r_list, c_list, twist_list, V0, omega, R, B) :
    n = len(r_list)
    a_list = [0]*n
    aprim_list = [0]*n
    pn_list = [0]*n
    pt_list = [0]*n
    for i in range(n) :
        # No loads at the end
        if i == n-1 :
            pn_list[i], pt_list[i] = 0, 0
        else :
            # We call the BEM algo at every spot of the blade
            a_list[i], aprim_list[i], pn_list[i], pt_list[i] = BEM_propeller(V0, omega, c_list[i], R, r_list[i], B, twist_list[i])
    Power, Thrust = Power_Thrust_calculus(r_list, pn_list, pt_list, omega, B)
    CP = CP_calculus(Power, V0, pi*R*R)
    CT = CT_calculus(Thrust, V0, pi*R*R)
    return(Power, Thrust, CP, CT, a_list, aprim_list, pn_list, pt_list)


""" Iteration on multiple wind speed """

# Function to iterate on V0
def Wind_iteration(r_list, c_list, twist_list, omega, R, B) :
    V_min = 1
    V_max = 40
    V_step = 1
    n = int((V_max-V_min+1)/V_step)
    V_list = [0]*n
    P_list = [0]*n
    T_list = [0]*n
    CP_list = [0]*n
    CT_list = [0]*n
    Eff_list = [0]*n
    for i in range(n):
        V_list[i] = V_min+i*V_step
        P_list[i], T_list[i], CP_list[i], CT_list[i], _, _, _, _ = one_speed_config(r_list, c_list, twist_list, V_list[i], omega, R, B)
        Eff_list[i] = (T_list[i]*V_list[i])/P_list[i]
    return (V_list, P_list, T_list, CP_list, CT_list, Eff_list)
    

""" Test of the functions """

# Examples with the course examples
Power, Thrust, CP, CT, a_list, aprim_list, pn_list, pt_list = one_speed_config(r_list, c_list, twist_list, 10, 2600*2*pi/60, 0.690, 2)
plt.plot(r_list, pn_list, label='pn')
plt.plot(r_list, pt_list, label='pt')
plt.legend()
plt.grid()
plt.show()

V_list, P_list, T_list, CP_list, CT_list, Eff_list = Wind_iteration(r_list, c_list, twist_list, 2600*2*pi/60, 0.690, 2)
plt.plot(V_list, P_list, label='P')
plt.plot(V_list, T_list, label='T')
plt.legend()
plt.grid()
plt.show()

plt.plot(V_list, Eff_list, label='Efficiency')
plt.legend()
plt.grid()
plt.show()

Drag_list = [0]*40
for i in range(40):
    Drag_list[i] = 2.44*(i+1)*(i+1)
plt.plot(V_list, Drag_list, '--', label='Drag')
plt.plot(V_list, T_list, label='T')
plt.legend()
plt.grid()
plt.show()

    
