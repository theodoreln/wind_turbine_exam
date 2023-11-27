import matplotlib.pyplot as plt
from Interpolation_of_coeffs import *
from Blade_geometry import *
from BEM import *


""" Structure of the algorithm for q1 """
# Objective : Find the maximum Cp by varrying the tip speed ratio and the pitch angle

# Informations : For every r we have c, twist, t/c (to find Cl and Cd frome alpha)
# + R, B 

# Given the BEM algorithm : BEM_algo(V0, omega, c, R, r, B, pitch, twist, Cl, Cd) we miss :
    # V0 and omega come from the tip speed ratio, we choose one V0 and make the omega vary
    # pitch will be one of the loop
    
# Choose V0 = 10 m/s for example (only the tip speed ratio matters)
# For every (tip speed ratio) lambda from 5 to 10 (step 0.1 ?) :
    # omega = lambda*V0/R
    # For every pitch from -4 to 3 (degrees ?) (step 0.1 ?) :
        # For every element of the blades (table with r) :
             # Compute the BEM algorithm
             # !!!! Change the algorithm which do not take Cl and Cd as entrance but t/c
             # !!!! Double interpolation with the angle of attack and t/c
             # !! Watch out for the load at the end
        # Compute power and thrust from load distributions 
        # Compute Cp(lambda, pitch) and store it somewhere
# Plot and determine Cp,max(lambda, pitch)


""" Computation of the power and thrust of the wind turbine """
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
        

""" Algorithm for the calculus of one configuration """

# Function which gives Power, Thrust, CP, pn and pt for one configuration of omega and pitch, using glauert correction
def oneconfig_calculus_glauert(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch) :
    n = 18
    pn_list = [0]*n
    pt_list = [0]*n
    for i in range(n) :
        # No loads at the end
        if i == n-1 :
            pn_list[i], pt_list[i] = 0, 0
        else :
            # We call the BEM algo at every spot of the blade
            pn_list[i], pt_list[i] = BEM_algo_glauert(V0, omega, c_list[i], R, r_list[i], B, pitch, twist_list[i], thick_list[i])
    Power, Thrust = Power_Thrust_calculus(r_list, pn_list, pt_list, omega, B)
    CP = CP_calculus(Power, V0, pi*R*R)
    CT = CT_calculus(Thrust, V0, pi*R*R)
    return(Power, Thrust, CP, CT)

# Function which gives Power, Thrust, CP, pn and pt for one configuration of omega and pitch, using ww correction
def oneconfig_calculus_ww(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch) :
    n = 18
    pn_list = [0]*n
    pt_list = [0]*n
    for i in range(n) :
        # We call the BEM algo at every spot of the blade
        pn_list[i], pt_list[i] = BEM_algo_ww(V0, omega, c_list[i], R, r_list[i], B, pitch, twist_list[i], thick_list[i])
    Power, Thrust = Power_Thrust_calculus(r_list, pn_list, pt_list, omega, B)
    CP = CP_calculus(Power, V0, pi*R*R)
    CT = CT_calculus(Thrust, V0, pi*R*R)
    return(Power, Thrust, CP, CT)

# Test 
Power, Thrust, CP, CT = oneconfig_calculus_glauert(r_list, c_list, twist_list, thick_list, 10, 8.029*2*pi/60, 89.17, 3, 0)
Power, Thrust, CP, CT = oneconfig_calculus_ww(r_list, c_list, twist_list, thick_list, 10, 8.029*2*pi/60, 89.17, 3, 0)


""" Algorithm for the iteration on lambda and the pitch """

# With Glauert correction
def assignement_1_Q1_glauert(V0,R,B,r_list,c_list,twist_list,thick_list):
    # Definition of maximum, minimum and step for lambda
    lambda_min=5
    lambda_max=10
    l_step=0.1
    # Initialisation of lambda_list and the number of element
    i=int((lambda_max-lambda_min)/l_step)
    lambda_list=[0]*i
    # Definition of maximum, minimum and step for pitch
    pitch_min=-4
    pitch_max=3
    p_step=0.1
    # Initialisation of pitch list and the numer of element
    j=int((pitch_max-pitch_min)/p_step)
    pitch_list=[0]*j
    # Creation of the different table for CP, the Power and the Thrust
    CP_table = np.zeros([i,j])
    CT_table = np.zeros([i,j])
    P_table = np.zeros([i,j])
    T_table = np.zeros([i,j])
    
    # Initialisation of lambda for the iteration
    lamb = lambda_min
    for k in range(i):
        # Calculus of omega based on the value of lambda and adding the value of lambda in the table
        omega = lamb*V0/R
        lambda_list[k] = round(lamb,1)
        
        # Initialisation of pitch for the iteration
        pitch = pitch_min       
        for l in range(j):
            # We just want to store the pitch one time and not at each iteration
            if k==0 :
                pitch_list[l] = round(pitch,1)
            # Computation of the Power, the Thrust and CP for a given configuration
            Power,Thrust,CP,CT = oneconfig_calculus_glauert(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch)
            CP_table[k,l] = CP 
            CT_table[k,l] = CT
            P_table[k,l] = Power 
            T_table[k,l] = Thrust
            # Next iteration on the pitch
            pitch = pitch + p_step
        # Next iteration on lambda
        lamb = lamb + l_step
        print (k*100/i,'%')
    # Find the maximum value and the lambda and pitch coefficients corresponding to it
    CP_max = np.amax(CP_table)
    CP_argmax = np.argmax(CP_table)
    nb_row = floor(CP_argmax/j)
    lambda_max = round(lambda_list[nb_row],2)
    pitch_max = round(pitch_list[CP_argmax-nb_row*j],2)
    
    # Find the value of CT, Power and Thrust
    CT_max = CT_table[nb_row, CP_argmax-nb_row*j]
    Power_max = P_table[nb_row, CP_argmax-nb_row*j]
    Thrust_max = T_table[nb_row, CP_argmax-nb_row*j]
    return CP_table, CT_table, lambda_list, pitch_list, CP_max, lambda_max, pitch_max, CT_max, Power_max, Thrust_max

# With WW correction 
def assignement_1_Q1_ww(V0,R,B,r_list,c_list,twist_list,thick_list):
    # Definition of maximum, minimum and step for lambda
    lambda_min=5
    lambda_max=10
    l_step=0.1
    # Initialisation of lambda_list and the number of element (here 50)
    i=int((lambda_max-lambda_min)/l_step)
    lambda_list=[0]*i
    # Definition of maximum, minimum and step for pitch
    pitch_min=-4
    pitch_max=3
    p_step=0.1
    # Initialisation of pitch list and the numer of element (here 70)
    j=int((pitch_max-pitch_min)/p_step)
    pitch_list=[0]*j
    # Creation of the different table for CP, the Power and the Thrust
    CP_table = np.zeros([i,j])
    CT_table = np.zeros([i,j])
    P_table = np.zeros([i,j])
    T_table = np.zeros([i,j])
    
    # Initialisation of lambda for the iteration
    lamb = lambda_min
    for k in range(i):
        # Calculus of omega based on the value of lambda and adding the value of lambda in the table
        omega=lamb*V0/R
        lambda_list[k]= round(lamb,1)
        
        # Initialisation of pitch for the iteration
        pitch = pitch_min       
        for l in range(j):
            # We just want to store the pitch one time and not at each iteration
            if k==0 :
                pitch_list[l]= round(pitch,1)
            # Computation of the Power, the Thrust and CP for a given configuration
            Power,Thrust,CP,CT = oneconfig_calculus_ww(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch)
            CP_table[k,l] = CP 
            CT_table[k,l] = CT
            P_table[k,l] = Power 
            T_table[k,l] = Thrust
            # Next iteration on the pitch
            pitch = pitch + p_step
        # Next iteration on lambda
        lamb = lamb + l_step
        print (k*100/i,'%')
    # Find the maximum value and the lambda and pitch coefficients corresponding to it
    CP_max = np.amax(CP_table)
    CP_argmax = np.argmax(CP_table)
    nb_row = floor(CP_argmax/j)
    lambda_max = round(lambda_list[nb_row],2)
    pitch_max = round(pitch_list[CP_argmax-nb_row*j],2)
    
    # Find the value of CT, Power and Thrust
    CT_max = CT_table[nb_row, CP_argmax-nb_row*j]
    Power_max = P_table[nb_row, CP_argmax-nb_row*j]
    Thrust_max = T_table[nb_row, CP_argmax-nb_row*j]
    return CP_table, CT_table, lambda_list, pitch_list, CP_max, lambda_max, pitch_max, CT_max, Power_max, Thrust_max

# CP_table, CT_table -> table with every value of CP and CT, row -> lambda, column -> pitch
# lambda_list, pitch_list -> list with every value of lambda and pitch used, correspond to row and column in table
# CP_max, lambda_max, pitch_max, CT_max, Power_max, Thrust_max -> all values at the maximum CP point


""" Results of the computation for the assignement and step = 0.01 """
# Glauert (step = 0.01)
lambda_max = 8.01
pitch_max = 0.09
CP_max = 0.46539499834694914
CT_max = 0.8383602630214817
# Power_max (10 m/s) = 7120577.384235464
# Thrust_max (10 m/s) = 1282697.310868423

# WW (step = 0.01)
# lambda_max = 8.71
# pitch_max = -2.32
# CP_max = 0.5255144563745352
# CT_max = 1.0950742198824845
# Power_max (10 m/s) = 8040409.472470731
# Thrust_max (10 m/s) = 1675471.535330399


""" Please outcomment all the rest when lauching next questions """
""" Computation of the function for the exam """
cp_table_glauert, ct_table_glauert, lambda_list_glauert, p_list_glauert, cp_max_glauert, lambda_max_glauert, pitch_max_glauert, ct_max_glauert, power_max_glauert, thrust_max_glauert = assignement_1_Q1_glauert(10, 89.17, 3, r_list, c_list, twist_list, thick_list)
cp_table_ww, ct_table_ww, lambda_list_ww, p_list_ww, cp_max_ww, lambda_max_ww, pitch_max_ww, ct_max_ww, power_max_ww, thrust_max_ww = assignement_1_Q1_ww(10, 89.17, 3, r_list, c_list, twist_list, thick_list)
print('glauert correction gives us :')
print('cp_max = ',cp_max_glauert)
print('lambda_max = ',lambda_max_glauert)
print('pitch_max = ',pitch_max_glauert)
print('ww correction gives us :')
print('cp_max = ',cp_max_ww)
print('lambda_max = ',lambda_max_ww)
print('pitch_max = ',pitch_max_ww)


""" Plot CP and CT for a specifid pitch, depending on lambda """

# Search the index of a specific pitch
pitch_searched = 0
pitch_index = p_list_glauert.index(pitch_searched)

# Plot multiple CP, depending on the value of the pitch, with lambda on the axis
plt.plot(lambda_list_glauert, cp_table_glauert[:,pitch_index], label = 'pitch = {0}°'.format(pitch_searched))
plt.ylabel('CP')
plt.xlabel('lambda')
plt.legend()
plt.grid()
plt.savefig("plot_CP_lambda_pitch.pdf", format="pdf", bbox_inches="tight")
plt.show() 

# Plot multiple CT, depending on the value of the pitch, with lambda on the axis
plt.plot(lambda_list_glauert, ct_table_glauert[:,pitch_index], label = 'pitch = {0}°'.format(pitch_searched))
plt.ylabel('CT')
plt.xlabel('lambda')
plt.legend()
plt.grid()
plt.savefig("plot_CT_lambda_pitch.pdf", format="pdf", bbox_inches="tight")
plt.show() 


""" Searched for specific values of CP and CT, with lambda and pitch """

# Search the index of a specific pitch
pitch_searched = 0
pitch_index = p_list_glauert.index(pitch_searched)
# Search the index of a specific lambda
lambda_searched = 7.5
lambda_index = lambda_list_glauert.index(lambda_searched)

CP_searched = cp_table_glauert[lambda_index, pitch_index]
CT_searched = ct_table_glauert[lambda_index, pitch_index]


""" Compute the loads on one example for the report """

# The only difference between this one and the previous function is that this one is returning also the list of loads on the blade
def oneconfig_calculus_glauert_ex(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch) :
    n = 18
    pn_list = [0]*n
    pt_list = [0]*n
    for i in range(n) :
        # No loads at the end
        if i == n-1 :
            pn_list[i], pt_list[i] = 0, 0
        else :
            pn_list[i], pt_list[i] = BEM_algo_glauert(V0, omega, c_list[i], R, r_list[i], B, pitch, twist_list[i], thick_list[i])
    Power, Thrust = Power_Thrust_calculus(r_list, pn_list, pt_list, omega, B)
    CP = CP_calculus(Power, V0, pi*R*R)
    CT = CT_calculus(Thrust, V0, pi*R*R)
    return(Power, Thrust, CP, CT, pn_list , pt_list)

P_ex, T_ex, CP_ex, CT_ex, pn_list_ex, pt_list_ex = oneconfig_calculus_glauert_ex(r_list, c_list, twist_list, thick_list, 10, 7*10/89.17, 89.17, 3, -2)

# Plot the loads distribution 
plt.plot(r_list, pn_list_ex, 'x-', label = 'pn')
plt.plot(r_list, pt_list_ex, 'x-', label = 'pt')
plt.ylabel('Loads (N/m)')
plt.xlabel('Radial position (m)')
plt.legend()
plt.grid()
plt.savefig("plot_loads_ex.pdf", format="pdf", bbox_inches="tight")
plt.show() 

