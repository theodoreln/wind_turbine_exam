from Assignement_1_Q1 import *
# Value of lambda_max and pitch_max come from the file of the first question

r = 71.97

""" We want to find the Cl, Cd and angle of attack corresponding to this exact condition """

# Algorithm to retrieve Cl, Cd and the angle of attack
def Cl_Cd_alpha_q2(V0, lambda_max, c, R, r, B, pitch_max, twist, thick) :
    omega = lambda_max * V0 / R
    a, aprim, Cl, Cd = a_glauert_algo(V0, omega, c, R, r, B, pitch_max, twist, thick)
    flow = flow_calculus(a, aprim, V0, omega, r)
    # Degrees to radians
    pitch = pitch_max * pi / 180
    twist = twist * pi / 180 
    alpha = alpha_calculus(flow, theta_calculus(pitch, twist))
    return (Cl, Cd, alpha)

print('Cl, Cd and angle of attach corresponding to Q1')

Cl_max, Cd_max, alpha_max = Cl_Cd_alpha_q2(10, lambda_max, 2.91, 89.17, r, 3, pitch_max, -1.11, 24.10)
print('Cl_max = ', Cl_max)
print('Cd_max = ', Cd_max)
print('alpha_max = ', alpha_max*180/pi)


""" Now we need a modified version of the a and aprim calculus algorithm to calcul the optimal chord"""

# Flow calculus with lambda instead of omega
def flow_lambda_calculus(a, aprim, lamb, R, r) :
    flow = ((1-a)*R)/((1+aprim)*r*lamb)
    return flow

# Theta calculus with angle of attack and flow angle
def theta_alpha_calculus (flow, alpha) :
    theta = flow - alpha
    return theta

# Algorithm to compute a and aprim alpha_max already in radians, only the glauert correction for this one
def a_algo_q2(lambda_max, c, R, r, B, pitch_max, alpha_max, Cl_max, Cd_max) :
    a0 = -1
    aprim0 = -1
    a1 = 0
    aprim1 = 0 
    eps = 0.00001
    n=0
    while abs(a1-a0)>eps or abs(aprim1-aprim0)>eps :
        a0 = a1
        aprim0 = aprim1
        flow = flow_lambda_calculus(a0, aprim0, lambda_max, R, r)
        theta = theta_alpha_calculus(flow, alpha_max)
        Cn = Cn_calculus(Cl_max, Cd_max, flow)
        Ct = Ct_calculus(Cl_max, Cd_max, flow)
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
    twist = theta - pitch_max
    return (a1, aprim1, twist)

""" We need an algorithm to calculate the local dCP """

def local_CP_calculus(lambda_max, flow, c, R, B, Cl_max, Cd_max, a) :
    Ct = Ct_calculus(Cl_max, Cd_max, flow)
    CP = (B*c*lambda_max*Ct*(1-a)*(1-a))/(2*pi*R*sin(flow)*sin(flow))
    return (CP)

""" And finally, the total algorithm to iterate on c between 0 and 3m """

# Compute the best chord solution and also return CP_table and c_list
def assignement_1_Q2(lambda_max, pitch_max, alpha_max, Cl_max, Cd_max, R, r, B ) :
    # Definition of maximum, minimum and step for the chord
    c_min=0
    c_max=3
    c_step=0.01
    # Initialisation of c_list and the number of element (here 300)
    i=int((c_max-c_min)/c_step)
    c_list=[0]*i
    # Creation of the table for the local CP and the twist
    CP_table = np.zeros([2,i])
    # Initialisation of the chord for the iteration
    c = c_min
    for k in range(i):
        # We store the value of the chord c
        c_list[k] = c
        # Computation of a and aprim, the twist and then CP
        a, aprim, twist = a_algo_q2(lambda_max, c, R, r, B, pitch_max, alpha_max, Cl_max, Cd_max)
        flow = flow_lambda_calculus(a, aprim, lambda_max, R, r)
        CP = local_CP_calculus(lambda_max, flow, c, R, B, Cl_max, Cd_max, a)
        # We keep CP and the twist in the table
        CP_table[0,k] = CP
        CP_table[1,k] = twist
        c = c + c_step
    # Find the maximum value of local CP and the chord and twist corresponding to it
    CP_max = np.amax(CP_table[0,:])
    CP_argmax = np.argmax(CP_table[0,:])
    c_max = round(c_list[CP_argmax],2)
    twist_max = round(CP_table[1,CP_argmax],2)
    return (CP_table, c_list, CP_max, c_max, twist_max)

""" Results """
# CP_max = 0.5257695280083223
# c_max = 2.90
# twist_max = -0.11


""" Computation of the function """
CP_table, c_list, CP_max, c_max, twist_max = assignement_1_Q2(lambda_max, pitch_max, alpha_max, Cl_max, Cd_max, 89.17, r, 3 )
print('Results of the calculation of optimal chord')
print('CP_max = ',CP_max)
print('c_max = ',c_max)
print('twist_max = ',twist_max)

""" Plot """
# Plot general local CP
plt.plot(c_list, CP_table[0,:], label = 'CP ')
plt.xlabel('Chord (m)')
plt.ylabel('CP')
plt.legend()
plt.grid()
plt.savefig("plot_Cp_chord_gen.pdf", format="pdf", bbox_inches="tight")
plt.show()

# Plot focus local CP
plt.plot(c_list, CP_table[0,:], label = 'CP ')
plt.xlim(2.5,3)
plt.ylim(0.510,0.530)
plt.xlabel('Chord (m)')
plt.ylabel('CP')
plt.legend()
plt.grid()
plt.savefig("plot_Cp_chord_focus.pdf", format="pdf", bbox_inches="tight")
plt.show()


""" Variation where we take alpha as Cl/Cd is maximum """
# Iterate on alpha to find the best Cl and Cd
def alpha_iterate_1(thick,aoa_tab,cl_tab,cd_tab,cm_tab) :
    # Range for the iteration
    alpha_min = -30
    alpha_max= 30
    alpha_step = 0.01
    i=int((alpha_max-alpha_min)/alpha_step)
    # Creation of the table to store element
    Table = np.zeros([4,i])
    alpha = alpha_min
    for k in range(i):
        # Computation of Cl and Cd and the report
        Cl, Cd, _ = force_coeffs_10MW(alpha,thick,aoa_tab,cl_tab,cd_tab,cm_tab)
        rep = Cl/Cd
        # Store everything
        Table[0,k], Table[1,k], Table[2,k], Table[3,k] = alpha, Cl, Cd, rep
        # Iterate
        alpha = alpha + alpha_step
    # Find the maximum
    argmax = np.argmax(Table[3,:])
    alpha = Table[0,argmax]
    Cl = Table[1,argmax]
    Cd = Table[2,argmax]
    return (Cl, Cd, alpha)

print('\nVariation where we take Cl/Cd maximum')

Cl_max, Cd_max, alpha_max = alpha_iterate_1(24.10, aoa_tab, cl_tab, cd_tab, cm_tab)
print('alpha = ', alpha_max)
print('Cl_max = ', Cl_max)
print('Cd_max = ', Cd_max)

CP_table, c_list, CP_max, c_max, twist_max = assignement_1_Q2(lambda_max, pitch_max, alpha_max*pi/180, Cl_max, Cd_max, 89.17, r, 3 )
print('Results of the calculation of optimal chord')
print('CP_max = ',CP_max)
print('c_max = ',c_max)
print('twist_max = ',twist_max)



""" To study the impact of a different angle of attack on the result """

# Iterate on alpha to study the impact on optimal chord length and twist angle
def alpha_iterate_2(thick,aoa_tab,cl_tab,cd_tab,cm_tab, lambda_max, pitch_max, R, r, B ):
    # Range for the iteration
    alpha_min = 6.3
    alpha_max = 15
    alpha_step = 0.01
    i=int((alpha_max-alpha_min)/alpha_step)
    # Creation of the table to store
    alpha_list = []
    Table = np.zeros([4,i])
    # Starte
    alpha = alpha_min
    for k in range (i):
        # Computation of Cl and Cd
        Cl, Cd, _ = force_coeffs_10MW(alpha,thick,aoa_tab,cl_tab,cd_tab,cm_tab)
        Table[0,k], Table[1,k] = Cl, Cd
        # Computation of optimal c and twist
        alpha_rad = alpha * pi/180
        _, _, _, c, twist = assignement_1_Q2(lambda_max, pitch_max, alpha_rad, Cl, Cd, R, r, B )
        Table[2,k], Table[3,k] = c, twist
        alpha_list.append(alpha)
        alpha += alpha_step
    return(alpha_list, Table)

alpha_list, Table = alpha_iterate_2(24.10,aoa_tab,cl_tab,cd_tab,cm_tab, lambda_max, pitch_max, 89.17, r, 3 )

# Plot the optimal chord$
plt.plot(alpha_list, Table[2,:], label = 'Optimal Chord Length')
plt.plot(6.59, 2.90, marker="o", markersize=5, markeredgecolor="red", markerfacecolor="green")
plt.plot(10, 2.19, marker="o", markersize=5, markeredgecolor="red", markerfacecolor="green")
plt.text(6.59 + 0.5, 2.90, 'Glauert optimum', fontsize=12, verticalalignment='center')
plt.text(10 + 0.5, 2.19, 'Cl/Cd maximum', fontsize=12, verticalalignment='center')
plt.xlabel('Angle of attack (°)')
plt.ylabel('Chord Length (m)')
plt.legend()
plt.grid()
plt.savefig("plot_Opti_c_alpha.pdf", format="pdf", bbox_inches="tight")
plt.show()

# Plot the optimal twist
plt.plot(alpha_list, Table[3,:], label = 'Optimal Twist Angle')
plt.plot(6.59, -0.11, marker="o", markersize=5, markeredgecolor="red", markerfacecolor="green")
plt.plot(10, -0.17, marker="o", markersize=5, markeredgecolor="red", markerfacecolor="green")
plt.text(6.59 + 0.5, -0.11 + 0.005, 'Glauert optimum', fontsize=12, verticalalignment='center')
plt.text(10 + 0.5, -0.17 + 0.005, 'Cl/Cd maximum', fontsize=12, verticalalignment='center')
plt.xlabel('Angle of attack (°)')
plt.ylabel('Twist angle (°)')
plt.legend()
plt.grid()
plt.savefig("plot_Opti_twist_alpha.pdf", format="pdf", bbox_inches="tight")
plt.show()





