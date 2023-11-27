from Assignement_1_Q1 import *
from Assignement_1_Q3 import *

R=89.17
B=3

""" Question 4 to compute all values we want """

def assignement_1_Q4(lambda_max, CP_max, CT_max, P_mech, pitch_max, r_list, c_list, twist_list, thick_list, R, B) :
    # Get the rated wind, omega_max, Power and omega from V0=0 to rated wind
    wind_rated, omega_max, P_table, wind_list = assignement_1_Q3(lambda_max, CP_max, CT_max, P_mech, R)
    # Definition of maximum, minimum and step for the V0
    wind_min=0
    wind_max=25
    wind_step=(wind_max-wind_min)/len(wind_list)
    # Retrieve the number of wind step
    i=len(wind_list)
    # Definition of maximum, minimum and step for the iteration on the pitch after rated wind speed
    p_min = pitch_max
    p_max = 30
    # Step for the pitch (here 0.1, in other case it's too long)
    p_step=0.1
    # Number of iteration by varying the pitch
    j=int((p_max-p_min)/p_step)
    # Initialisation of the final pitch list, Power list, Thrust list, CP and CT list
    pitch_list=[0]*i
    omega_list = [0]*i
    Power_list=[0]*i
    CP_list=[0]*i
    Thrust_list=[0]*i
    CT_list=[0]*i
    # We want to do this algorithm for each wind between 0 and 25
    for k in range (i):
        V0 = wind_list[k]
        # If V0 < 4, the turbine is not running 
        if V0 < 4 :
            pitch_list[k], omega_list[k], Power_list[k], CP_list[k], Thrust_list[k], CT_list[k] = pitch_max, P_table[1,k], P_table[0,k], 0, P_table[2,k], 0
        # Between 4m/s and V0_rated, we used values from the previous question
        elif V0 >= 4 and V0 < wind_rated :
            pitch_list[k], omega_list[k], Power_list[k], CP_list[k], Thrust_list[k], CT_list[k] = pitch_max, P_table[1,k], P_table[0,k], CP_max, P_table[2,k], CT_max
        # Else, we compute vary the pitch to find the place where we have the right power
        else :
            if k == 0 :
                pitch = p_min
            else :
                pitch = pitch_list[k-1]
            Pmax = 0
            for n in range (j):
                Power, Thrust, CP, CT = oneconfig_calculus_glauert(r_list, c_list, twist_list, thick_list, V0, omega_max, R, B, pitch)
                if Power > Pmax and Power < P_mech:
                    Pmax = Power
                    CPmax = CP
                    Tmax = Thrust
                    CTmax = CT
                    pitch_opt = pitch
                pitch += p_step
            pitch_list[k], omega_list[k], Power_list[k], CP_list[k], Thrust_list[k], CT_list[k] = pitch_opt, omega_max, Pmax, CPmax, Tmax, CTmax
            j = int((p_max-pitch_opt)/p_step)
        print (k*100/i,'%')
        V0 += wind_step
    return (wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list)


""" Compute the results and plot it """

wind_list,pitch_list,omega_list,Power_list,CP_list,Thrust_list,CT_list=assignement_1_Q4(lambda_max, CP_max, CT_max, P_mech, pitch_max, r_list, c_list, twist_list, thick_list, R, B)

plt.plot(wind_list, pitch_list, label = 'Pitch')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Pitch setting (deg°)')
plt.legend()
plt.grid()
plt.savefig("q4_pitch.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot(wind_list, Power_list, color='r', label = 'Power')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Power output (W)')
plt.legend()
plt.grid()
plt.savefig("q4_power.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot(wind_list, Thrust_list, color='r', label = 'Thrust')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('Thrust force (N)')
plt.legend()
plt.grid()
plt.savefig("q4_thrust.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot(wind_list, CT_list, color='g', label = 'Thrust Coefficient')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('CT')
plt.legend()
plt.grid()
plt.savefig("q4_ct.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot(wind_list, CP_list, color='g', label = 'Power Coefficient')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('CP')
plt.legend()
plt.grid()
plt.savefig("q4_cp.pdf", format="pdf", bbox_inches="tight")
plt.show()

fig, ax1 = plt.subplots()
line1, = ax1.plot(wind_list, Power_list, color='r', label = 'Power')
ax1.set_xlabel('Wind speed (m/s)')
ax1.set_ylabel('Power output (W)')
ax2 = ax1.twinx()
line2, = ax2.plot(wind_list, Thrust_list, color='b', label = 'Thrust')
ax2.set_ylabel('Thrust force (N)')
line1.set_zorder(2)
line2.set_zorder(1)
ax1.legend(loc='upper left')
ax2.legend(loc='lower right')
ax1.grid(True)
plt.savefig("q4_power_thrust.pdf", format="pdf", bbox_inches="tight")
plt.show()

plt.plot(wind_list, CP_list, color='r', label = 'Power Coefficient')
plt.plot(wind_list, CT_list, color='b', label = 'Thrust Coefficient')
plt.xlabel('Wind speed (m/s)')
plt.ylabel('CP')
plt.legend()
plt.grid()
plt.savefig("q4_cp_ct.pdf", format="pdf", bbox_inches="tight")
plt.show()

