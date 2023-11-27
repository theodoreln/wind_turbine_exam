from BEM import *
from Interpolation_of_coeffs import *
from Blade_geometry import *
from Assignement_1_Q1 import *
from Assignement_1_Q3 import *
from Assignement_1_Q4 import *


def pn_pt_calculus_glauert(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch) :
    n = 18
    pn_list = [0]*n
    pt_list = [0]*n
    for i in range(n) :
        if i == n-1 :
            pn_list[i], pt_list[i] = 0, 0
        else :
            pn_list[i], pt_list[i] = BEM_algo_glauert(V0, omega, c_list[i], R, r_list[i], B, pitch, twist_list[i], thick_list[i])
    return(pn_list, pt_list)


def assignement_1_Q5 (V0, r_list, c_list, twist_list, thick_list, wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list) :
    wind_index = wind_list.index(V0)
    pitch = pitch_list[wind_index]
    omega = omega_list[wind_index]
    Power = Power_list[wind_index]
    Thrust = Thrust_list[wind_index]
    CP = CP_list[wind_index]
    CT = CT_list[wind_index]
    pn_list, pt_list = pn_pt_calculus_glauert(r_list, c_list, twist_list, thick_list, V0, omega, R, B, pitch)
    return pitch, Power, Thrust, CP, CT, pn_list, pt_list 

""" Computation for the values we need for comparison with Ashes """ 

pitch5, Power5, Thrust5, CP5, CT5, pn_list5, pt_list5 = assignement_1_Q5 (5, r_list, c_list, twist_list, thick_list, wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list)
pitch9, Power9, Thrust9, CP9, CT9, pn_list9, pt_list9 = assignement_1_Q5 (9, r_list, c_list, twist_list, thick_list, wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list)
pitch11, Power11, Thrust11, CP11, CT11, pn_list11, pt_list11 = assignement_1_Q5 (11, r_list, c_list, twist_list, thick_list, wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list)
pitch20, Power20, Thrust20, CP20, CT20, pn_list20, pt_list520 = assignement_1_Q5 (20, r_list, c_list, twist_list, thick_list, wind_list, pitch_list, omega_list, Power_list, CP_list, Thrust_list, CT_list)

