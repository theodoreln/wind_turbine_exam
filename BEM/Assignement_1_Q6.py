from math import *

A=9 #m/s
k=1.9

""" This function calculate the probability for different cases using Weibull probability density function"""
def Weibull_function (Vmin,Vmax,A,k):
    if Vmin==None and Vmax==None:
        return "error"
    else:
        if Vmin==None :
            # f(V0<=Vmax)
            f=1-exp(-(Vmax/A)**k)
        if Vmax==None :
            # f(V0>Vmin)
            f=exp(-(Vmin/A)**k)
        else:
            # f(Vmin<V0<Vmax)
            f=exp(-(Vmin/A)**k)-exp(-(Vmax/A)**k)
        return(f)


"""This function calculate the annual energy production for a specific value of V0"""    
def annual_energy_production(Power_list,wind_list,A,k,N):

    Sum_energy=0
    for i in range(0,N-2):
        Sum_energy+=1/2*(Power_list[i]+Power_list[i+1])*8760*Weibull_function(wind_list[i],wind_list[i+1],A,k)
    return Sum_energy

# You need wind_list and Power_list from q4

N=len(wind_list)
print("Energy at V0=25m/s =",annual_energy_production(Power_list,wind_list,A,k,N))


