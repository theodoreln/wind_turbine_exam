""" This file is used to search for the right value of Cl and Cd by doing a double interpolation on angle of attach and thickness """
""" This is also the file where the blades profils are imported"""

# Import the numpy module
import numpy as np

# The path were the different blades profils are
path = 'C:/Users/theod/Documents/Documents importants/DTU/Courses/46300 - Wind Turbine and Aerodynamics/Assignements/Blades Profils/'

# We charge the different document
files=[path + 'FFA-W3-241.txt',
       path + 'FFA-W3-301.txt',
       path + 'FFA-W3-360.txt',
       path + 'FFA-W3-480.txt',
       path + 'FFA-W3-600.txt',
       path + 'cylinder.txt']

#Initializing tables were we are going to store our values from blades profils
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])

#Reading of tables and storing the values
for i in range(np.size(files)):
     aoa_tab[:],cl_tab[:,i],cd_tab[:,i],cm_tab[:,i] = np.loadtxt(files[i], skiprows=0).T

# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;


# Function whick return Cl, Cd and Cm for a given angle of attack and a given thickness
def force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    

    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])


    return cl, cd, cm 



""" If you want to use the function do this : """

"""
angle_of_attack=-10 # in degrees
thick = 24.1 # in percent !
[clift,cdrag,cmom]=force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab)

print('cl:',clift)
print('cd:',cdrag)
print('cm:',cmom)
"""
