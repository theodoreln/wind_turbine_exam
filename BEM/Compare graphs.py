import matplotlib.pyplot as plt 

#compared to all

# Power 
x1 = [5,9,11,20] 

y1 = [890.072173,5190.900913,9477.488498,10615.91345]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.') 
  
# line 2 points 
y2 = [887.446,5175.59,9449.53,10395.3]
# plotting the line 2 points  
plt.plot(x1, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('V0 wind speed [m/s]') 
# naming the y axis 
plt.ylabel('Power [kW]') 
# giving a title to my graph 
plt.title('Power to Wind Speed') 
  
# show a legend on the plot 
plt.legend() 

plt.savefig("Power.pdf", format="pdf", bbox_inches="tight")

  
# function to show the plot 
plt.show() 

# Thrust 
x1 = [5,9,11,20] 

y1 = [320.6743277,1038.984822,1552.063746,674.7767693]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.') 
  
# line 2 points 
y2 = [319.019,1033.62,1544.05,660.389]
# plotting the line 2 points  
plt.plot(x1, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('V0 wind speed [m/s]') 
# naming the y axis 
plt.ylabel('Thrust [kN]') 
# giving a title to my graph 
plt.title('Thrust to Wind Speed') 
  
# show a legend on the plot 
plt.legend() 
  
plt.savefig("Thrust.pdf", format="pdf", bbox_inches="tight")

# function to show the plot 
plt.show() 

# CP 
x1 = [5,9,11,20] 

y1 = [46.53949983,46.53949983,46.53949983,8.673090599]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.') 
  
# line 2 points 
y2 = [46.4921,46.492,46.492,8.50929]
# plotting the line 2 points  
plt.plot(x1, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('V0 wind speed [m/s]') 
# naming the y axis 
plt.ylabel('[%]') 
# giving a title to my graph 
plt.title('CP to Wind Speed') 
  
# show a legend on the plot 
plt.legend() 
  
plt.savefig("Cp.pdf", format="pdf", bbox_inches="tight")

# function to show the plot 
plt.show() 

# CT
x1 = [5,9,11,20] 

y1 = [83.8360263,83.8360263,83.8360263,11.0257117]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.') 
  
# line 2 points 
y2 = [83.5649,83.5648,83.5648,10.8115]
# plotting the line 2 points  
plt.plot(x1, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('V0 wind speed [m/s]') 
# naming the y axis 
plt.ylabel('[%]') 
# giving a title to my graph 
plt.title('CT to Wind Speed') 
  
# show a legend on the plot 
plt.legend() 

plt.savefig("Ct.pdf", format="pdf", bbox_inches="tight")

  
# function to show the plot 
plt.show() 

