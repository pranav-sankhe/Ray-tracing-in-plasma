import numpy as np
import simpy 
import matplotlib.pyplot as plt
import sys
from scipy.integrate import odeint
import math


	

def pend(y, theta, freq,sun_radii,mass_electron,permittivity):
	r, v = y
	val = (2*permittivity*mass_electron*math.pow(freq,2)*( 47.84*math.pow(sun_radii,16)/math.pow(r,17) + 9.3*math.pow(sun_radii,6)/math.pow(r,7) )*math.pow(10,14))/( permittivity*mass_electron*math.pow(freq,2) - (2.99*math.pow(sun_radii,16)/math.pow(r,16) + 1.55*math.pow(sun_radii,6)/math.pow(r,6) )*math.pow(10,14)) 
	dvdt = [v, (math.pow(v,2)/r) - ((math.pow(v,2) + math.pow(r,2))/v)*val]
	return dvdt

freq = 5*math.pow(10,14)
sun_radii = 696342000
mass_electron = 9.10938356*math.pow(10,-31) 
permittivity = 8.85418782*math.pow(10,-12)

Rt =  10*sun_radii
c = sun_radii*0 
var = math.pow(Rt,2) + math.pow(c,2)
ini_cond  = [math.sqrt(var) , math.sqrt(var)*(c/Rt)]

step = 0.00001*np.pi

theta  = np.linspace(0,step,2*np.pi)

sol = odeint(pend, ini_cond , theta , args=(freq,sun_radii,mass_electron,permittivity))



ax = plt.subplot(111, polar=True, axisbg='Azure')    # add subplot in polar coordinates 
ax.set_rmax(12*sun_radii)  
ax.grid(True)  
#ax.legend(loc='best')
ax.plot(theta,sol[:, 0])

plt.show()
