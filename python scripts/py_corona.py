from __future__ import division 
import numpy as np
import simpy 
import matplotlib.pyplot as plt
import sys
import math 
from sympy import *
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy import integrate
from scipy.optimize import fsolve

sun_radii = 696342000.0   
freq = 18*np.power(10,6)


def sqr_ref_index(r,freq):
	# r = r/sun_radii
	#val = 1 - (38.27/(freq*freq*np.power(r,6)))*(1 + 1.93/np.power(r,10))
	val = 1 - 38.27/np.power(r,6)
	return val

def singularity():
	no_of_samples = 0
	d = 6
	find = 4
	dist = np.linspace(0.1,4,10000)
	for distance in dist:
		if(sqr_ref_index(distance,freq) >0):  
			no_of_samples = no_of_samples + 1
			if(distance < d):
				d = distance		
	return  d 	 


def DthetaDr(r):
	val = -ray_param/r
	val1 = sqr_ref_index(r,freq)*r*r - ray_param*ray_param
	if val1 >0:
		val1 = np.sqrt(val1)
		val = val/val1 
	return val 

rC = singularity()
r = np.linspace(rC-1,5,5000)
ray_param = 2

for x in r:
	print x, sqr_ref_index(x,freq)

sol1 = []
r1 = []
sol2 = []
r2 = []
thetaC = abs(integrate.quad(lambda x: DthetaDr(x) , rC, 5)[0])
print rC ,thetaC 


for i in range(len(r)):
	result = integrate.quad(lambda x: DthetaDr(x) , r[i], 5)
	sol1.append(abs(result[0]))
	r1.append(r[i])
	
theta = np.linspace(thetaC,np.pi,10000)


for ang in theta:
	sol2.append(ang)
	r2.append( (rC*np.sin(thetaC) + np.tan(2*thetaC)*np.cos(thetaC)*rC )/(np.sin(ang) + np.cos(ang)*np.tan(2*thetaC)) )


ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
ax.plot(sol1, r1,sol2,r2, color='b', linewidth=3)
ax.set_rmax(4)  
ax.grid(True)  
#ax.legend(loc='best')

plt.show()




