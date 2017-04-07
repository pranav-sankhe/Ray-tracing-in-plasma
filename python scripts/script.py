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

freq = 18#np.power(10,6)

def sqr_ref_index(r,freq):
	val = 1 - 38.27/np.power(r,6)
	return val

def singularity(r_list):
	for r in r_list:
		val = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
		if val > 0:
			return r	

def DthetaDr(r):
	val = -ray_param/r
	val1 = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
	val1 = np.sqrt(val1)
	val = val/val1 
	return val 

ray_param = 3

r_list = np.linspace(1,5,10000)
rCritical = singularity(r_list)
thetaC = integrate.quad(lambda x: DthetaDr(x) , rCritical, 5)[0]

radial_dist = np.linspace(rCritical,5,10000)

sol1 = []
for i in range(len(radial_dist)):
	result = integrate.quad(lambda x: DthetaDr(x) , radial_dist[i], 5)
	print i, ' : ', result
	sol1.append(result[0])

sol2 = []

ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
ax.plot(sol1, radial_dist ,color='b', linewidth=3)
ax.set_rmax(4)  
ax.grid(True)  
#ax.legend(loc='best')

plt.show()
