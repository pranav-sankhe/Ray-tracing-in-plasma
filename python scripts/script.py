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

freq = 10#np.power(10,6)

def sqr_ref_index(r,freq):
	val = 1 - (12400/(freq*freq))/np.power(r,6)
	return val

def singularity(r_list,ray_param):
	for r in r_list:
		val = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
		if val > 0:
			return r	

def DthetaDr(r,ray_param):
	val = -ray_param/r
	val1 = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
	val1 = np.sqrt(val1)
	val = val/val1 
	return val 

def simulation(ray_param):
	r_list = np.linspace(1,5,10000)
	rCritical = singularity(r_list,ray_param)

	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param) , rCritical, 5)[0]

	radial_dist_1 = np.linspace(rCritical,5,10000)

	sol1 = []
	for i in range(len(radial_dist_1)):
		result = integrate.quad(lambda x: DthetaDr(x,ray_param) , radial_dist_1[i], 5)
		print 'calculating for ray parameter = ',ray_param,'::', i, ' : ', result 
		sol1.append(result[0])

	sol2 = []
	for i in range(len(radial_dist_1)):
		sol2.append(2*thetaC - sol1[i] )		

	
	return [radial_dist_1,sol1,sol2]

ray_param_list = np.linspace(-4.95,4.95,20)

ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
ax.set_rmax(4)  
ax.grid(True)  

for i in range(len(ray_param_list)):
	array = simulation(ray_param_list[i])
	radial_dist = array[0]
	sol1 = array[1]
	sol2 = array[2]
	
	ax.plot(sol1,radial_dist,sol2,radial_dist ,color='b', linewidth=3)
plt.show()


