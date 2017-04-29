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

Frequency = 100#np.power(10,6)              # Frequency of the incident ray.        


#------------------------------------calculate refractive index-------------------------------------------------
def sqr_ref_index(r,freq):
	val = 1 - (12400/(freq*freq))/np.power(r,6)        
	return val


#-----------------------------------calculate value of r at the singularity--------------------------------------
def singularity(r_list,ray_param,freq):
	for r in r_list:
		val = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
		if val > 0: 
			return r	

#-----------------------------------differential equation governing the trajectory of the ray-------------------
def DthetaDr(r,ray_param,freq):
	val = -ray_param/r
	val1 = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
	val1 = np.sqrt(val1)
	val = val/val1 
	return val 

#----------------------------------solve the differential equation and return the arrays of theta and r---------
def simulation(ray_param,freq):
	
	print "calculating r and theta at the singularity point"
	r_list = np.linspace(1,5,10000)
	rCritical = singularity(r_list,ray_param,freq)
	print rCritical
	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq) , rCritical, 215)[0]
	print "r at singularity:", rCritical, " theta at singularity:", thetaC

	radial_dist_1 = np.linspace(rCritical,215,100000) # define the range of values of r

	sol1 = []
	for i in range(len(radial_dist_1)):
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq) , radial_dist_1[i], 215)
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] 
		sol1.append(result[0])

	sol2 = []
	for i in range(len(radial_dist_1)):
		sol2.append(2*thetaC - sol1[i] )		

	
	return [radial_dist_1,sol1,sol2]

ray_param_list = np.linspace(-4.95,4.95,30)  #define the value of perpendicular distance between the starting point and sun's equator.              

ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
ax.set_rmax(4)  
ax.grid(True)  

for i in range(len(ray_param_list)):
	array = simulation(ray_param_list[i],Frequency)
	radial_dist = array[0]
	sol1 = array[1]
	sol2 = array[2]
	  
	ax.plot(sol1,radial_dist,sol2,radial_dist ,color='b', linewidth=1)
	sun_theta = np.linspace(0,2*np.pi,1000)
	sun_r = [1]*len(sun_theta)
	ax.set_rmax(10)
	ax.plot(sun_theta, sun_r ,color='r', linewidth=1)

plt.show()