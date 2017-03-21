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
freq = 18#np.power(10,6)
sintheta = 0

def Sqr_refIndex(r):
	return  1 - 12400/(np.square(freq)*np.power(r,6))

def der_refIndex(r):
	return 6*12400/(np.square(freq)*np.power(r,7))

def singularity():
	d = 0
	dist = np.linspace(0.1,5,10000)
	count = 0 
	for distance in dist:
		if(Sqr_refIndex(distance) > 0 and count == 0):
			d = distance
			count = 1		
	
	return  d 	 
	

def dydt(y,theta,b,c):
	r, der_r = y 
	y_der = [der_r, (74400*(1+r*r))/(r*der_r*(freq*freq*np.power(r,6)-12400)) + der_r/r ]
	
	return y_der

b=0
c=0 

y0 = [ 5*np.sin(sintheta), 0.1 ]

theta = np.linspace(0.001,np.pi,1000)

sol = odeint(dydt, y0, theta, args=(b, c)) 

var = []
# for x in sol:
# 	print x
dist = np.linspace(singularity(),5,10000)
for x in dist:
	var.append(der_refIndex(x))


plt.plot(dist, var)
plt.xlabel('x')
plt.ylabel('y')
plt.show()
	

