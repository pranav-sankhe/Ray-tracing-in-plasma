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
import pyttsx
import cv2


Frequency = 100              # Frequency of the incident ray.        
Solar_radius = 695700000                   # in metres  
Te = np.power(10,6)

#------------------------------------Voice debugging-------------------------------------------------
def speak(speech):
    engine = pyttsx.init()
    engine.say(speech)
    engine.runAndWait()
	

#------------------------------------Electron density-------------------------------------------------

def eDensity(r):
	val = (1.55*np.power(10,14))/np.power(r,6)
	val = val*(1 + 1.93/np.power(r,10))
	return val 

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
def ray_trajectory(ray_param,freq):
	
	print "calculating r and theta at the singularity point"
	r_list = np.linspace(1,5,10000)
	rCritical = singularity(r_list,ray_param,freq)
	print rCritical
	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq) , rCritical, 215)[0]
	print "r at singularity:", rCritical, " theta at singularity:", thetaC

	radial_dist_1 = np.linspace(rCritical,215,10000) # define the range of values of r

	sol1 = []
	for i in range(len(radial_dist_1)):
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq) , radial_dist_1[i], 215)
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency
		sol1.append(result[0])

	sol2 = []
	for i in range(len(radial_dist_1)):
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq), rCritical ,radial_dist_1[i])
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency 
		sol2.append(result[0] + thetaC)
	
	# sol3 = []
	# for i in range(len(radial_dist_1)):

	# 	sol3.append(2*thetaC - sol1[i] )		

	return [radial_dist_1,sol1,sol2]

#-----------------------------------calculate integral for optical depth-----------------------------

def lineIntegral(r,ray_param,freq):
	val = np.square(eDensity(r))*r
	val1 = sqr_ref_index(r,freq)*r*r - ray_param*ray_param	
	val1 = np.sqrt(val1)
	val = val/val1
	return val 	

#-----------------------------------calculate optical depth------------------------------------------
def opticalDepth(r1,r2,freq,ray_param):
	integral = integrate.quad(lambda x: lineIntegral(x,ray_param,freq*np.power(10,6)), r1 ,r2)[0]
	tau = (Solar_radius*integral)/(np.square(freq*np.power(10,6))*np.power(Te,1.5)*np.power(10,11)) 
	return tau	

#-----------------------------------calculate electron Temperature------------------------------------------
def electron_Temp(r):
	if r == 1:
		#speak("on surface")
		print "ON SURFACE"
		temp = 6200
	if r > 1 and r < 1.0005:
		#speak("photosphere")
		print "PHOTOSPHERE"
		temp = 6200 - np.power(5,6)*(r-1)
	if r >1.0005 and r < 1.00301:
		#speak("Chromosphere")
		print "CHROMOSPHERE"
		temp = 3700 + 1593625.498*(r-1.0005)
	if r >1.00301 and r < 1.00316:
		#speak("transition region")
		print "TRANSITION REGION"
		temp = 7700 + 6615333333*(r-1.00301)
	if r > 1.00316:
		temp = np.power(10,6)
	return temp 		  	 	 	

#-----------------------------------calculate brightness temperature------------------------------------------
def brightness_temp(freq,ray_param):
	print "calculating r and theta at the singularity point"
	r_list = np.linspace(1,5,10000)
	rCritical = singularity(r_list,ray_param,freq)
	print rCritical
	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq) , rCritical, 215)[0]
	print "r at singularity:", rCritical, " theta at singularity:", thetaC

	radial_dist1 = np.linspace(rCritical,5,10000) # define the range of values of r
	radial_dist2 = np.linspace(5,rCritical,10000) # define the range of values of r

	observed_temp = 0 	
	background_Temp = 0

	for i in range(len(radial_dist2)-1):	 
		tau = opticalDepth(radial_dist2[i+1],radial_dist2[i],freq,ray_param)
		observed_temp = background_Temp*np.exp(-tau) + electron_Temp(radial_dist2[i+1])*(1-np.exp(-tau)) 
		print i, ": running from 5 to Rc. Currently r1 = ",radial_dist2[i+1],"r2=",radial_dist2[i],'ray parameter=',ray_param, "tau=",tau,"temp =", observed_temp 
		#print background_Temp,electron_Temp(radial_dist2[i+1]),np.exp(-tau),electron_Temp(radial_dist2[i+1])*(1-np.exp(-tau)) 
		background_Temp = observed_temp	 
	 
	for i in range(len(radial_dist1)-1):
		print i, ": running from Rc to 5. Currently r1 = ",radial_dist1[i],"r2=",radial_dist1[i+1],'ray parameter = ', ray_param, "tau=",tau, "temp =", observed_temp
		tau = opticalDepth(radial_dist1[i],radial_dist1[i+1],freq,ray_param)
		observed_temp = background_Temp*np.exp(-tau) + electron_Temp(radial_dist1[i+1])*(1-np.exp(-tau)) 
		background_Temp = observed_temp	 

	print "Brightness temperature" ,observed_temp
	return observed_temp
	

#=================================================================================================================
def plot_trajectory(freq,param_list):
	
	ray_param_list = param_list
	ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
	ax.set_rmax(4)  
	ax.grid(True)  

	for i in range(len(ray_param_list)):
		array = ray_trajectory(ray_param_list[i],freq)
		radial_dist = [value for value in array[0] if value < 10]
		sol1 = array[1][0:len(radial_dist)]
		sol2 = array[2][0:len(radial_dist)]
		ax.plot(sol1,radial_dist,sol2,radial_dist ,color='b', linewidth=1)
		sun_theta = np.linspace(0,2*np.pi,1000)
		sun_r = [1]*len(sun_theta)
		ax.plot(sun_theta, sun_r ,color='r', linewidth=1)

	plt.show()


def tempProfile(freq,param_list):
	ray_param_list = param_list
	observed_temp = []
	for a in ray_param_list:
		observed_temp.append(brightness_temp(freq,a))
	plt.plot(param_list, observed_temp)
	plt.show()
	return observed_temp	
	# img = np.zeros((512,512,3), np.uint8)
	# for radius in radii
	# 	cv2.circle(img,(447,63), 63, (0,0,255), -1)

ray_param_list = np.linspace(0,4.95,100)  #define the value of perpendicular distance between the starting point and sun's equator           
#ray_param_list = [0]
tempProfile(Frequency,ray_param_list)
#plot_trajectory(Frequency,ray_param_list)