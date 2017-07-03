from __future__ import division 
import numpy as np
#import simpy 
import matplotlib.pyplot as plt
import sys
import math 
from sympy import *
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy import integrate
#import pyttsx
#import cv2

Frequency = 10*np.power(10,6)            # Frequency of the incident ray.        
Solar_radius = 695700000                   # in metres  
Te = np.power(10,6)

#------------------------------------Voice debugging-------------------------------------------------
# def speak(speech):
#     engine = pyttsx.init()
#     engine.say(speech)
#     engine.runAndWait()
	
#------------------------------------Electron density-------------------------------------------------

# def eDensity(r):
	
# 	if flag == 0:  
# 		val = (1.55*np.power(10,14))/np.power(r,6)
# 		val = val*(1 + 1.93/np.power(r,10))
# 		return val 

# 	if flag == 1:
			
def eDensity(r,flag,cone_angle,angle):
	if flag == 0:  
		val = (1.55*np.power(10,14))/np.power(r,6)
		val = val*(1 + 1.93/np.power(r,10))
		return val 
	if flag == 1:
		val = (1.55*np.power(10,14))/np.power(r,6)
		val = val*(1 + 1.93/np.power(r,10))
		val = val*np.sin(abs(angle))/np.sin(abs(cone_angle/2))
		return val

#------------------------------------calculate refractive index-------------------------------------------------
def sqr_ref_index(r,freq,flag,cone_angle,angle):
	val = 1 - eDensity(r,flag,cone_angle,angle)*179.64/(freq*freq)
	return val

#-----------------------------------calculate value of r at the singularity--------------------------------------
def singularity(r_list,ray_param,r,freq,flag,cone_angle,angle):
	for r in r_list:
		val = sqr_ref_index(r,freq,flag,cone_angle,angle)*r*r - ray_param*ray_param
		if val > 0:
			return r

#-----------------------------------differential equation governing the trajectory of the ray-------------------
def DthetaDr(r,ray_param,freq,flag,cone_angle,angle):
	val = -ray_param/r
	val1 = sqr_ref_index(r,freq,flag,cone_angle,angle)*r*r - ray_param*ray_param
	val1 = np.sqrt(abs(val1))
	val = val/val1 
	#print "final",val
	return val 

#----------------------------------solve the differential equation and return the arrays of theta and r---------
def ray_trajectory(ray_param,freq,cone_angle):
	
	print "calculating r and theta at the singularity point"
	r_list = np.linspace(1,5,10000)
	angle = 0
	flag = 0
	
	rCritical = singularity(r_list,ray_param,0,freq,flag,cone_angle,angle)
	print "rCritical", rCritical
	
	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , rCritical, 215)[0]
	print "Theta critical", thetaC	
	if abs(thetaC) > cone_angle/2:
		pass
	if abs(thetaC) < cone_angle/2:
		flag = 1
		rCritical = singularity(r_list,ray_param,rCritical,freq,flag,cone_angle,angle)
		thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , rCritical, 215)[0]

	print "r at singularity:", rCritical, " theta at singularity:", thetaC

	radial_dist_1 = np.linspace(rCritical,215,10000) # define the range of values of r

	sol1 = []
	for i in range(len(radial_dist_1)):
		flag = 0
		angle = 0
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist_1[i], 215)
		
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency
		angle = result[0]

		if abs(angle) > cone_angle/2:
			pass
		if abs(angle) < cone_angle/2:
			flag = 1		   
			print angle
			result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist_1[i], 215)
			print 'calculating for ray parameter with flag 1 = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency			

		sol1.append(result[0])
	

	sol2 = []
	for i in range(len(radial_dist_1)):
		
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle), rCritical ,radial_dist_1[i])
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency 
		
		angle = result[0]

		if abs(angle) > cone_angle/2:
			pass
		if abs(angle) < cone_angle/2:
			flag = 1		   
			result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist_1[i], 215)
			print 'calculating for ray parameter with flag 1 = ',ray_param,'::', i,':' ,'r:',radial_dist_1[i] , 'theta :',result[0] , 'at frequency', Frequency	

		sol2.append(result[0] + thetaC)		

	return [radial_dist_1,sol1,sol2]

#-----------------------------------calculate integral for optical depth-----------------------------
def lineIntegral(r,ray_param,freq,flag,cone_angle,angle):
	val = np.square(eDensity(r,flag,cone_angle,angle))*r
	val1 = sqr_ref_index(r,freq,flag,cone_angle,angle)*r*r - ray_param*ray_param	
	val1 = np.sqrt(abs(val1))
	val = val/val1
	return val 	

#-----------------------------------calculate optical depth------------------------------------------
def opticalDepth(r1,r2,freq,ray_param,flag,cone_angle,angle):
	integral = integrate.quad(lambda x: lineIntegral(x,ray_param,freq, flag,cone_angle,angle), r1 ,r2,limit=50000)
	integral = integral[0]
	tau = (Solar_radius*integral)/(np.square(freq)*np.power(Te,1.5)*np.power(10,11))
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
#		speak("transition region")
		print "TRANSITION REGION"
		temp = 7700 + 6615333333*(r-1.00301)
	if r > 1.00316:
		temp = np.power(10,6)
	return temp 		  	 	 	


#-----------------------------------calculate brightness temperature------------------------------------------
def brightness_temp(freq,ray_param,cone_angle):
	r_list = np.linspace(1,5,10000)
	angle = 0
	flag = 0
	
	rCritical = singularity(r_list,ray_param,0,freq,flag,cone_angle,angle)
	
	thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , rCritical, 215)[0]
	

	if thetaC < cone_angle:
		flag = 1
		rCritical = singularity(r_list,ray_param,rCritical,freq,flag,cone_angle,angle)
		thetaC = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , rCritical, 215)[0]


	radial_dist1 = np.linspace(rCritical,215,10000) # define the range of values of r
	radial_dist2 = np.linspace(215,rCritical,10000) # define the range of values of r

	observed_temp = 0 	
	background_Temp = 0

	# print rCritical , thetaC

	for i in range(len(radial_dist2)-1):	 
		flag = 0
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist2[i], 215)
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist2[i] , 'theta :',result[0] , 'at frequency', Frequency
		angle = result[0]

		if angle < cone_angle:
			flag = 1		   
			result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist2[i], 215)
			print 'calculating for ray parameter with flag 1 = ',ray_param,'::', i,':' ,'r:',radial_dist2[i] , 'theta :',result[0] , 'at frequency', Frequency			

		angle = result[0]

		tau = opticalDepth(radial_dist2[i+1],radial_dist2[i],freq,ray_param,flag,cone_angle,angle)
		observed_temp = background_Temp*np.exp(-tau) + electron_Temp(radial_dist2[i+1])*(1-np.exp(-tau)) 
		
		print i, ": running from 5 to Rc. Currently r1 = ",radial_dist2[i+1],"r2=",radial_dist2[i],'ray parameter=',ray_param, "tau=",tau,"temp =", observed_temp 
		
		background_Temp = observed_temp	 
	 
	for i in range(len(radial_dist1)-1):

		flag = 0
		result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist1[i], 215)
		print 'calculating for ray parameter = ',ray_param,'::', i,':' ,'r:',radial_dist1[i] , 'theta :',result[0] , 'at frequency', Frequency
		angle = result[0]

		if angle < cone_angle:
			flag = 1		   
			result = integrate.quad(lambda x: DthetaDr(x,ray_param,freq,flag,cone_angle,angle) , radial_dist1[i], 215)
			print 'calculating for ray parameter with flag 1 = ',ray_param,'::', i,':' ,'r:',radial_dist1[i] , 'theta :',result[0] , 'at frequency', Frequency			

		angle = result[0]

		print i, ": running from Rc to 5. Currently r1 = ",radial_dist1[i],"r2=",radial_dist1[i+1],'ray parameter = ', ray_param, "tau=",tau, "temp =", observed_temp
		tau = opticalDepth(radial_dist1[i],radial_dist1[i+1],freq,ray_param,flag,cone_angle,angle)
		observed_temp = background_Temp*np.exp(-tau) + electron_Temp(radial_dist1[i+1])*(1-np.exp(-tau)) 
		background_Temp = observed_temp	 

	print "Brightness temperature" ,observed_temp
	return observed_temp
	


#=================================================================================================================
def plot_trajectory(freq,param_list,cone_angle):
	
	ray_param_list = param_list
	ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
	ax.set_rmax(4)  
	ax.grid(True)  

	for i in range(len(ray_param_list)):
		array = ray_trajectory(ray_param_list[i],freq,cone_angle)
		radial_dist = [value for value in array[0] if value < 10]
		sol1 = array[1][0:len(radial_dist)]
		sol2 = array[2][0:len(radial_dist)]
		sol1 = sol1[~np.isnan(sol1)]
		sol2 = sol2[~np.isnan(sol2)]
		ax.plot(sol2,radial_dist,color='b', linewidth=1)
		sun_theta = np.linspace(0,2*np.pi,1000)
		sun_r = [1]*len(sun_theta)
		ax.plot(sun_theta, sun_r ,color='r', linewidth=1)

	plt.show()


def tempProfile(freq,param_list,cone_angle):
	ray_param_list = param_list
	observed_temp = []
	for a in ray_param_list:
		observed_temp.append(brightness_temp(freq,a,cone_angle))
	plt.plot(param_list, observed_temp)
	plt.show()
	return observed_temp	
	# img = np.zeros((512,512,3), np.uint8)
	# for radius in radii
	# 	cv2.circle(img,(447,63), 63, (0,0,255), -1)


ray_param_list = np.linspace(0,4.95,50)  #define the value of perpendicular distance between the starting point and sun's equator           
# ray_param_list = [1]
cone_angle = 1
tempProfile(Frequency,ray_param_list,cone_angle)

# plot_trajectory(Frequency,ray_param_list,cone_angle)

# r_list = np.linspace(1,5,10000)
# ray_param = 1
# flag = 0
# angle = 0
# r = 1
# for r in r_list:
# print eDensity(2.62,flag,cone_angle,angle)
# print sqr_ref_index(r,Frequency,flag,cone_angle,angle)
# print singularity(r_list,ray_param,r,Frequency,flag,cone_angle,angle)
# print DthetaDr(r,Frequency,flag,cone_angle,angle,ray_param) 
