from __future__ import division 
import numpy as np
import simpy 
import matplotlib.pyplot as plt
import sys
import math 
from scipy.optimize import fsolve

sun_radii = 696342000.0   #gloabla constant


in_xCor =  np.sqrt(2)*sun_radii  #float
in_yCor =  np.sqrt(2)*sun_radii  #float
in_r =   np.power( np.power(in_xCor,2) + np.power(in_yCor,2) ,0.5)

in_theta = 1.5*np.pi - math.atan(in_xCor/in_yCor)
step_theta = -0.001 
fin_theta = np.pi*0.5	
number_of_samples = abs(int((in_theta - fin_theta)/step_theta ))
n = number_of_samples
theta = np.linspace(in_theta, fin_theta ,num=number_of_samples )
delta_theta = theta[3] - theta[2]

print 'no of iterations:' , n  


#0.273727727728

def sqr_ref_index(r,freq):
	r = r/sun_radii
	val = 1 - (12400/(freq*freq*np.power(r,6)))*(1 + 1.93/np.power(r,10))
	return val
	


# def get_Rmaxmin():
# 	in_r = 0
# 	dist = np.linspace(0.274*sun_radii , 4*sun_radii,n)
# 	for i in range(n):
# 		if ref_index(dist[i]) > 0.9999999999999:
# 			in_r = dist[i]
# 			break

# 	min_r = 0.273727727728*sun_radii
# 	return [min_r, in_r]

# # print get_Rmaxmin()[0]/sun_radii, get_Rmaxmin()[1]/sun_radii
# # print ref_index( get_Rmaxmin()[1]), ref_index( get_Rmaxmin()[0])
# r_min = get_Rmaxmin()[0]
# r_max = get_Rmaxmin()[1]


# r'(theta) = omega(theta)
#r = distance/sun_radius
def der_omega(omega,r,freq):
	val1 = (1 + r*r)/omega
	val2 = (3*12400)/(freq*freq*np.power(r,4))
	val3 = (1.93*5*12400)/(freq*freq*np.power(r,14))
	val4 = sqr_ref_index(r,freq)

	val = omega + (val1*(val2+val3))/val4
	return val

def pend(y,theta,freq,one):
    r, omega = y
    dydt = [omega, der_omega(omega,r,freq)*one*(theta/theta) ]
    return dydt

freq = 5*np.power(10,6)
one = 1
r0 = [in_theta,in_r/np.tan(in_theta)]  #first element  = theta . 2nd element = omega
# initial value of omega can be calculated using polar differntial eqns 
 


from scipy.integrate import odeint
sol = odeint(pend, r0, theta, args=(freq,one))

for solution,ang in zip(sol[:,0], theta):
	print ang, solution



#plt.plot(t, sol[:, 0], 'b', label='theta(t)')	
ax = plt.subplot(111, projection='polar')    # add subplot in polar coordinates 
ax.set_rmax(5)  
ax.grid(True)  
#ax.legend(loc='best')
ax.plot(theta,abs(sol[:,0]),color='b', linewidth=3)
plt.show()





'''
for i in range(n-1):
	#1st iteration will yield dr1 and so on.....  
	
	func = lambda delta_r : np.arctan(-(delta_theta/(delta_r))*(delta_r + r) ) - delta_theta - np.arcsin( np.sin(prev_alpha)*ref_index(r - dr)/ref_index(r) )

	deltaR_initial_guess = 100000
	dr = fsolve(func, deltaR_initial_guess)
	prev_alpha = math.atan( -(delta_theta/(dr))*(dr + r) )	
	
	# solve the equation. get the value of dr 
	r = r + dr
	r_array[i+1] = r
	print i, "distance",  r ,r/sun_radii,"increment", dr ,'theta' , theta[i] 
	if (r - r_min) < 10e-10:
		break
'''