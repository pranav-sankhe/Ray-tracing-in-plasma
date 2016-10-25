import numpy as np
import simpy 
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(1500)

n = 1000
vec_x = np.arange(n)
#vec_y = np.arange(n)

r_index = np.ones(n)

for i in range(1000) :
	r_index[i] = r_index[i]*2 

in_angle  = np.pi/4
in_height = 10 
vec_y = np.array([])

def ang_of_incidence(n):
    if n == 1:
    	incident_ang = r_index[0]/r_index[1] 
        return incident_ang
    else:
    	incident_ang = (r_index[n-1]/r_index[n])*ang_of_incidence(n-1) 
        return  incident_ang

def path_eq(n):
	path_of_ray = np.array([]) 
	if n  == 0:
		path_of_ray[0] = in_height
		return path_of_ray
	else:
		path_of_ray[n] = path_eq(n-1) + np.tan(np.arcsin(np.sin(in_angle)*ang_of_incidence(n))) 
		return path_of_ray

vec_y = path_eq(n)

plt.plot(vec_x, vec_y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()




