""" Maria Lavrovskaya: Numerical Methods Assignment 1 - Question 2.2 """


import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from scipy.stats import norm
import timeit
import time
import math
from scipy.optimize import curve_fit

from ex2_1 import matrix_B

#Finding the time it takes to calculate the matmul for (a) and (b) relatively to the number of entries in matrices

nvec1 = np.array([10,20,50,100,200,500,1000, 2000, 5000])
CPUt1 = np.zeros(len(nvec1))
seed = 10
np.random.seed(seed)
for index in range(len(nvec1)): 
    
    A=np.random.rand(nvec1[index],nvec1[index])
    B=np.random.rand(nvec1[index],nvec1[index])
    start_time = time.time()
    Ca = np.matmul(A,B)
    end_time = time.time()
    CPUt1[index]=end_time-start_time

nvec2 = np.array([10, 20, 50, 100, 200, 500, 1000])
CPUt2 = np.zeros(len(nvec2))

for index in range(len(nvec2)): 

    A=np.random.rand(nvec2[index],nvec2[index])
    B=np.random.rand(nvec2[index],nvec2[index])
    start_time = time.time()
    matrix_B(A,B)
    end_time = time.time()
    CPUt2[index]=end_time-start_time


#Curve fitting
#Objective function to get the straight line between inputs and outputs
def objective(x,a,b):
    return a*x+b
#Out values on log10 scale, since we are about to plot loglog graph
x_values = np.log10(nvec2)
y_values = np.log10(CPUt2)

#Fitting curve and getting the parameters
param, _ = curve_fit(objective, x_values, y_values)
#New values we are trying to fit
x_line = x_values #np.linspace(nvec2[0], nvec2[6], num=7)
a,b = param
#Calculating new values for our data 
y_new = 10**objective(x_line, a,b)

#------ Calculating for matmul
x_nvec1 = np.log10(nvec1)
y_CPUt1  = np.log10(CPUt1)

param1, _ = curve_fit(objective, x_nvec1, y_CPUt1)
x_line1 = x_nvec1
a1,b1 = param1
y_new1 = objective(x_line1, a1,b1)

plt.figure(figsize = (12,8))
plt.loglog(10**x_line, y_new, '--', color='green', label = 'Triple loop')
plt.loglog(10**x_line1, 10**y_new1, '--', color='b', label = 'Numpy Matmul')
plt.xlabel('n')
plt.ylabel('CPU t/s')
plt.legend(loc = 'best')
plt.title('Matrix multiplication with triple loop', fontsize =10, fontweight= 'bold', fontname = 'DejaVu Sans')

plt.savefig('artworks/q2.png')
