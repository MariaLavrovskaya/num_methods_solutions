""" Maria Lavrovskaya: Numerical Methods Assignment 1 - Question 2.1 """



import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from scipy.stats import norm
import timeit
import time
import math
from scipy.optimize import curve_fit

np.random.seed(5)
A=np.random.rand(3,3)
B=np.random.rand(3,3)


#Compute the matrix product using the built-in matrix product function (a)
Ca = np.matmul(A,B)
print(Ca)

#Compute C=AB from the definition (b)
def matrix_B(matrix1, matrix2):
    A=matrix1
    B=matrix2
    Cb=np.zeros((A.shape[0],B.shape[1]))
    #iterate through rows of X
    for j in range(len(B[0])):
        #iterate through the columns of Y
        for i in range(len(A)):
            #iterate through rows of Y
            for k in range(len(B)):
                Cb[i][j]+=A[i][k]*B[k][j]
    return Cb


#Repeat (b) swapping the two outerloops over i and j (c)

def matrix_C(matrix1, matrix2):
    A = matrix1
    B = matrix2
    #Repeat (b) swapping the two outerloops over i and j (c)
    Cc=np.zeros((A.shape[0],B.shape[1]))
    #iterate through rows of X
    for j in range(len(B[0])):
        #iterate through the columns of Y
        for i in range(len(A)):
            #iterate through rows of Y
            for k in range(len(B)):
                Cc[i][j]+=A[i][k]*B[k][j]
    return Cc

#Repeat (b) with two for loops over i and j and sum() in place of k (d)

def matrix_D(matrix1, matrix2):
    A=matrix1
    B=matrix2
    Cd=np.zeros((Ca.shape))
    for i in range(len(A)):
        #iterate through the columns of Y
        for j in range(len(B[0])):
            #iterate through rows of Y
            Cd[i][j]+=np.sum(A[i,:]*B[:,j])
    return Cd


#Repeat (d) transposing A to increase the data locality and thus the cache hits (e)


#check if that works if we multiply not square matrix
def matrix_E(matrix1, matrix2):
    A=matrix1
    B=matrix2.T
    Ce=np.zeros((A.shape[0],B.shape[1]))
    for i in range(len(A)):
        #iterate through the columns of Y
        for j in range(len(B[0])):
            #iterate through rows of Y
            Ce[i][j]+=np.sum(A[i,:]*B[j,:])
    return Ce


print('Matrix Ca: ', Ca, '\n')
print('Matrix Cb: ' , matrix_B(A,B), '\n')
print('Matrix Cc: ', matrix_C(A,B), '\n') 
print('Matrix Cd:', matrix_D(A,B), '\n')
print('Matrix Ce:', matrix_E(A,B), '\n')
