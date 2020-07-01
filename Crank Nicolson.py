import numpy as np
import scipy as sp
from scipy import linalg

def CrankNicolson(l, N, m, T, alpha, sigma, f):
    '''
    Trapezoidal(Crank Nicolson) in time, Central difference in space
    l: total number of space steps
    N: endpoint in space
    m: total number of time steps
    T: endpoint in time
    alpha: initial condition
    sigma: given in heat equation
    f: intial condition
    '''
    k= T/m
    print('k = ', k)
    h= N/l
    print('h = ', h)
    lambd = sigma*(k/(2*(h**2)))
    print('lambda=', lambd)
    numList = []
    A = np.zeros([l+1, l+1])
    B = np.zeros([l+1, l+1])
    
    for i in range(1, l):
        A[i,i-1] = -lambd
        A[i,i+1] = -lambd
        A[i,i] = 1 + (2*lambd)
        B[i,i-1] = lambd
        B[i,i+1] = lambd
        B[i,i] = 1 - (2*lambd)
    A[0,0] = A[-1,-1] = 1
    B[0,0] = B[-1,-1] = 1
    #print('A=', A)
    #print('B=', B)
    
    U = alpha #Vector valued
    
    for j in range(0, m): 
        #Through Time
        
        #Through Space
        newU = sp.linalg.solve(A, np.matmul(B, U))
        newU[0] = newU[-1] = 0
        #numList.append(newU)
        #print('newU=', newU)
        U=newU
        
    #print('U=', U)    
    return newU #numList
