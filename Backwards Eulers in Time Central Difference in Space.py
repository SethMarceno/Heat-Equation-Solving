import numpy as np
import scipy as sp
from scipy import linalg

def BECS(l, N, m, T, alpha, sigma, f):
    '''
    Backwards Eulers in time, Central difference in space
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
    lambd = sigma*(k/(h**2))
    numList = []
    A = np.zeros([l+1, l+1])
    A[0,0] = A[l,l] = 1
    
    for i in range(1, l):
        A[i,i-1] = -lambd
        A[i,i+1] = -lambd
        A[i,i] = 1 + 2*lambd
    #print('A=',A)    


    
    U = alpha #Vector valued
    Uvec = np.zeros([l+1])
    
    for j in range(0, m): 
        #Through Time
        
        for i in range(0, l+1): 
            #Through space
            if i == 0:
                Uvec[i] = 0
            elif i == 100:
                Uvec[i] = 0
            else:
                Uvec[i] = -U[i]
                

        #numList.append(Uvec)
        U[:]=sp.linalg.solve(A, Uvec)
        U[0]=U[-1]=0
        
    return U #numList
