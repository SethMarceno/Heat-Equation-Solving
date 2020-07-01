import numpy as np
import scipy as sp
from scipy import linalg


f = lambda x, y:  #Function of (x, y) given by initial condition

def ADI(l, N, m, T, sigma, f, time1, time2, time3):
    
    '''
    l=endpoint in space in both x and y directions
    N=total number of steps in space
    m=endpoint in time
    T=total number of steps in time
    alpha=initial condition given u(0, x, y)
    sigma=given in 2-D Heat equation
    Time1, Time2, Time3: Threee times in which to return the approximations
    '''
    
    h = l/N
    k = m/T
    lambd = (sigma*k)/(2*(h**2))
    print('h=', h)
    print('k=', k)
    print('lambda=', lambd)
    numList = []
    Ustar = np.zeros([N+1, N+1])
    U1 = np.zeros([N+1, N+1])
    
    
    A = np.zeros([N+1, N+1])
    alpha= np.zeros([N+1, N+1])
    
    A[0,0] = A[N, N] = 1
    
    for i in range(1, N):
        A[i,i-1] = -lambd
        A[i,i+1] = -lambd
        A[i,i] = 1 + 2*lambd
    #print('A=', A)    
    
    
    #initializing our first matrix
    for y in np.arange(1, N):
        for x in np.arange(1, N):
            alpha[y, x] = f((x*h), (y*h))
    #print('alpha=', alpha)
    numList.append(alpha)
    
    Uvec=np.zeros([N+1])
    
    #Through Time
    for time in np.arange(0, m, k):
        
        
        # Step-1:
        U0 = alpha
        for col in range(1, N):
            b = np.zeros([N+1])
            for row in range(1, N):
                b[row] = (lambd*U0[row, col+1]) + ((1-(2*lambd))*U0[row, col]) + (lambd*U0[row, col-1])


            Ustar[:, col] =  sp.linalg.solve(A, b) #assigning a column to U1

        #print('Ustar', Ustar)

        # Step-2:
        for row in range(1, N):
            b = np.zeros([N+1])
            for col in range(1, N):
                b[col] = (lambd*Ustar[row+1, col]) + ((1-(2*lambd))*Ustar[row, col]) + (lambd*Ustar[row-1, col])

            U1[row, :] = sp.linalg.solve(A, b) #assigning a row to U1
         
        #print('U1', U1)
        
        
        if np.round(time, 3) == time1:
            T1 = U1.copy()
        if np.round(time, 3) == time2:
            T2 = U1.copy()
        if np.round(time, 3) == time3:
            T3 = U1.copy()

        
        alpha = U1
    
    
    return T1, T2, T3
