# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'shell dual'
startp = 1e-4*np.ones(15)
startp[6] = 60.0
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 7
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)


def feval(x):
    x = x.reshape(-1,1)
    A = np.array([[-16.0,2.0,0.0,1.0,0.0],
        [0.0,-2.0,0.0,4.0,2.0],
        [-3.5,0.0,2.0,0.0,0.0],
        [0.0,-2.0,0.0,-4.0,-1.0],
        [0.0,-9.0,-2.0,1.0,-2.8],
        [2.0,0.0,-4.0,0.0,0.0],
        [-1.0,-1.0,-1.0,-1.0,-1.0],
        [-1.0,-2.0,-3.0,-2.0,-1.0],
        [1.0,2.0,3.0,4.0,5.0],
        [1.0,1.0,1.0,1.0,1.0]])
 
    b = np.array([[-40.0],[-2.0],[-0.25],[-4.0],[-4.0],[-1.0],[-40.0],[-60.0],[5.0],[1.0]])
    
    C = np.array([[30.0,-20.0,-10.0,32.0,-10.0],
         [-20.0,39.0,-6.0,-31.0,32.0],
         [-10.0,-6.0,10.0,-6.0,-10.0],
         [32.0,-31.0,-6.0,39.0,-20.0],
         [-10.0,32.0,-10.0,-20.0,30.0]])
     
    d = np.array([[4.0],[8.0],[10.0],[6.0],[2.0]])
    
    ee = np.array([-15.0,-27.0,-36.0,-18.0,-12.0])
    
    J10 = np.arange(1,11)
    J5  = np.arange(1,6)
    
    Q = np.sum(np.minimum(0,x))
    P = np.matmul(A.T,x[J10-1]) - 2*C*x[J5+10-1] - 3*d*x[J5+10-1]**2 - ee.T
    X = np.tile(x[J5+10-1],(1,5))
    
    y = 2*np.abs(np.sum(d*x[J5+10-1]**3)) + np.sum(np.sum(C*X*X.T)) - np.sum(b*x[J10-1]) + 100*(np.sum(np.maximum(0,P)) - Q)
    return y
