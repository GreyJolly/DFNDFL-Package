# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'gill'
startp = -0.1*np.ones(10)
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 5
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    #   print(x)
    fx = np.zeros(3)
    fx[0] = np.sum((x-1.0)**2) + 0.001*np.sum((x**2 - 0.25)**2)
    
    j = np.arange(2,11)
    i = np.arange(2,31)
    A = np.tile((x[j-1]*(j-1)).reshape(-1,1),(1,29))*np.tile(((i-1)/29),(9,1))**np.tile((j-2).reshape(-1,1),(1,29))
    #print(A)
    
    j = np.arange(1,11)
    B = np.tile(((i-1)/29),(10,1))**np.tile((x*(j-1)).reshape(-1,1),(1,29))
    fx[1] = np.sum((np.sum(A,0) - (np.sum(B,0))**2 - 1.0)**2) + x[0]**2 + (x[1] - x[0]**2 - 1.0)**2
    
    i = np.arange(2,11)
    fx[2] = np.sum(100.0*(x[i-1] - x[i-2]**2)**2 + (1.0 - x[i-1])**2)
    
    y = np.max(fx)
    return y

