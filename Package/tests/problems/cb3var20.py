# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'cb3(20)'
n      = 20
startp = 2.0*np.ones(n)
lb     = startp - 10.0
ub     = startp + 10.0
nint   = 10
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    y = 0.0
    for i in range(n-1):
        fx1 = x[i]**4+x[i+1]**2
        fx2 = (2-x[i])**2+(2-x[i+1])**2
        fx3 = 2.0*np.exp(-x[i]+x[i+1])
        y += np.max([fx1,fx2,fx3])
        
    return y
