# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'watson'
startp = np.zeros(20)
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 10
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    x = x.reshape(-1,1)
    f = np.zeros(31)
    I = np.tile(np.arange(3,32),(20,1))
    J = np.tile(np.arange(1,21).reshape(-1,1),(1,29)) 
    X = np.tile(x,(1,29))
    f[0:29] = np.sum((J-1)*X*((I-2)/29)**(J-2),axis=0) - np.sum(X*((I-2)/29)**(J-1),axis=0)**2
    f[29] = x[0]
    f[30] = x[1] - x[0]**2 - 1
    y = np.max(np.abs(f))
    return y
