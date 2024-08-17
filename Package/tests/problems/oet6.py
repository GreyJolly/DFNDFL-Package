# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'oet6'
startp = np.array([1.0, 1.0, -3.0, -1.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 2
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    I = np.arange(1,22).reshape(21,-1)
    t = -0.5 + (I-1)/20;
    f = x[0]*np.exp(x[2]*t) + x[1]*np.exp(x[3]*t) - (t+1)**(-1)
    y = np.max(np.abs(f))
    return y
