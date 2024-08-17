# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'exp'
startp = np.array([0.5, 0.0, 0.0, 0.0, 0.0])
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
    t = -1 + (np.arange(1,22)-1)/10     
    f = (x[0] + x[1]*t)/(1.0 + x[2]*t + x[3]*t**2 + x[4]*t**3) - np.exp(t)
    y = np.max(f)
    return y
