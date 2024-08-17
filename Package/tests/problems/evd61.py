# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'evd61'
startp = np.array([2.0, 2.0, 7.0, 0.0, -2.0, 1.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 3
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    t = 0.1*(np.arange(1,52) - 1)
    z = 0.5*np.exp(-t) - np.exp(-2.0*t) + 0.5*np.exp(-3.0*t) + 1.5*np.exp(-1.5*t)*np.sin(7*t) + np.exp(-2.5*t)*np.sin(5*t) 
    f = x[0]*np.exp(-x[1]*t)*np.cos(x[2]*t + x[3]) + x[4]*np.exp(-x[5]*t) - z
    y = np.max(np.abs(f))
    return y
