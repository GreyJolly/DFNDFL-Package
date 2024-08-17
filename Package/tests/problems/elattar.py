# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'elattar'
startp = np.array([2.0,2.0,7.0,0.0,-2.0,1.0])
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
    t = np.zeros(51)
    y = np.zeros(51)
    f = np.zeros(102)
    for i in range(1,52):
        t[i-1] = (i - 1.0) / 10.0
        y[i-1] = np.exp(t[i-1])/2.0 - np.exp(-2.0*t[i-1])
        y[i-1] = y[i-1] + np.exp(-3.0*t[i-1])/2.0
        y[i-1] = y[i-1] + 1.5*np.exp(-1.5*t[i-1])*np.sin(7.0*t[i-1])
        y[i-1] = y[i-1] +     np.exp(-2.5*t[i-1])*np.sin(5.0*t[i-1])
        f[i-1] = x[0]*np.exp(-x[1]*t[i-1])*np.cos(x[2]*t[i-1] + x[3])
        f[i-1] = f[i-1] + x[4]*np.exp(-x[5]*t[i-1]) - y[i-1]

    for i in range(52,103):
        f[i-1] = -f[i-52]
    
    y = np.max(f)
    
    return y
