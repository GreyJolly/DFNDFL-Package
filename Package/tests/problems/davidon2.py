# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'davidon2'
startp = np.array([25.0,5.0,-5.0,-1.0])
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
    t = np.zeros(21)
    f = np.zeros(21)
    for i in range(1,22):
        t[i-1] = 0.25 + 0.75*(i-1)/20.0
        f[i-1] = x[3] - (x[0]*t[i-1]**2 + x[1]*t[i-1] + x[2])**2 - np.sqrt(t[i-1])
    y = np.max(np.abs(f));
    return y
