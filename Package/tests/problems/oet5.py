# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'oet5'
startp = np.ones(4)
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
    t = 0.25 + 0.75*(I-1)/20;
    f = x[3] - (x[0]*t**2 + x[1]*t + x[2])**2 - np.sqrt(t);
    y = np.max(np.abs(f));
    return y
