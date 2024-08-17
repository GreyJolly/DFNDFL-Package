# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'hs78'
startp = np.array([-2.0, 1.5, 2.0, -1.0, -1.0])
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
    fx = np.zeros(3)
    fx[0] = np.sum(x**2) - 10.0
    fx[1] = x[1]*x[2] -5.0*x[3]*x[4]
    fx[2] = x[0]**3 + x[1]**3 + 1.0
    y= np.prod(x) + 10.0*np.sum(np.abs(fx));
    return y
