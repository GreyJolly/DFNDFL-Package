# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'polak 6'
startp = np.zeros(4)
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
    f = np.zeros(4)
    a    = x[0] - (x[3]+1.0)**4;
    f[0] = a**2 + (x[1]-a**4)**2 + 2.0*x[2]**2 + x[3]**2 - 5.0*a - 5.0*(x[1]-a**4) - 21.0*x[2] + 7.0*x[3]
    f[1] = f[0] + 10.0*(a**2 + (x[1]-a**4)**2 + x[2]**2 + x[3]**2 + a - (x[1]-a**4) + x[2] - x[3] - 8.0)
    f[2] = f[0] + 10.0*(a**2 + 2.0*(x[1]-a**4)**2 + x[2]**2 + 2.0*x[3]**2 - a - x[3] - 10.0)
    f[3] = f[0] + 10.0*(a**2 + (x[1]-a**4)**2 + x[2]**2 + 2.0*a - (x[1]-a**4) -x[3] - 5.0)
    y = np.max(f)
    return y
