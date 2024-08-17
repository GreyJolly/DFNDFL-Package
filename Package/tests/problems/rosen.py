# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'rosen-suzuki'
startp = np.array([0.0,0.0,0.0,0.0])
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
    f[0] =  x[0]**2 + x[1]**2 + 2.0*x[2]**2 + x[3]**2 
    f[0] =  f[0] - 5.0*( x[0] + x[1] ) - 21.0*x[2] + 7.0*x[3] 
    f[1] =  x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2
    f[1] =  f[1] + x[0] - x[1] + x[2] - x[3] - 8.0
    f[2] =  x[0]**2 + 2.0*x[1]**2 + x[2]**2 + 2.0*x[3]**2
    f[2] =  f[2] - x[0] - x[3] - 10.0
    f[3] =  x[0]**2 + x[1]**2 + x[2]**2 + 2.0*x[0] 
    f[3] =  f[3] - x[1] - x[3] - 5.0

    f[1] = 10.0*f[1] + f[0]
    f[2] = 10.0*f[2] + f[0]
    f[3] = 10.0*f[3] + f[0]
    
    y = np.max(f)
    return y
