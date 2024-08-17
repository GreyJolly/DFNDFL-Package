# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'pbc1'
startp = np.array([0.0, -1.0, 10.0, 1.0, 10.0])
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
    t = -1 + 2*(np.arange(1,31)-1)/29
#    f = np.divide((x[0] + x[1]*t + x[2]*t**2),(1.0 + x[3]*t + x[4]*t**2)) - np.divide((np.sqrt((8*t - 1.0)**2 + 1.0)*np.arctan(8*t)),(8.0*t))
    f = (x[0] + x[1]*t + x[2]*t**2)/(1.0 + x[3]*t + x[4]*t**2) - (np.sqrt((8*t - 1.0)**2 + 1.0)*np.arctan(8*t))/(8.0*t)

    y = np.max(np.abs(f))
    return y
