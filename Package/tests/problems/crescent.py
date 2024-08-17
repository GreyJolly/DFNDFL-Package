# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'crescent'
n      = 2
startp = np.array([-1.5, 2.0])
lb     = startp - 10.0
ub     = startp + 10.0
nint   = 2
ncont  = n-nint
lbmix  = np.zeros(n)
ubmix  = 100.0*np.ones(n)
x_initial = 50*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):
    f1 =   x[0]**2 + (x[1]-1.0)**2 + x[1] - 1.0
    f2 = - x[0]**2 - (x[1]-1.0)**2 + x[1] + 1.0

    return np.maximum(f1,f2)
