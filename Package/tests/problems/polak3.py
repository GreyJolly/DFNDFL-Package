# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'polak 3'
startp = np.ones(11)
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 5
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    J = np.tile(np.arange(0,11).reshape(-1,1),(1,10))
    I = np.tile(np.arange(1,11),(11,1))
    f = np.sum((J+I)*np.exp((np.tile(x.reshape(-1,1),(1,10)) - np.sin(I-1+2*J))**2),0)
    y = np.max(f)
    return y
