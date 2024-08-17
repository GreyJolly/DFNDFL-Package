# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'polak 2'
startp = 0.1*np.ones(10)
startp[0] = 100.0
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
    e2 = np.zeros(10); e2[1] = 1.0
    w = x + 2*e2
    f1 = np.exp(1e-8*w[0]**2 + w[1]**2 + w[2]**2 + 4*w[3]**2 + w[4]**2 + w[5]**2 + w[6]**2 + w[7]**2 + w[8]**2 + w[9]**2)
    w = x - 2*e2
    f2 = np.exp(1e-8*w[0]**2 + w[1]**2 + w[2]**2 + 4*w[3]**2 + w[4]**2 + w[5]**2 + w[6]**2 + w[7]**2 + w[8]**2 + w[9]**2)    
    y = np.maximum(f1,f2)    
    return y
