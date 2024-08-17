# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'maxq(40)'
startp = np.concatenate([np.array(np.arange(1,21)),np.array(np.arange(-21,-41,-1))])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 20
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    f = np.zeros(40)
    for i in range(0,40):
        f[i] = x[i]**2
    y = np.max(f);
    return y
