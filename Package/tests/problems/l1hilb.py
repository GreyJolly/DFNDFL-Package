# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'l1hilb'
startp = np.ones(50)
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 25
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    i = np.arange(1,51)
    j = np.arange(1,51)
    I = np.tile(i,(50,1));
    J = np.tile(j.reshape(-1,1),(1,50));
    X = np.tile(x.reshape(-1,1),(1,50));
    y = np.sum(np.abs(np.sum(X/(I+J-1))))
    return y
