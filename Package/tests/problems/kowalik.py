# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'kowalik-osborne'
startp = np.array([0.25, 0.39, 0.415, 0.39])
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
    z = np.array([[0.1957,4.0],
     [0.1947,2.0],
     [0.1735,1.0],
     [0.1600,0.5],
     [0.0844,0.25],
     [0.0627,0.1670],
     [0.0456,0.1250],
     [0.0342,0.1],
     [0.0323,0.0833],
     [0.0235,0.0714],
     [0.0246,0.0625]])
    u = z[:,1]
    f = (x[0]*(u**2 + x[1]*u))/(u**2 + x[2]*u + x[3]) - z[:,0]
    y = np.max(np.abs(f));
    return y
