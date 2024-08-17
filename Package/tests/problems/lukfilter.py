# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'filter'
startp = np.array([0.0, 1.0, 0.0, -0.15, 0.0, -0.68, 0.0, -0.72, 0.37])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 4
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    t = np.zeros(41)
    t[0:6] = 0.01*(np.arange(1,7)-1)
    t[6:20] = 0.07 + 0.03*(np.arange(7,21)-7);
    t[20]   = 0.5;
    t[21:35] = 0.54 + 0.03*(np.arange(22,36)-22);
    t[35:41] = 0.95 + 0.01*(np.arange(36,42)-36);

    z = np.abs(1-2*t) 
    eta = np.pi*t    
    A = ((x[0] + (1+x[1])*np.cos(eta))**2 + ((1-x[1])*np.sin(eta))**2)/((x[2] + (1+x[3])*np.cos(eta))**2 + ((1-x[3])*np.sin(eta))**2)        
    B = ((x[4] + (1+x[5])*np.cos(eta))**2 + ((1-x[5])*np.sin(eta))**2)/((x[6] + (1+x[7])*np.cos(eta))**2 + ((1-x[7])*np.sin(eta))**2)    
    f = x[8]*np.sqrt(A)*np.sqrt(B) - z;    
    y = np.max(np.abs(f));
    return y
