# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'cb2'
n      = 2
startp = np.array([2.0, 2.0])
lb     = startp - 10.0
ub     = startp + 10.0
nint   = 2
ncont  = n-nint
lbmix  = np.zeros(n)
ubmix  = 100.0*np.ones(n)
x_initial = 50*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):
    fx1 = x[0]**2+ x[1]**4
    fx2 = (2-x[0])**2 + (2-x[1])**2
    fx3 = 2*np.exp(x[1]-x[0])
    return np.max([fx1,fx2,fx3])
