# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'banex'
n      = 2
startp = np.zeros(n)
lb     = np.array([-50.0, 0.0])
ub     = np.array([10.0, 100.0])
nint   = 2
ncont  = n-nint
lbmix  = np.zeros(n)
ubmix  = 100.0*np.ones(n)
x_initial = 50*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):
    return (x[0]-1.0)**2 + 100.0*(x[0]**2-x[1])**2
