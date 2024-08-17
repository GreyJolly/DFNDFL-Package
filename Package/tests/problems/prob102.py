#**************************************************
# prob. n.2 described in paper:
# J. Müller, C.A. Shoemaker, R. Piché
# SO-I: a surrogate model algorithm for expensive nonlinear
# integer programming problems including global optimization applications
# Journal of Global Optimization, 59(4):865-889 (2014)
# 
# N.B. variables u_i are x(i), i = 1..nu
#      variables x_i are x(r+i), i = 1..nx
#**************************************************
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'SO-I prob. 2'
#devono essere continue le prime 3
n      = 5
nint   = 2
ncont  = n-nint
lb     = -100.0*np.ones(n)
ub     =  100.0*np.ones(n)
lbmix  = -100.0*np.ones(n)
ubmix  =  100.0*np.ones(n)
startp =   50.0*np.ones(n) 
x_initial =50.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  #*x[0] *x[1] *x[2] x[3] x[4]
    f = x[3]*np.sin(x[3]) + 1.7*x[4]*np.sin(x[3]) - 1.5*x[0] - 0.1*x[1]*np.cos(x[1] + x[2] - x[3]) + 0.2*x[2]**2 - x[4] -1.0
    return f
