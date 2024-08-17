#**************************************************
# prob. n.13 described in paper:
# J. Müller, C.A. Shoemaker, R. Piché
# SO-I: a surrogate model algorithm for expensive nonlinear
# integer programming problems including global optimization applications
# Journal of Global Optimization, 59(4):865-889 (2014)
# 
# N.B. variables u_i are x(i), i = 1..nu
#      variables x_i are x(r+i), i = 1..nx
#**************************************************
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'SO-I prob. 13'
n      = 10
nint   = 5
ncont  = n-nint
lb     = 3.0*np.ones(n)
ub     =99.0*np.ones(n)
lbmix  = 3.0*np.ones(n)
ubmix  =99.0*np.ones(n)
startp =   51.0*np.ones(n) 
x_initial =51.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  #*x[0] *x[1] *x[2] x[3] x[4]
    f = np.sum(np.log(x-2)**2)*np.sum(np.log(100-x)**2) - np.prod(x**0.2)
    return f

