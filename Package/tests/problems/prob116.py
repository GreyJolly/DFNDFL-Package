#**************************************************
# prob. n.16 described in paper:
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

name      = 'SO-I prob. 16'
n      = 8
nint   = 4
ncont  = n-nint
lb     =-10.0*np.ones(n)
ub     = 10.0*np.ones(n)
lbmix  =-10.0*np.ones(n)
ubmix  = 10.0*np.ones(n)
startp =   0.0*np.ones(n) 
x_initial =0.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  #*x[0] *x[1] *x[2] x[3] x[4]
    f = 3.1*x[0]**2 + 7.6*x[1]**2 + 6.9*x[2]**2 + 0.004*x[3]**2 + 19*x[4]**2 + 3*x[5]**2 + x[6]**2 + 4*x[7]**2
    return f

