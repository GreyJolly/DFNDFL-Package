#**************************************************
# prob. n.6 described in paper:
# J. Müller
# MISO: Mixed-Integer Surrogate Optimization Framework
# Optimization and Engineering, 17(1):177-203 (2016)
# 
# N.B. variables u_i are x(i), i = 1..nu
#      variables x_i are x(r+i), i = 1..nx
#**************************************************
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'MISO prob. 6'
n      = 15
nint   = 7
ncont  = n-nint
lb     =-15.0*np.ones(n)
ub     = 30.0*np.ones(n)
lbmix  =-15.0*np.ones(n)
ubmix  = 30.0*np.ones(n)
startp =   0.0*np.ones(n) 
x_initial =7.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  #*x[0] *x[1] *x[2] x[3] x[4]
    f = - 20*np.exp(-0.2*np.sqrt(np.sum(x**2)/15)) - np.exp(np.sum(np.cos(2*np.pi*x))/15)
    return f

