#**************************************************
# prob. n.8 described in paper:
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

name      = 'MISO prob. 8'
n      = 15
nint   = 10
ncont  = n-nint
lb     =-15.0*np.ones(n)
ub     = 30.0*np.ones(n)
lbmix  =-15.0*np.ones(n)
ubmix  = 30.0*np.ones(n)
startp =   7.0*np.ones(n) 
x_initial =7.0*np.ones(n) 
xmix   = np.zeros(n)

def feval(x):  
    # continue 5  da 11 a 15  --> (per noi) da 0 a 4
    # discrete 10 da  1 a 10  --> (per noi) da 5 a 14  
    #f = (x(1) - 1)^2 + sum([2:15]'.*(2*x(2:15).^2 - x(1:14)).^2);
    f = (x[5]-1)**2 + np.sum(np.arange(1,14)*(2*x[1:14]**2 - x[0:13]**2))
    return f
