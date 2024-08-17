# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'wong1'
startp = np.array([1.0,2.0,0.0,4.0,0.0,1.0,1.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 3
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    x = x.reshape(-1,1)
    f = np.zeros(5)
    f[0] = (x[0]-10.0)**2 +5.0*(x[1]-12.0)**2 + x[2]**4 + 3.0*(x[3]-11.0)**2 + 10.0*x[4]**6+ 7.0*x[5]**2 + x[6]**4 - 4.0*x[5]*x[6] - 10.0*x[5] - 8.0*x[6]
    f[1] = f[0] + 10.0*(2.0*x[0]**2 + 3.0*x[1]**4 + x[2] + 4.0*x[3]**2 + 5.0*x[4] -127.0)
    f[2] = f[0] + 10.0*(7.0*x[0] + 3.0*x[1] + 10.0*x[2]**2 + x[3] - x[4] - 282.0)
    f[3] = f[0] + 10.0*(23.0*x[0] + x[1]**2 + 6.0*x[5]**2 - 8.0*x[6] -196.0)
    f[4] = f[0] + 10.0*(4.0*x[0]**2 + x[1]**2 - 3.0*x[0]*x[1] + 2.0*x[2]**2 + 5.0*x[5] - 11.0*x[6])
    y = np.max(f)
    return y
