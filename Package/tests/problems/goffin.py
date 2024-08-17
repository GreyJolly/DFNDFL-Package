# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'goffin'
startp = np.arange(1,51)-25.5
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 25
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
#x_initial = np.array([-14.5, -13.5, -12.5, -11.5, -10.5, -9.5, -8.5, -7.5, 
#                       -6.5,  -5.5,  -4.5,  -3.5,  -2.5, -1.5, -0.5,  0.5, 
#                        1.5,   2.5,   3.5,   4.5,   5.5,  6.5,  7.5,  8.5, 
#                        9.5,   100,   100,   100,   100,  100,  100,  100, 
#                        100,   100,    99,    95,    90,   85,   80,   75, 
#                         70,    65,    60,    55,    50,   45,   40,   35, 30, 25])
xmix   = np.zeros(n)

def feval(x):
    f = np.zeros(50)
    somma = np.sum(x)
    for i in range(0,50):
        f[i] = 50.0*x[i] - somma
    y = np.max(f);
    return y