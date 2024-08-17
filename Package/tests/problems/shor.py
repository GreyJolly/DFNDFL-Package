# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'shor'
startp = np.array([0.0,0.0,0.0,0.0,1.0])
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 2
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)


def feval(x):
    x = x.reshape(-1,1)
    A = np.array([[0.0,0.0,0.0,0.0,0.0],
     [2.0,1.0,1.0,1.0,3.0],
     [1.0,2.0,1.0,1.0,2.0],
     [1.0,4.0,1.0,2.0,2.0],
     [3.0,2.0,1.0,0.0,1.0],
     [0.0,2.0,1.0,0.0,1.0],
     [1.0,1.0,1.0,1.0,1.0],
     [1.0,0.0,1.0,2.0,1.0],
     [0.0,0.0,2.0,1.0,0.0],
     [1.0,1.0,2.0,0.0,0.0]]) 
    b = np.squeeze(np.array([[1.0],[5.0],[10.0],[2.0],[4.0],[3.0],[1.7],[2.5],[6.0],[3.5]]))
    X = np.tile(x.T,(10,1))
    fx = b*np.sum((X - A)**2,1)
    y= np.max(fx)
    return y
