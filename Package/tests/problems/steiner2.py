# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'steiner 2'
aa = np.array([[0.0,2.0],[2.0,3.0],[3.0,-1.0],[4.0,-0.5],[5.0,2.0],[6.0,2.0]])
m=6
xbar = np.zeros(2*m)
xbar[0] = (aa[0,0]+aa[1,0])/3
j = np.arange(2,m)
xbar[j-1] = (xbar[j-2] + aa[j-1,0] + aa[j,0])/3
xbar[m-1] = (xbar[m-2] + aa[m-1,0] + 5.5)/3
xbar[m] = (aa[0,1] + aa[1,1])/3
xbar[j+m-1] = (xbar[j-2+m] + aa[j-1,1] + aa[j,1])/3
xbar[2*m-1] = (xbar[2*m-2] + aa[m-1,1] - 1.0)/3
startp = xbar
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 6
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    x = x.reshape(-1,1)
    m = 6
    a = np.array([[0.0,2.0],
         [2.0,3.0],
         [3.0,-1.0],
         [4.0,-0.5],
         [5.0,2.0],
         [6.0,2.0]])
    p = np.array([[2.0],[1.0],[1.0],[5.0],[1.0],[1.0]])
    ptilde = np.array([[1.0],[1.0],[2.0],[3.0],[2.0]])   
    jm = np.arange(1,m+1)
    jm1= np.arange(1,m)
    y = np.sqrt(x[0]**2 + x[m]**2) + np.sqrt((5.5 - x[m-1])**2 + (-1.0 - x[2*m-1])**2) + np.sum(p*np.sqrt((a[jm-1,0].reshape(-1,1) - x[jm-1])**2 + (a[jm-1,1].reshape(-1,1) - x[jm+m-1])**2),axis=0) + np.sum(ptilde*np.sqrt((x[jm1-1] - x[jm1])**2 + (x[jm1+m-1] - x[jm1+m])**2),axis=0)
    return np.asscalar(y)
