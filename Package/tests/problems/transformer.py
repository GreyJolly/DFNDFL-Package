# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'transformer'
startp = np.array([0.8,1.5,1.2,3.0,0.8,6.0])
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
    t = np.array([0.5,0.6,0.7,0.77,0.9,1.0,1.1,1.23,1.3,1.4,1.5]).reshape(-1,1)
    f = np.zeros(11)
    v = np.array([0.0+0.0j,0.0+0.0j,0.0+0.0j,0.0+0.0j])
    w = np.array([0.0+0.0j,0.0+0.0j,0.0+0.0j,0.0+0.0j])
    for i in range(0,11):
        theta = np.pi*t[i]/2
        v[3] = np.array([1.0+0.0j])
        w[3] = np.array([10.0+0.0j])
        for k in range(2,-1,-1):
            a = np.cos(theta*x[2*k-1])             
            b = np.sin(theta*x[2*k-1])/x[2*k]            
            c = np.sin(theta*x[2*k-1])*x[2*k]            
            v[k] = np.array([a*v[k+1].real - b*w[k+1].imag + (a*v[k+1].imag + b*w[k+1].real)*1j])
            w[k] = np.array([a*w[k+1].real - c*v[k+1].imag + (a*w[k+1].imag + b*v[k+1].real)*1j])        
        f[i] = np.abs(1.0-2.0*v[0]/(w[0]+v[0]))
    y = np.max(f)
    return y
