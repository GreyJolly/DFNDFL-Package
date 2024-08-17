# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:11 2020

@author: giamp
"""

import numpy as np

name      = 'maxquad'
startp = np.zeros(10)
lb     = startp - 10.0
ub     = startp + 10.0
n      = len(lb)
nint   = 5
ncont  = n-nint
lbmix  = np.zeros(n); lbmix[:ncont] = lb[:ncont]
ubmix  = 100*np.ones(n); ubmix[:ncont] = ub[:ncont]
x_initial = 50*np.ones(n); x_initial[:ncont] = (ub[:ncont] + lb[:ncont])/2 
xmix   = np.zeros(n)

def feval(x):
    j = np.arange(1,11,dtype=float).reshape(-1,1)
    JK = np.triu(np.matmul(j,j.T),1) + np.triu(np.matmul(j,j.T),1).T
    JonK = np.triu(np.matmul(j,(j**(-1)).T),1) + np.triu(np.matmul(j,(j**(-1)).T),1).T
    A1 = np.exp(JonK)*np.cos(JK)*np.sin(1); A1 = A1 + np.diag(np.squeeze(j/10*np.sin(1) + np.sum(np.sum(np.abs(A1)))))
    A2 = np.exp(JonK)*np.cos(JK)*np.sin(2); A2 = A2 + np.diag(np.squeeze(j/10*np.sin(2) + np.sum(np.sum(np.abs(A2)))))
    A3 = np.exp(JonK)*np.cos(JK)*np.sin(3); A3 = A3 + np.diag(np.squeeze(j/10*np.sin(3) + np.sum(np.sum(np.abs(A3)))))
    A4 = np.exp(JonK)*np.cos(JK)*np.sin(4); A4 = A4 + np.diag(np.squeeze(j/10*np.sin(4) + np.sum(np.sum(np.abs(A4)))))
    A5 = np.exp(JonK)*np.cos(JK)*np.sin(5); A5 = A5 + np.diag(np.squeeze(j/10*np.sin(5) + np.sum(np.sum(np.abs(A5)))))
    
    b1 = np.exp(j/1)*np.sin(j*1);
    b2 = np.exp(j/2)*np.sin(j*2);
    b3 = np.exp(j/3)*np.sin(j*3);
    b4 = np.exp(j/4)*np.sin(j*4);
    b5 = np.exp(j/5)*np.sin(j*5);
    
    fx = np.zeros(5)
    
    fx[0] = np.matmul(x.T,np.matmul(A1,x)) - np.matmul(x.T,b1)
    fx[1] = np.matmul(x.T,np.matmul(A2,x)) - np.matmul(x.T,b2)
    fx[2] = np.matmul(x.T,np.matmul(A3,x)) - np.matmul(x.T,b3)
    fx[3] = np.matmul(x.T,np.matmul(A4,x)) - np.matmul(x.T,b4)
    fx[4] = np.matmul(x.T,np.matmul(A5,x)) - np.matmul(x.T,b5)
    
    y = np.max(fx);
    return y
