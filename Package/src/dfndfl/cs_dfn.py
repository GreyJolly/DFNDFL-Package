#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 18:05:46 2020

@author: giampo
"""

import numpy as np
import ghalton
import sobol_seq
import logging

logger = logging.getLogger("dfndfl")

class cs_dfn:

    x = []
    xint = []
    F = []
    eps = []
    funct = 0
    functpen = 0
    V = 0
    J = 0
    n = 0
    m = 0
    iprint = 0
    istop = 0
    nf_max = 0
    index_halton = 1000
    index_sobol  = 10000
    hschoice = 2
    num_fal = 0
    soglia = 1.e-3
    tol = 1.0
    print_format = '(c)| %5d | %5d | %5d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d       |   '

    flag_fail = []
    fstop = []
    xfstop= []

    sequencer = []
    Phalton   = []

    alpha_d = []
    alpha_dense = []
    alpha_diag   = []
    alpha_coord  = []
    alpha_max    = np.inf

    d_dense = []
    d1 = []
    direzioni = []
    d = []

    i_corr  = 0
    i_dense = 0
    j_dense = 0
    tipo_direzione = 0

    f = 0
    ni = 0
    fstop   = []
    xfstop  = []
    z = []
    z1= []
    z2= []
    bl = []
    bu = []
    fz1 = 0
    fz2 = 0

    def __init__(self,functpen,f,J,V,eps,xtot,nf_max,tol,bl,bu,iprint):
        
        """
        Parameters
        ----------
        n : integer
            number of CONTINUOUS variables
        nint : integer
            number of DISCRETE variables
        xtot : array
            complete vector of variables:
                first nint DISCRETE
                then  n    CONTINUOUS
        nf_max : intero
            DESCRIPTION.
        iprint : intero
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.n = V.ncont
        n = self.n
        self.nint = V.nint
        self.m = V.m[J]
        self.x = np.copy(xtot[:n])
        self.xint = np.copy(xtot[n:])
        self.z = np.copy(self.x)
        self.bl = np.copy(bl)
        self.bu = np.copy(bu)
        print_format = '(c)| %05d | %05d | %05d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d       |   '
        self.iprint = iprint
        self.tol = tol
        self.nf_max = nf_max
        self.functpen = functpen
        self.f = f
        self.V = V
        self.J = J
        self.eps = eps
        self.flag_fail = [False for i in range(n)]
        self.fstop = np.zeros(2*n)
        self.xfstop= np.zeros((n,2*n))

        self.sequencer = ghalton.Halton(n)
        self.Phalton   = self.sequencer.get(nf_max + 1000)

        self.alpha_d = np.zeros(n)
        logger.info(self.alpha_d.shape)
        self.alpha_dense = np.zeros(n)
        for i in range(n):
            self.alpha_d[i]   = np.maximum(np.double(self.soglia*10.0),np.minimum(np.double(1.0),np.abs(self.x[i])))
            self.alpha_d[i]   = (bu[i]-bl[i])/2
            for j in range(n):
                self.alpha_dense[j] += self.alpha_d[i]
            if self.iprint >= 2:
                logger.info("alphainit(%d)=%f" % (i,self.alpha_d[i]))

        self.alpha_dense /= np.double(n)
        self.alpha_diag   = np.copy(self.alpha_d)
        self.alpha_coord  = np.copy(self.alpha_d)
        self.alpha_max    = np.max(self.alpha_d)

        if n > 1:
            if self.hschoice == 1:
                self.d_dense = np.asarray(self.Phalton[self.index_halton-1], dtype = np.double)
            else:
                self.d_dense, self.index_sobol = sobol_seq.i4_sobol(n, self.index_sobol)

        self.d1 = np.zeros(n)
        self.direzioni = np.zeros((n,n))
        self.d = np.ones(n)
        for i in range(n):
            self.direzioni[i,i] = np.double(1.0)

        self.fstop   = np.zeros(2*n+1)
        self.xfstop  = np.zeros((n,2*n+1))
        for j in range(2*n+1):
            self.fstop[j] = self.f
        for i in range(n):
            for j in range(2*n+1):
                self.xfstop[i,j] = self.x[i]

    def viol_constr(self,x,J,V):
        xtot = np.concatenate((x,self.xint),axis=0)
        if J == 'a':
            gloc = V.fconstr_a(xtot)
        elif J == 'b':
            gloc = V.fconstr_b(xtot)
        elif J == 'c':
            gloc = V.fconstr_c(xtot)
        elif J == 'd':
            gloc = V.fconstr_d(xtot)
        elif J == 'e':
            gloc = V.fconstr_e(xtot)
        elif J == 'f':
            gloc = V.fconstr_f(xtot)
        else:
            gloc = V.fconstr_z(xtot)

        viol = np.maximum(0.0,np.max(gloc))

        return viol, gloc

    def func_cont(self,x,eps,J,V,CACHE):
        funct_value = self.functpen(np.concatenate((x,self.xint),axis=0),eps,J,V,CACHE)
        return funct_value

    def alpha_init(self,n,x):
        self.alpha_d = np.zeros(n)
        self.alpha_dense = np.zeros(n)
        for i in range(n):
            self.alpha_d[i]   = np.maximum(np.double(self.soglia*10.0),np.minimum(np.double(1.0),np.abs(self.x[i])))
            for j in range(n):
                self.alpha_dense[j] += self.alpha_d[i]
            if self.iprint >= 2:
                logger.info("alphainit(%d)=%f" % (i,self.alpha_d[i]))

        self.alpha_dense /= np.double(n)
        self.alpha_diag   = np.copy(self.alpha_d)
        self.alpha_coord  = np.copy(self.alpha_d)
        self.alpha_max    = np.max(self.alpha_d)

    def step(self,alg,xint,f,eps,V,CACHE,nf):
        self.f    = f
        self.xint = np.copy(xint)
        self.eps  = np.copy(eps)
        cambio_eps = False
        while True:
            if self.n > 1:
                self.alpha_max = np.max([np.max(self.alpha_coord),np.max(self.alpha_diag),np.max(self.alpha_dense)])
            else:
                self.alpha_max = np.max(self.alpha_d)

            istop = stop(self.alpha_d,self.alpha_max,nf,self.ni,self.tol,self.nf_max)

            self.alpha_max = np.max(self.alpha_d)

            if istop >= 1:
                self.eps = np.copy(eps)
                #self.f   = f
                if self.iprint >= 2:
                    logger.info(self.alpha_d)
                    logger.info(self.alpha_max)
                return nf, cambio_eps
                #return nf

            if self.i_corr == 0:
                self.dconv = np.zeros(self.n)
                for i in range(self.n):
                    self.dconv += -self.direzioni[:,i]

            if self.iprint >= 2:
                logger.info("----------------------------------------------")
            #if self.iprint >= 0:
            #	logger.info(" ni=%4d  nf=%5d   f=%12.5e   alphamax=%12.5e" % (self.ni,self.nf,self.f,self.alpha_max))
            if self.iprint >= 2:
                #if (self.m > 0):
                #	logger.info(self.print_format %(self.ni, self.nf, 0, self.f, self.f, 0.0, 0.0, self.alpha_max, 0))
                #else:
                viol, g = self.viol_constr(self.x,self.J,self.V)
                logger.info(self.print_format %(self.ni, nf, CACHE.hits, self.f, self.f, np.sum(np.where(g<0,0,g)), self.alpha_max, self.i_corr))
            if self.iprint >= 2:
                for i in range(self.n):
                    logger.info(" x(",i,")=",self.x[i])

            self.d = np.copy(self.direzioni[:,self.i_corr])

            if self.tipo_direzione == 0:
                self.alpha, self.fz, self.z1, self.fz1,self.z2, self.fz2, self.i_corr_fall, nf = linesearchbox_cont(
                                    self.func_cont,self.x,self.f,self.d,
                                    self.alpha_d,self.z,self.i_corr,self.alpha_max,
                                    self.bl,self.bu,nf,self.iprint,eps,self.J,self.V,CACHE)

                if self.alpha >= 1.e-12:
                    self.x[self.i_corr] = self.x[self.i_corr]+self.alpha*self.d[self.i_corr]
            else:

                self.alpha, self.alphatilde, self.fz, self.d, nf = linesearchbox_dense(self.func_cont,self.x,self.f,self.d,
                                                                              self.alpha_d[self.i_corr],self.alpha_max,
                                                                              self.bl,self.bu,nf,self.iprint,
                                                                              eps,self.J,self.V,CACHE)
                self.alpha_d[self.i_corr]      = self.alphatilde

                if self.alpha >= 1.e-12:
                    self.x = np.maximum(self.bl,np.minimum(self.bu,self.x+self.alpha*self.d))

            self.direzioni[:,self.i_corr] = np.copy(self.d)

            if self.alpha >= 1.e-12:
                self.flag_fail[self.i_corr] = False
                self.f = self.fz
                self.num_fal = 0
            else:
                self.flag_fail[self.i_corr] = True
                if self.i_corr_fall == 0:
                    self.fstop[self.i_corr]   = self.fz1
                    self.fstop[2*self.i_corr] = self.fz2
                    for j in range(self.n):
                        self.xfstop[j,self.i_corr]   = self.z1[j]
                        self.xfstop[j,2*self.i_corr] = self.z2[j]
                    self.num_fal += 1

            self.ni += 1
            self.z = np.copy(self.x)

            if self.iprint >= 1:
                viol, g = self.viol_constr(self.x,self.J,self.V)
                logger.info(self.print_format %(self.ni, nf, CACHE.hits, self.f, self.f, np.sum(np.where(g<0,0,g)), self.alpha_max, self.i_corr))

            if self.i_corr < self.n-1:
                self.i_corr += 1
            else:
                if alg == 'DFN_DFL' and np.max(self.alpha_d) <= self.soglia and self.n > 1:
                    if self.tipo_direzione == 0:
                        fmin = self.fstop[0]
                        fmax = self.fstop[0]
                        imin = 0
                        imax = 0
                        doldalphamin = self.alpha_d[0]
                        doldalphamax = self.alpha_d[0]
                        iminalpha = 0
                        imaxalpha = 0
                        for i in range(1,self.n):
                            if self.alpha_d[i] < doldalphamin:
                                doldalphamin = self.alpha_d[i]
                                iminalpha = i
                            if self.alpha_d[i] > doldalphamax:
                                doldalphamax = self.alpha_d[i]
                                imaxalpha = i
                        rapalpha = 3.0
                        if doldalphamax/doldalphamin > rapalpha:
                            for i in range(self.n):
                                self.d1[i] = self.dconv[i]
                            self.dnr = np.sqrt(np.double(self.n))
                        else:
                            for i in range(1,2*self.n):
                                    if self.fstop[i] < fmin:
                                        fmin = self.fstop[i]
                                        imin = i
                                    if self.fstop[i] > fmax:
                                        fmax = self.fstop[i]
                                        imax = i

                            self.dnr = np.double(0.0)
                            doldalphamedio = (doldalphamin+doldalphamax)/2.0
                            for i in range(self.n):
                                self.d1[i] = self.xfstop[i,imin]-self.xfstop[i,imax]
                                self.dnr += self.d1[i]*self.d1[i]
                            self.dnr = np.sqrt(self.dnr)
                            if self.dnr <= 1.e-24:
                                for i in range(self.n):
                                    self.d1[i] = self.dconv[i]
                                self.dnr = np.sqrt(np.double(self.n))

                        self.direzioni = gen_base(self.d1)
                        self.direzioni = gram_schmidt(self.direzioni)
                        self.tipo_direzione = 1

                        self.alpha_coord = np.copy(self.alpha_d)
                        if doldalphamax/doldalphamin > rapalpha:
                            self.alpha_d = np.copy(self.alpha_diag)
                        else:
                            self.dnr = np.sum(self.alpha_d)/np.double(self.n)
                            for i in range(self.n):
                                self.alpha_d[i] = self.dnr

                        if self.iprint >= 2:
                            logger.info("END DIR. COORDINATES")

                    elif self.tipo_direzione == 1:
                        self.direzioni = gen_base(self.d_dense)
                        self.direzioni = gram_schmidt(self.direzioni)
                        self.index_halton += 2*self.n
                        if self.hschoice == 1:
                            self.d_dense = np.asarray(self.Phalton[self.index_halton-1], dtype = np.double)
                        else:
                            self.d_dense, self.index_sobol = sobol_seq.i4_sobol(self.n, self.index_sobol)

                        self.tipo_direzione = 2
                        self.alpha_diag = np.copy(self.alpha_d)

                        self.dnr = np.sum(self.alpha_d)/np.double(self.n)
                        for i in range(self.n):
                            self.alpha_d[i] = np.double(10.0)*self.dnr

                        if self.iprint >= 2:
                            logger.info("END DIR. N+1")

                    elif self.tipo_direzione == 2:

                        self.direzioni = np.zeros((self.n,self.n))
                        for i in range(self.n):
                            self.direzioni[i,i] = np.double(1.0)

                        self.tipo_direzione = 0

                        self.alpha_dense = np.copy(self.alpha_d)
                        self.alpha_d = np.copy(self.alpha_coord)

                        if self.iprint >= 2:
                            logger.info("END DIR. DENSE")

                    self.i_corr = 0
                    break

                self.i_corr = 0
                break

        #####################
        #check eps
        #####################
        cambio_eps = False
        if self.m > 0:
            viol, constr = self.viol_constr(self.x,self.J,self.V)
            #logger.info(viol,constr)
            if viol > 0.0:
                maxeps = np.max(eps)
                for i in range(self.m):
                    if(eps[i]*constr[i] > 1.0*np.max(self.alpha_d)):
                        eps[i]=1.e-2*eps[i]
                        if self.iprint >= 0:
                            logger.info('**************************************')
                            logger.info('*********** updating eps(',i,')=',eps[i],' *************')
                            logger.info('**************************************')
                        cambio_eps = True

                if cambio_eps:
                    self.f = self.func_cont(self.x,eps,self.J,self.V,CACHE)
                    #nf += 1

                    for i in range(self.n):
                        self.alpha_d[i]   = np.maximum(np.double(self.soglia*100.0),self.alpha_d[i])
                        self.alpha_d[i]   = (self.bu[i]-self.bl[i])/2

        self.eps = np.copy(eps)
        #self.f   = f

        return nf, cambio_eps
        #return nf

def gen_base(d):
    n = len(d)
    ind = np.argmax(np.abs(d))
    H = np.zeros((n,n))
    H[:,0] = d
    #logger.info(ind)
    for i in range(1,ind+1):
        H[i-1,i] = np.double(1.0)
    for i in range(ind+1,n):
        H[i,i] = np.double(1.0)

    return H

def gram_schmidt(H):
    (n,n) = H.shape
    for i in range(1,n):
        proj = np.double(0.0)
        for j in range(i):
            proj += (np.dot(H[:,i],H[:,j])/np.dot(H[:,j],H[:,j]))*H[:,j]
        H[:,i] -= proj

    for i in range(n):
        H[:,i] = H[:,i]/np.linalg.norm(H[:,i])

    return H

def stop(alpha_d,alpha_max,nf,ni,tol,nf_max):
    istop = 0
    if alpha_max <= tol:
        istop = 1

    if nf > nf_max:
        istop = 2

    if ni > nf_max:
        istop = 3

    return istop

def linesearchbox_cont(funct,x,f,d,alpha_d,z,j,alpha_max,bl,bu,nf,iprint,eps,J,V,CACHE):
    gamma = np.double(1.e-6)
    delta = np.double(0.5)
    delta1= np.double(0.5)
    ifront= 0
    i_corr_fall = 0
    n = len(x)

    z1 = np.zeros(n)
    z2 = np.zeros(n)
    fz1 = 0.0
    fz2 = 0.0

    if iprint >= 2:
        logger.info("continuous variable j =",j,"    d(j) =",d[j]," alpha=",alpha_d[j])

    if np.abs(alpha_d[j]) <= 1.e-3*np.minimum(1.0,alpha_max):
        alpha = np.double(0.0)
        fz = 0.0
        if iprint >= 2:
            logger.info(" small alpha")
            logger.info(" alpha_d(j)=",alpha_d[j],"    alphamax=",alpha_max)
        return alpha,fz, z1, fz1, z2, fz2, i_corr_fall, nf

    ifront = 0

    for ielle in [1,2]:
        if d[j] > 0.0:
            if alpha_d[j] - (bu[j]-x[j]) < -1.e-6:
                alpha = np.maximum(1.e-24,alpha_d[j])
            else:
                alpha = bu[j]-x[j]
                ifront = 1
                if iprint >= 2:
                    logger.info(" point expan. on edge *")
        else:
            if alpha_d[j] - (x[j]-bl[j]) < -1.e-6:
                alpha = np.maximum(1.e-24,alpha_d[j])
            else:
                alpha = x[j]-bl[j]
                ifront = 1
                if iprint >= 2:
                    logger.info(" point expan. on edge *")

        if np.abs(alpha) <= 1.e-3*np.minimum(1.0,alpha_max):
            d[j] = -d[j]
            i_corr_fall += 1
            ifront = 0
            if iprint >= 2:
                logger.info(" opposite direction for small alpha")
                logger.info(" j =",j,"    d(j) =",d[j])
                logger.info(' alpha=',alpha,'    alphamax=',alpha_max)
            alpha = np.double(0.0)
        else:
            alphaex = alpha
            z[j] = x[j] + alpha*d[j]
            fz = funct(z,eps,J,V,CACHE)
            nf += 1

            if ielle == 1:
                z1 = np.copy(z)
                fz1 = fz
            else:
                z2 = np.copy(z)
                fz2 = fz

            if iprint >= 2:
                logger.info(' fz =',fz,'   alpha =',alpha)
            if iprint >= 3:
                for i in range(n):
                    logger.info(' z(',i,')=',z[i])

            fpar = f - gamma * alpha**2
            if fz < fpar:
                while True:
                    if ifront == 1:
                        if iprint >= 2:
                            logger.info(' accept point on edge fz =',fz,'   alpha =',alpha)

                        alpha_d[j] = delta*alpha

                        return alpha, fz, z1, fz1, z2, fz2, i_corr_fall, nf

                    if d[j] > 0.0:
                        if alpha/delta1 - (bu[j]-x[j]) < -1.e-6:
                            alphaex = alpha/delta1
                        else:
                            alphaex = bu[j]-x[j]
                            ifront = 1
                            if iprint >= 2:
                                logger.info(' point expan. on edge')
                    else:
                        if alpha/delta1 - (x[j]-bl[j]) < -1.e-6:
                            alphaex = alpha/delta1
                        else:
                            alphaex = x[j]-bl[j]
                            ifront = 1
                            if iprint >= 2:
                                logger.info(' point expan. on edge')

                    z[j] = x[j] + alphaex*d[j]
                    fzdelta = funct(z,eps,J,V,CACHE)
                    nf += 1

                    if iprint >= 2:
                        logger.info(' fzex=',fzdelta,'  alphaex=',alphaex)
                    if iprint >= 3:
                        for i in range(n):
                            logger.info(' z(',i,')=',z[i])

                    fpar = f - gamma * alphaex**2
                    if fzdelta < fpar:
                        fz = fzdelta
                        alpha = alphaex
                    else:
                        alpha_d[j] = delta*alpha
                        if iprint >= 2:
                            logger.info(' accept point fz =',fz,'   alpha =',alpha)
                        return alpha, fz, z1, fz1, z2, fz2, i_corr_fall, nf

            else:
                d[j] = -d[j]
                ifront = 0
                if iprint >= 2:
                    logger.info(' opposite direction')
                    logger.info(' j =',j,'    d(j) =',d[j])

    if not i_corr_fall==2:
        alpha_d[j] = delta*alpha_d[j]

    alpha = 0.0

    if iprint >= 2:
        logger.info(' direction fail')

    return alpha, fz, z1, fz1, z2, fz2, i_corr_fall, nf


def linesearchbox_dense(funct,x,f,d,alpha_d,alpha_max,bl,bu,nf,iprint,eps,J,V,CACHE):
    gamma = np.double(1.e-6)
    delta = np.double(0.5)
    delta1= np.double(0.5)
    ifront= 0

    if iprint >= 2:
        logger.info("halton direction, alpha=",alpha_d)

    for ielle in [1, 2]:
        alpha   = alpha_d
        alphaex = alpha
        z      = x + alpha*d
        z      = np.maximum(bl,np.minimum(bu,z))
        fz     = funct(z,eps,J,V,CACHE)
        nf    += 1

        if iprint >= 2:
            logger.info(" fz =",fz,"   alpha =",alpha)
        if iprint >= 3:
            for i in range(n):
                logger.info(" z(",i,")=",z[i])

        fpar = f - gamma*alpha**2
        if fz < fpar:
            while True:
                alphaex = alpha/delta1
                z      = x + alphaex*d
                z      = np.maximum(bl,np.minimum(bu,z))
                fzdelta= funct(z,eps,J,V,CACHE)
                nf    += 1

                if iprint >= 2:
                    logger.info(" fzex=",fzdelta,"   alphaex=",alphaex)
                if iprint >= 3:
                    for i in range(n):
                        logger.info(" z(",i,")=",z[i])

                fpar = f - gamma*alphaex**2
                if fzdelta < fpar:
                    fz   = fzdelta
                    alpha = alphaex
                else:
                    alpha_d = alpha
                    if iprint >= 2:
                        logger.info(" dense: accept point fz =",fz,"   alpha =",alpha)
                    return alpha, alpha_d, fz, d, nf
        else:
            d      = -d
            ifront =  0

            if iprint >= 2:
                logger.info("dense:  opposite direction")

    alpha_d = delta*alpha_d
    alpha   = np.double(0.0)

    if iprint >= 2:
        logger.info("dense: direction fail")

    return alpha, alpha_d, fz, d, nf
