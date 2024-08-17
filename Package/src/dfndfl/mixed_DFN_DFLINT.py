
##########################################################################
# VARIABLES MEANING
# The first NINT variables are discrete and the following N-NINT are continues
# THe i-th discrete variable has the values:
#	      LBMIX[i], LBMIX[i]+1, ..., UBMIX[i]
##########################################################################

import numpy as np
import ghalton
from . import cs_dfn
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def box_DFN_DFL(alg, M, J, V, CACHE, max_fun, outlev):
    #
    # Written by T. Giovannelli, G. Liuzzi, S. Lucidi, F. Rinaldi, 2020.
    #

    nf        = 0
    lb        = V.lbmix
    ub        = V.ubmix
    x_initial = np.copy(V.x_initial)
    ncont     = V.ncont
    nint      = V.nint
    x_initial[ncont:] = np.round(V.x_initial[ncont:])

    ubint = np.copy(ub)
    lbint = np.copy(lb)
    ubint[ncont:] = np.round(ubint[ncont:])
    lbint[ncont:] = np.round(lbint[ncont:])

    m = V.m[J]
    if M < 1:
        M = 1

    #if ncont > 0:
    eps         = 0.1*np.ones((m,1))

    if(nint <= 1):
        logger.error('Number of integer variables must be > 1.')
        x       = x_initial
        f       = np.inf
        stopfl  = 99
        Dused   = []
        return x,f, stopfl, Dused, nf

    if (np.sum(abs(V.ubmix[ncont:]-ubint[ncont:])) != 0) or (np.sum(abs(V.lbmix[ncont:]-lbint[ncont:])) != 0):
        logger.error('Upper and/or lower bound on some variable is NOT integer.')
        x       = x_initial
        f       = np.inf
        stopfl  = 99
        Dused   = []
        return x,f, stopfl, Dused, nf

    iter        = 0                 # iteration counter
    n 			= nint              # dimension of the problem
    stop        = 0                 # stopping condition
    W           = np.empty(M)       # memory of function value for NM linesearch
    xW          = np.empty((n+ncont,M))

    sequencer   = ghalton.Halton(n)       ## (input) space size
    Phalton     = sequencer.get(1000000)   ## (input) number of points to generate in the space
    ihalton     = 7

    #logger.info(Phalton[0])
    #logger.info(CS.Phalton[0])


    eta         = 1.5
    allones     = 0
    tol         = 1.e-4

    rho         = 1.0
    xi          = 1.0

    chk_feas = (x_initial >= V.lbmix) & (x_initial <= V.ubmix)
    if min(chk_feas) == 0:
        logger.error('Initial point does not satisfy the bound constraints!')
        x       = x_initial
        f       = np.inf
        stopfl  = 99
        Dused   = []
        return x,f, stopfl, Dused, nf

    # D           denotes the set of search direction (one per each column)
    # alpha_tilde is a row vector of stepsizes (one per each direction in D)

    x  = x_initial
    f  = functpen(x,eps,J,V,CACHE)
    nf += 1

    CS = cs_dfn.cs_dfn(functpen,f,J,V,eps,x_initial,max_fun,tol,lbint[:ncont],ubint[:ncont],outlev)

    if m > 0:
        g = CACHE.xf[0,(n+ncont+1):]
        rho = max(1.0,np.sum(np.where(g<0,0,g)))
#        logger.info(g.shape)
#        logger.info(eps.shape)
        eps[(g > 0)<1.0] = 1e-3
    else:
        g = np.array([-1.0])

    W[0] = f
    xW[:,0] = x
    bestf = f
    bestx = x

    D = np.concatenate((np.zeros((ncont,nint)),np.identity(nint)),axis=0)
    successes = np.zeros(nint) # components of this vector counts number of successes each direction has had
                            # where "success" means that alpha > 0 is returned
                            # by LS along the direction itself


    alpha_tilde   = np.round((ubint[ncont:]-lbint[ncont:])/2.0)

    old_maxalpha  = np.inf
    ndir          = np.shape(D)[1]
    
    # logger.info("%04d-%02d-%02d" % (year, month, day), file=log)

    if outlev > 0:
        print_format = '(d)| %5d | %5d | %5d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d/%5d |   '
        logger.debug('   |  iter |    nf | cache |        fpen     |        f_ref    |         viol    |    max_alpha    |  ndir |')
        logger.info('   |  iter |    nf | cache |        fpen     |        f_ref    |         viol    |    max_alpha    |  ndir |')
    else:
        logger.debug('   fun.evals =      ')
#        logger.info('   fun.evals =      ')

    while stop != 1:

        iter += 1

#########################################################
#########################################################
#       explore continuous variables
#########################################################
#########################################################
        if ncont > 0:
            CS.ni = iter
            nf, change_eps = CS.step(alg,x[ncont:],f,eps,V,CACHE,nf)
            iter = CS.ni
            y = np.concatenate((CS.x,CS.xint),axis=0)
            fy = CS.f

            eps = np.copy(CS.eps)
            if change_eps:
                #W       = np.nan*np.ones((1,M))
                #xW      = np.nan*np.ones((n+ncont,M))
                W       = -np.inf*np.ones(M)       # memory of function value for NM linesearch
                xW      = np.nan*np.empty((n+ncont,M))
                W[0]    = fy
                xW[:,0] = y
                if outlev >= 2:
                    logger.info('(0) W = ',W)
                bestf   = fy
                bestx   = y
                logger.info('after cs-dfn fy = ',fy,max(W))
            else:
                if m == 0:
                    W = np.roll(W,1)
                    xW = np.roll(xW,1,axis=1)
                    W[0] = fy
                    xW[:,0] = y
                if(fy < bestf):
                    bestf = fy
                    bestx = y

            if (m > 0) and (np.linalg.norm(x-y) > 1.e-9):
                W = np.roll(W,1)
                xW = np.roll(xW,1,axis=1)
                W[0] = fy
                xW[:,0] = y
                if outlev >= 2:
                    logger.info('(6) W = ',W)
                if(fy < bestf):
                    bestf = fy
                    bestx = y
        else:
            y  = x
            fy = f

        for idir in range(ndir):

#########################################################
#########################################################
            if nf >= max_fun:
                return x, f, stop, D, nf
#########################################################
#########################################################


            d = D[:,idir]

            if iter == 1:
                f_ref = W[0]
            else:
                f_ref = max(W)

            if outlev >= 2:
                logger.info('discrete_search along d',f_ref,max(W),W)
            alpha, x_trial, f_trial, nf = nm_discrete_linesearch(y,d,alpha_tilde[idir],lb,ub,f_ref,eps,xi,m,ncont,nf,J,V,CACHE,outlev)

#########################################################
#########################################################
            if nf >= max_fun:
                return x, f, stop, D, nf
#########################################################
#########################################################
            if outlev >= 2:
                logger.info('discrete_search returned alpha = ',alpha)

            if alpha <= 0:
                d = -d

                if outlev >= 2:
                    logger.info('discrete_search along -d')

                alpha, x_trial, f_trial, nf = nm_discrete_linesearch(y,d,alpha_tilde[idir],lb,ub,f_ref,eps,xi,m,ncont,nf,J,V,CACHE,outlev)
                if outlev >= 2:
                    logger.info('discrete_search returned alpha = ',alpha)
                if alpha > 0:
                    successes[idir] = successes[idir]+1
                    if ncont > 0:
                        CS.alpha_init(ncont,x_trial)
                    if allones >= 1:
                        allones = 0
                    D[:,idir] = d
                    y  = x_trial
                    fy = f_trial

                    alpha_tilde[idir] = alpha
                    W = np.roll(W,1)
                    xW = np.roll(xW,1,axis=1)
                    W[0] = fy
                    xW[:,0] = y
                    if outlev >= 2:
                        logger.info('(1) W = ',W)
                    if(fy < bestf):
                        bestf = fy
                        bestx = y
                else:
                    alpha_tilde[idir] = max(1,np.floor(alpha_tilde[idir]/2))
            else:
                successes[idir] = successes[idir]+1
                if ncont > 0:
                    CS.alpha_init(ncont,x_trial)
                if allones >= 1:
                    allones = 0

                y  = x_trial
                fy = f_trial
                alpha_tilde[idir] = alpha
                W = np.roll(W,1)
                xW = np.roll(xW,1,axis=1)
                W[0] = fy
                xW[:,0] = y
                if outlev >= 2:
                    logger.info('(2) W = ',W)
                if(fy < bestf):
                    bestf = fy
                    bestx = y

            if outlev > 0:
                viol, g = viol_constr(y,J,V)
                logger.debug(print_format %(iter, nf, CACHE.hits, fy, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), idir,ndir))
                logger.info(print_format %(iter, nf, CACHE.hits, fy, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), idir,ndir))

            if m > 0:
                if (allones >= 1):
                   break
            else:
                if (allones > 1):
                   break

        if(m > 0):
            sxf = np.shape(CACHE.xf)[0]
            diff = np.square(CACHE.xf[:,:(n+ncont)]-np.tile(y,(sxf,1)))
            mn = min(np.sum(diff,axis=1))
            ind = np.argmin(np.sum(diff,axis=1))
            g = CACHE.xf[ind,(n+ncont+1):]

        if alg == 'DFN_DFL':
            if (np.linalg.norm(y-x) <= 1e-14) and (max(alpha_tilde) == 1) and (old_maxalpha == 1):

                xi=xi/2

                allones=allones+1

                iexit = 0

                while iexit == 0:
                    # enrich set D
                    D, successes, alpha_tilde, iexit, ihalton = generate_dirs(ncont,nint,D,successes,alpha_tilde,eta,0,Phalton,ihalton)
                    if iexit == 0:
                        eta = eta + 0.5
                    if eta >= 0.5*(np.linalg.norm(ub -lb)/2):
                        #stop execution
                        if (bestf < fy):
                            y  = bestx
                            fy = bestf
                        else:
                            stopfl = 1
                            stop   = 1
                            Dused  = D
                        iexit = 1

                ndir = np.shape(D)[1]

        viol, g = viol_constr(y,J,V)
        if (m > 0) and (viol > 0.0):
            #check on the penalty parameters eps
            ind_change = np.where(g<0,0,g) > rho
            nchg = np.shape(ind_change)[0]
            eps_changed = 0
            for i in range(nchg):
                if (ind_change[i] == True) and (eps[i] > 1e-10):
                    eps[i] = eps[i]/2
                    allones = 0
                    eps_changed = 1

            if eps_changed == 1:
                sxf= np.shape(CACHE.xf)[0]
                diff=np.square(CACHE.xf[:,:(n+ncont)]-np.tile(y,(sxf,1)))
                mn=min(np.sum(diff,axis=1))
                ind=np.argmin(np.sum(diff,axis=1))
                if mn <= 1e-16:
                    fval = CACHE.xf[ind,n+ncont]
                    gval = CACHE.xf[ind,n+ncont+1:]
                    m = np.shape(gval)[0]
                    fy = fval + np.sum(np.divide(np.where(gval<0,0,gval),eps))
                    CACHE.hits += 1
                else:
                    fy = functpen(y,eps,J,V,CACHE)
                    nf += 1

                W       = np.nan*np.ones(M)
                xW      = np.nan*np.ones((n+ncont,M))
                W[0]    = fy
                xW[:,0] = y
                if outlev >= 2:
                    logger.info('(3) W = ',W)
                bestf   = fy
                bestx   = y

        if m > 0:
            rho = max(1e-8,rho*0.5)
        else:
            rho = rho/2.0

        #if (np.linalg.norm(y-x) <= 1e-14) and (ndir >=5000):
        if (np.linalg.norm(y-x) <= 1e-14) and (((alg == 'DFN_DFL') and (ndir >=5000)) or (alg == 'DFL' and (max(alpha_tilde) == 1) and (old_maxalpha == 1) and (CS.alpha_max <= tol))):
        #if (np.linalg.norm(y-x) <= 1e-14) and (((alg == 'DFN_DFL') and (ndir >=5000)) or (alg == 'DFL')):
            stopfl = 1
            stop   = 1
            Dused  = D
            x      = bestx
            f      = bestf

        x = y
        f = fy

        old_maxalpha = max(alpha_tilde)

        if outlev > 0:
            viol, g = viol_constr(x,J,V)
            logger.debug(print_format %(iter, nf, CACHE.hits, f, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), ndir,ndir))
            logger.info(print_format %(iter, nf, CACHE.hits, f, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), ndir,ndir))

        if nf >= max_fun:
            stopfl = 2
            stop   = 1
            Dused  = D
            x      = bestx
            f      = bestf

    if outlev > 0:
            viol, g = viol_constr(x,J,V)
            logger.debug(print_format %(iter, nf, CACHE.hits, f, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), ndir,ndir))
            logger.info(print_format %(iter, nf, CACHE.hits, f, f_ref, np.sum(np.where(g<0,0,g)), max(alpha_tilde), ndir,ndir))
    if outlev == 0:
        logger.debug('\n')

    return x, f, stopfl, Dused, nf
##########################################################################
# END OF CODE box_DFL
##########################################################################

def viol_constr(xtot,J,V):
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

##########################################################################
# functpen
##########################################################################
def functpen(xint,eps,J,V,CACHE):
    floc = V.func_f(xint)
    if J == 'a':
        gloc = V.fconstr_a(xint)
    elif J == 'b':
        gloc = V.fconstr_b(xint)
    elif J == 'c':
        gloc = V.fconstr_c(xint)
    elif J == 'd':
        gloc = V.fconstr_d(xint)
    elif J == 'e':
        gloc = V.fconstr_e(xint)
    elif J == 'f':
        gloc = V.fconstr_f(xint)
    else:
        gloc = V.fconstr_z(xint)

    nx = np.shape(xint)[0]
    ng = np.shape(gloc)[0]

    CACHE.xf[CACHE.nnf,:(nx)] = xint
    CACHE.xf[CACHE.nnf, (nx)] = floc
    if V.m[J] > 0:
        CACHE.xf[CACHE.nnf,(nx+1):] = gloc
        #logger.info('FUNCTPEN:   g1 = ',np.maximum(0.0,np.max(gloc)),J)
        #logger.info('FUNCTPEN:    g = ',gloc)
        #logger.info('FUNCTPEN: viol = ',np.sum(np.where(gloc<0,0,gloc)))
        fpen = floc + np.sum(np.divide(np.where(gloc<0,0,gloc),eps))
        if np.sum(np.where(gloc<0,0,gloc)) > 1.e-4:
            if CACHE.nnf < CACHE.maxf:
                CACHE.F.append(np.inf)
                CACHE.nnf += 1
        else:
            if CACHE.nnf < CACHE.maxf:
                CACHE.F.append(floc)
                CACHE.nnf += 1
    else:
        fpen = floc
        if CACHE.nnf < CACHE.maxf:
            CACHE.F.append(floc)
            CACHE.nnf += 1

    return fpen
##########################################################################
# END OF CODE functpen
##########################################################################


def nm_discrete_linesearch(y,d,alpha_tilde,lb,ub,f_ref,eps,xi,m,ncont,nf,J,V,CACHE,outlev):

    #
    # Function nm_discrete_linesearch
    #
    # Purpose:
    #
    # This function performs a nonmonotone discrete linesearch
    # along a given direction d (d \in Z^n)
    #
    # Inputs:
    #
    # y            : starting point for the linesearch
    #
    # d            : search direction
    #
    # alpha_tilde  : starting stepsize
    #
    # lb, ub       : lower and upper bounds
    #
    # f_ref        : reference o.f. value
    #
    # Output:
    #
    #
    # alpha        : 1) alpha > 0 if linesearch finds a point guaranteeing
    #                simple decrease: f(y+alpha d)<f_ref
    #                2) alpha = 0 failure
    #
    # x            : best point found in the linesearch
    #
    # f            : o.f. value related to x
    #

    # calculate dimension of the problem
    n = len(d)

    # initialize vector alpha_max
    alpha_max = np.inf * np.ones(n)

    # calculate max alpha
    indices = ( d > 0 )

    alpha_max[indices]=np.divide(ub[indices] - y[indices],d[indices])


    indices = ( d < 0 )

    alpha_max[indices]=np.divide(lb[indices] - y[indices],d[indices])

    #compute starting alpha
    alpha_bar  = np.floor( min(alpha_max) )
    alpha_init = min(alpha_tilde, alpha_bar)

    if outlev >= 2:
        logger.info('nm_discrete_search: alpha_init = ',alpha_init)

    #Build first point for starting linesearch
    if (alpha_init > 0):
        y_trial = y + alpha_init * d
        sxf=np.shape(CACHE.xf)[0]
        diff=np.square(CACHE.xf[:,:(n)]-np.tile(y_trial,(sxf,1)))
        mn=min(np.sum(diff,axis=1))
        ind=np.argmin(np.sum(diff,axis=1))
        #diff
        #keyboard
        if (mn<=1e-16):
            fval = CACHE.xf[ind,n]
            gval = CACHE.xf[ind,n+1:]
            if V.m[J] > 0:
                f_trial = fval + np.sum(np.divide(np.where(gval<0,0,gval),eps))
            else:
                f_trial = fval
            CACHE.hits += 1
            if outlev >= 2:
                logger.info('nm_discrete_search: f_trial computed by cache')
        else:
            f_trial = functpen(y_trial,eps,J,V,CACHE)
            nf += 1
            if outlev >= 2:
                logger.info('nm_discrete_search: f_trial computed by functpen',f_trial)
                logger.info('                  : fob  = ',CACHE.xf[CACHE.nnf-1,n])
                gval = CACHE.xf[CACHE.nnf-1,n+1:]
                logger.info('                  : viol = ',np.sum(np.where(gval<0,0,gval)))
                logger.info('                  : gval = ',gval)
                logger.info('                  : eps  = ',eps)
    else:
        f_trial = np.inf


    # cicle for updating alpha
    if outlev >= 2:
        logger.info('nm_discrete_search: ftrial = ',f_trial)
        logger.info('                  : f_ref  = ',f_ref)
        logger.info('                  : alpha_init = ',alpha_init)

    if (alpha_init > 0) and (f_trial <= f_ref - xi):

        # initialize alpha and best point
        if outlev >= 2:
            logger.info('nm_discrete_search: f_trial <= f_ref - xi')
        alpha=alpha_init
        x = y_trial
        f = f_trial

        #calculate trial point
        if alpha < alpha_bar:
            y_trial = y + min(alpha_bar,2*alpha)* d
            sxf = np.shape(CACHE.xf)[0]
            diff = np.square(CACHE.xf[:,:(n)]-np.tile(y_trial,(sxf,1)))
            mn = min(np.sum(diff,axis=1))
            ind = np.argmin(np.sum(diff,axis=1))
            #diff
            #keyboard
            if (mn<=1e-16):
                fval = CACHE.xf[ind,n]
                gval = CACHE.xf[ind,n+1:]
                if V.m[J] > 0:
                    f_trial = fval + np.sum(np.divide(np.where(gval<0,0,gval),eps))
                else:
                    f_trial = fval
                CACHE.hits += 1
                if outlev >= 2:
                    logger.info('nm_discrete_search: f_trial computed by cache')
            else:
                f_trial = functpen(y_trial,eps,J,V,CACHE)
                nf += 1
                if outlev >= 2:
                    logger.info('nm_discrete_search: f_trial computed by functpen',f_trial)
                    logger.info('                  : fob  = ',CACHE.xf[CACHE.nnf-1,n])
                    gval = CACHE.xf[CACHE.nnf-1,n+1:]
                    logger.info('                  : viol = ',np.sum(np.where(gval<0,0,gval)))
                    logger.info('                  : gval = ',gval)
                    logger.info('                  : eps  = ',eps)
        else:
            f_trial = np.inf


        # expansion step (increase stepsize)
        while (alpha<alpha_bar) and (f_trial <= f_ref - xi):

            if outlev >= 2:
                logger.info('nm_discrete_search: f_trial <= fref - xi (while loop)')
            # alpha calulation and best point updatingd
            alpha=min(alpha_bar, 2*alpha)

            # best point updating
            x = y_trial
            f = f_trial

            #next point to be tested
            if(alpha < alpha_bar):
                y_trial = y + min(alpha_bar, 2* alpha) * d
                sxf = np.shape(CACHE.xf)[0]
                diff = np.square(CACHE.xf[:,:(n)]-np.tile(y_trial,(sxf,1)))
                mn = min(np.sum(diff,axis=1))
                ind = np.argmin(np.sum(diff,axis=1))
                #diff
                #keyboard
                if (mn<=1e-16):
                    fval = CACHE.xf[ind,n]
                    gval = CACHE.xf[ind,n+1:]
                    if V.m[J] > 0:
                        f_trial = fval + np.sum(np.divide(np.where(gval<0,0,gval),eps))
                    else:
                        f_trial = fval
                    CACHE.hits += 1
                    if outlev >= 2:
                        logger.info('nm_discrete_search: f_trial computed by cache')
                else:
                    f_trial = functpen(y_trial,eps,J,V,CACHE)
                    nf += 1
                    if outlev >= 2:
                        logger.info('nm_discrete_search: f_trial computed by functpen',f_trial)
                        logger.info('                  : fob  = ',CACHE.xf[CACHE.nnf-1,n])
                        gval = CACHE.xf[CACHE.nnf-1,n+1:]
                        logger.info('                  : viol = ',np.sum(np.where(gval<0,0,gval)))
                        logger.info('                  : gval = ',gval)
                        logger.info('                  : eps  = ',eps)
            else:
                f_trial = np.inf


    else:
        alpha = 0
        x = y
        f = np.inf

    return alpha, x, f, nf


##########################################################################
# END OF CODE nm_discrete_linesearch
##########################################################################

def prime_vector(d):
    n = len(d)
    flag = 0
    if(n==1):
        flag = True
        return flag
    temp = np.gcd(np.array(abs(d[0]),dtype=int),np.array(abs(d[1]),dtype=int))
    if(n==2):
        flag = (temp == 1)
        return flag
    for i in np.arange(2,n,1):
        temp = np.gcd(temp,np.array(abs(d[i]),dtype=int))
        #temp = numpy_gcd(temp,abs(d[i]))
        if temp == 1:
            flag = True
            return flag
    if temp != 1:
        flag = False
        return flag


##########################################################################
# END OF CODE prime_vector
##########################################################################


def generate_dirs(ncont,n,D,succ,alpha_tilde,eta,betaLS,Phalton,ihalton):
    #
    # Function generate_dirs
    #
    # Purpose:
    #
    # This function generate new integer directions which are added to set D
    #
    # Inputs:
    #
    # n            : dimension of the problem
    #
    # D            : matrix of current directions (one per each column)
    #
    # alpha_tilde  : array of stepsizes along direction already in D
    #
    # Output:
    #
    # Dout         : [new_direction D]
    #
    # succout      : [0 succ]
    #
    # alpha        : array of stepsizes along the directions in Dout
    #                alpha = [new_step_sizes alpha_tilde]
    #

    # d = rand(n,1)
    # d = d./norm(d)
    #
    # Q = [null(d') d]
    #
    # Dout = [Q D]
    # alpha = [ones(1,n) alpha_tilde]


    mD = np.shape(D)[1]

    for j in range(1000):
        #keyboard
        v = 2*np.asarray(Phalton[ihalton-1], dtype = np.float64) - np.ones(n)
        ihalton += 1
        v = eta*(v/np.linalg.norm(v))

        if (np.linalg.norm(v) < 1e-16):
            break

        #d = abs(round(v)) good if H=norm(d)^2*eye(n,n) - 2*d*d' used
        d = np.round(v)

        #now check whether d is a prime vector
        if prime_vector(d) == True:
            trovato = False
            #check whether d is already in D
            d = np.reshape(d,(len(d),1))
            d = np.concatenate((np.zeros((ncont,1)),d),axis=0)
            DIFF1 = D - np.tile(d,(1,mD))
            DIFF2 = D + np.tile(d,(1,mD))
            if( min ( np.sum(abs(DIFF1),axis=0)) == 0 ) or ( min ( np.sum(abs(DIFF2),axis=0)) == 0 ):
                trovato = True

            if trovato == False:
                H       = d.copy() #norm(d)^2*eye(n,n) - 2*d*d'
                Dout    = np.hstack((H,D))
                succout = np.hstack((np.array(0),succ))
                alpha   = np.hstack((np.array(max(betaLS,max(alpha_tilde))),alpha_tilde))
                iexit   = 1
                return Dout, succout, alpha, iexit, ihalton

    Dout    = D
    succout = succ
    alpha   = alpha_tilde
    iexit   = 0

    return Dout, succout, alpha, iexit, ihalton

##########################################################################
# END OF CODE generate_dirs
##########################################################################
