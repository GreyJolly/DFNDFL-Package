############################################################################################
#
#description  (problems for mixed-integer problems with AT LEAST
#              2 discrete variables and 2 continuous variables)
#
# a problem object is a structured type that has the following attributes:
# name   : string - name of the problem
# startp : numpy array - the starting point for the continuous problem
# lb     : numpy array - the lower bounds of the continuous problem
# ub     : numpy array - the upper bounds of the continuous problem
# nint   : int - the number of discrete variables (>= 2)
# ncont  : int - the number of continuous variables (>= 2), N.B. n = nint+ncont
#        : BEWARE the variables are so intended
#        :        x[0] ... x[ncont-1] are the continuous variables
#        :    x[ncont] ... x[n-1]     are the discrete variables
# lbmix  : numpy array - the actual lower bounds of the mixed integer problem
# ubmix  : numpy array - the actual upper bounds of the mixed integer problem
# x_initial : numpy array - the actual initial point of the mixed integer problem
# feval  : function handle - function to compute the objective function value
#        : N.B. the point must be reconstructed through the use of reconstruct_xmix
#        :      before calling feval!
# m      : dictionary with entries char => number of constraints
#        : it is equal to {'a': n-2, 'b': n-2, 'c': n-1, 'd': n-1, 'e': n-2, 'f': 1, 'z': 0}
#        : it is used to record the number of constraints for the given problem and for
#        : each of the six families of constraints a,b,c,d,e,f. 'z' means no constraints
#
############################################################################################

import numpy as np

class Cache:
    # nnf = number of entries used in numpy array xf
    # xf  = numpy array where computed values are stored
    # hits= number of times item found in cache
    # maxf= maximum number of function values that can be stored in F
    # F   = numpy array of function values used for performance profiles
    nnf  = 0
    hits = 0
    xf   = []
    maxf = 0
    F    = []
    def __init__(self,max_fun,dim):
        self.nnf  = 0
        self.hits = 0
        self.F    = []
        self.xf   = np.inf*np.ones((10*max_fun,dim)) # dim = n+ncont+m+1
        self.maxf = max_fun

class Problem:
    name = ""
    startp = None
    lb = None
    ub = None
    n = 0
    nint = 0
    ncont = 0
    lbmix = None
    ubmix = None
    x_initial = None
    xmix = None
    feval = None
    m = {}

    def reconstruct_xmix(self,x):
        self.xmix[self.ncont:] = self.lb[self.ncont:] + ((self.ub[self.ncont:] - self.lb[self.ncont:])/(self.ubmix[self.ncont:] - self.lbmix[self.ncont:]))*(x[self.ncont:]-self.lbmix[self.ncont:])
        self.xmix[:self.ncont] = x[:self.ncont]

    def func_f(self, x):
        self.reconstruct_xmix(x)
        return self.feval(self.xmix)

    def fconstr_a(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return  (3-2*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] - 2*self.xmix[J+2] + 1

    def fconstr_b(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return  (3-2*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] - 2*self.xmix[J+2] + 2.5

    def fconstr_c(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-1)
        return self.xmix[J]**2 + self.xmix[J+1]**2 + self.xmix[J]*self.xmix[J+1] - 2*self.xmix[J] - 2*self.xmix[J+1] + 1

    def fconstr_d(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-1)
        return self.xmix[J]**2 + self.xmix[J+1]**2 + self.xmix[J]*self.xmix[J+1] - 1

    def fconstr_e(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return (3-0.5*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] -2*self.xmix[J+2] +1

    def fconstr_f(self,x):
        self.reconstruct_xmix(x)
        J = np.arange(len(self.xmix)-2)
        return np.array([np.sum((3-0.5*self.xmix[J+1])*self.xmix[J+1] - self.xmix[J] -2*self.xmix[J+2] +1)])

    def fconstr_z(self,x):
        return np.array([-1.0])

    def __init__(self,name,startp,lb,ub,nint,ncont,lbmix,ubmix,x_initial,feval):
        self.name = name
        self.startp = startp
        self.lb = lb
        self.ub = ub
        self.nint = nint
        self.ncont = ncont
        self.n = nint + ncont
        self.lbmix = lbmix
        self.ubmix = ubmix
        self.x_initial = x_initial
        self.xmix = np.zeros(self.n)
        self.feval = feval
        if self.n >= 3:
            self.m = {
                'a': self.n - 2,
                'b': self.n - 2,
                'c': self.n - 1,
                'd': self.n - 1,
                'e': self.n - 2,
                'f': 1,
                'z': 0
            }
        elif self.n <= 2:
            self.m = {
                'c': self.n - 1,
                'd': self.n - 1,
                'f': 1,
                'z': 0
            }
        else:
            raise ValueError("Unexpected value of n: {}".format(self.n))