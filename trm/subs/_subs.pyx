#cython: embedsignature=True
# distutils: language = c++
cimport numpy as np
import numpy as np

cdef extern from "trm/subs.h" namespace "Subs":
    double c_voigt "Subs::voigt" (double a, double v, double eps)
    double c_gammq "Subs::gammq" (double a, double x)
    
def voigt(a, v, eps=1.e-8):
    """returns Voigt function

    v can be a float or a 1D array"""
    v = np.atleast_1d(v)
    cdef unsigned n = v.shape[0]
    cdef unsigned int i
    v = np.ascontiguousarray(v)
    cdef np.ndarray out = np.empty(n, dtype=v.dtype)
    
    for i in range(n):
        out[i] = c_voigt(a,v[i],eps)
    if n==1:
        return out[0]  
    else:
        return out
    
def gammq(a, x):
    """returns incomplete gamma function
    
    x can be a float or a 1D array"""
    x = np.atleast_1d(x)
    cdef unsigned n = x.shape[0]
    cdef unsigned int i
    x = np.ascontiguousarray(x)
    cdef np.ndarray[double, ndim=1] out = np.empty(x, dtype=np.double)
    
    for i in range(n):
        out[i] = c_gammq(a,x[i])
        
    if n==1:
        return out[0]  
    else:
        return out