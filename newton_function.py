#For this script, f and fp have to be the function and the derivative
#of the function. Also, all parameters must be the same letter. Both functions
#must have the same number of parameters.

import numpy as np

def f(x_t,a):
    return a*x_t**(3.0) - 27

def fp(x_t,a):
    return 3*a*x_t**(2.0)


################################################################################
#These two following functions are for testing newton_method.py if you want
#to see how the function works different types of functions with different types
#of inputs
################################################################################
def forb(x_t,M,e):
    return x_t - (e*np.sin(x_t)) - M

def forbp(x_t,M,e):
    return 1 - (e*np.cos(x_t))

#This is in case you want to test with a simple polynomial that will converge
#in six steps.
def fpoly(x_t):
    return x_t**(2.0) - 9

def fpolyp(x_t):
    return 2*x_t
