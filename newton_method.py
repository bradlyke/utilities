"""
A generalized version of the Newton-Raphson method for finding roots. The program
requires some function, f(x), to have a non-zero derivative, and has a number
of options. By default the program will run for 10,000 steps attempting to find
a root whereby |x_i+1 - x_i| <= 1e-8. You can either define your own function
in the calling program, or define your function (and derivative) in
the accompanying newton_function.py. You need to provide a function and its
derivative. In addition, the function and derivative must have the same number
of parameters and the parameters must have the same variable name.

use in another program via
import newton_method as nm

Parameters
----------
func : :class:'function'
	The function that we are trying to find the root for.
fP : :class:'function'
	The derivative of the function that we are trying to find the root for.
x0 : :class:'float'
	The initial guess for the root. This value cannot make either f(x0) or fP(x0)
    equal to zero. If a function has more than one root you need to pick an x0
    that is close to the root you want.
max_steps : :class:'int'
    OPTIONAL - This is the maximum number of iterations the program will try before
    reporting a failure to converge. 10,000 is the default.
thresh : :class:'float'
    OPTIONAL - The value that it checks against for convergence such that
    |f(x_i+1) - f(x_i)| < thresh. It defaults to 1e-8.
check_flag : :class:'boolean'
    OPTIONAL - This determines how it treats a failure to converge. If TRUE, the
    program will return the last value for x it generated as x_final and a failure
    flag, fail_check = 1. If FALSE it returns x_final=np.inf and no flag,
    probably crashing the calling program. It defaults to FALSE.
verbose : :class:'boolean'
    OPTIONAL - If TRUE it will print every value for x_i that it finds at each
    iteration. Useful to seeing how fast a function converges, but do not use
    if you have a large number of max_steps and the function doesn't converge.
    It defaults to FALSE.
**kwargs : :class:'int' or 'float'
    These are the parameters that exist in func and fP such as the coefficients
    of an n-th degree polynomial. These need to be called the same name here
    as they are in func and fP. If no parameters are needed, then leave blank.


Returns
----------
IF CHECK_FLAG=FALSE
:class:'float'
	The root of func() if it converges, or np.inf if it diverges.

IF CHECK_FLAG=TRUE
:class:'float'
    The root of func() and fail_check=0 if it converges, or the last value for
    x_n and fail_check=1 if it diverges. If it diverges you cannot automatically
    trust the value it returns for x_final.
"""

import numpy as np

#This is the actual function for Newton's Method it uses to generate x_i+1
def x_next(func,fP,x_i,**kwargs):
    return x_i - (func(x_i,**kwargs)*fP(x_i,**kwargs)**(-1.0))

#This is the looper part of the function that also returns the final values.
def newt(func,fP,x0,max_steps=10000,thresh=1e-8,check_flag=False,verbose=False,**kwargs):
    counter = max_steps+1 #Number of steps for protection
    fail_check = 0 #Initialize for success.
    check_val = thresh #Threshhold for convergence. Set smaller for better solution
    x_p = x0 #Initialize x_i for the first step
    #And here we loop.
    for i in range(1,counter):
        x_n = x_next(func,fP,x_p,**kwargs) #The meat and potatoes. Get x_i+1
        if verbose==True:
            #If you want to see all of the values of x_i it calculates set
            #verbose to be True. Be careful if you choose a large number of
            #maximum steps.
            print('x_{} = {}'.format(i,x_n))
        if np.abs(x_n - x_p) <= check_val: #Check for final solution.
            #This is what happens if it converged within threshhold.
            x_final = x_n
            break
        else: #If it hasn't converged yet, go to the next iteration.
            x_p = x_n
        if i + 1 == counter:
            #This is the check in case it diverged by the step limit.
            print('Did not converge in {} steps'.format(counter-1))
            if check_flag==True:
                fail_check = 1
                x_final = x_n
            else:
                x_final = np.inf
    #Return some value for x_final (hopefully the convergent one) and the
    #fail_check condition if the check_flag option was used.
    if check_flag==True:
        return x_final,fail_check
    else:
        return x_final

#If you want to see what different arguments and outputs look like with a couple
#test functions, this is what runs. It requires you have newton_function.py in
#the same folder or your python path.
if __name__=='__main__':
    import newton_function as nf
    #This is for recording how you would call this in another function
    #You have to list all of the variables, other than x_t, as keyword arguments
    print('\n------------------\nNon-polynomial Test with Variable Parameters')
    print('------------------')
    x_f = newt(nf.forb,nf.forbp,1.2,M=1.5,e=0.02)
    print('X_f: {}'.format(x_f))

    print('\n------------------\nNon-polynomial Test with failed convergence')
    print('and check_flag=True')
    print('------------------')
    x_ff,fcf = newt(nf.forb,nf.forbp,1.2,max_steps=1,check_flag=True,M=1.5,e=0.02)
    print('X_f: {}\nFail_checkf: {}'.format(x_ff,fcf))

    print('\n------------------\nNon-polynomial Test with failed convergence')
    print('and check_flag=False')
    print('------------------')
    x_ff = newt(nf.forb,nf.forbp,1.2,max_steps=1,M=1.5,e=0.02)
    print('X_ff: {}'.format(x_ff))

    print('\n------------------\nPolynomial Test with Verbose=True')
    print('------------------')
    x_fp = newt(nf.fpoly,nf.fpolyp,1,verbose=True)
    print('Polynomial X_f: {}'.format(x_fp))
