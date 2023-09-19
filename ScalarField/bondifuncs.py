from scipy.interpolate import UnivariateSpline
import numpy as np

#trapezoidal rule
def trap_integrator(func,f0,x):
    integration=[f0]
    for i in range(0,len(x)-1):
        dx=x[i+1]-x[i]
    
        integral=dx/2*(func(x[i+1])+func(x[i]))
        
        integration.append(integration[i]+integral)
    return integration

"""def funcc(x):
    return x

def trap_integrator(func,f0,x0,x1):
    dx=(x1-x0)
    integral=dx/2*(func(x0)+func(x1))
    
    return f0+integral

trap_integrator(funcc,0,0,1)"""



#convert central time to bondi time
def converttobondi(time, betascri):

    spl = UnivariateSpline(time, np.multiply(np.exp(np.multiply(betascri,-2)),time), s=0)

    return trap_integrator(spl,0,time)
    
