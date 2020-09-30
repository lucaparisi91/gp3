import numpy as np

def laplacianSpherical(y, domain=[0,1],order=2):
    y2=y*0
    dx=(domain[1] - domain[0] )/len(y)
    for i in range(1,len(y)-1):
        y2[i] = ( y[i + 1] + y[i-1] -2*y[i] )/(dx**2)
    return y
    