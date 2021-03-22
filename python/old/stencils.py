import numpy as np
import scipy as sp 
from math import *

def stencil( stencilPoints , d ):
    stencilPoints=np.array(stencilPoints)
    N=len(stencilPoints)
    rows=[ stencilPoints**i for i in range(0,N)  ]
    A=np.array(rows).reshape(N,N)
    b=np.zeros(len(rows) )
    b[d]=factorial(d)

    a=np.linalg.solve(A,b)
    return a

if __name__ == "__main__":
    print(stencil([1,3,5] , 1))