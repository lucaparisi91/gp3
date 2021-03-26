import gp
import unittest 
import gpIO_c
import numpy as np
import matplotlib.pylab as plt 
from math import *


class testMultifab(unittest.TestCase):

    def test_saveMultifab(self):
        geo=gp.geometry((40,40,40),[(-2,2),(-2,2),(-2,2)])
        r2=geo.positions(0)**2 + geo.positions(1)**2 + geo.positions(2)**2
        y=(np.exp(-r2) + 0*1j).reshape(geo.shape[0],geo.shape[1],geo.shape[2],1)
        wave=gp.wavefunction(geo,[y] ,name="phi" )
        wave.save("phiTest")
        wave2=gp.loadWavefunction("phiTest")
        y2=wave2.phis[0]
        
        assert(np.max(y2 - y)==0 )


        


if __name__ == "__main__":
    unittest.main()    
