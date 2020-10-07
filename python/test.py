import gp
import gp_c
import numpy as np
from importlib import reload
from colorama import Fore, Back, Style,init
from math import *


init()

shape=(128,128,128)
domain=( (-5,5.),(-5,5), ( -5 , 5 ))
geo=gp.geometry( shape,domain)
#print(geo.positions(0))
y=np.exp(- (geo.positions(0)**2 + geo.positions(1)**2 + geo.positions(2) ** 2 ))

# setup of the simulation
settings = {
    "geometry" : {"shape" : shape , "domain" : domain},
     "algorithm" : {"stepper" : "RK4","stepsPerBlock":100,"blocks" : 1,"timeStep" : 1e-3, "imaginary" : True }
     }

gp_c.run(y,y*0, geo)
