import gp
import gp_c
import numpy as np
from importlib import reload
import system

potential=r'potential(psi_real,psi_imag,x){ 100. * (psi_real*psi_real + psi_imag*psi_imag) + 0.5 *( x[0]*x[0] + x[1]*x[1])}'

shape=(256,256)
domain=( (0.,5.),(-5,5))
geo=gp.geometry( shape,domain)
#print(geo.positions(0))
y=np.exp(-geo.positions(0)**2 - geo.positions(1)**2)
#print(y.shape)
system.compileLocalPotential(potential,"cylindrical")

reload(gp_c)
gp_c.run(y,geo)

