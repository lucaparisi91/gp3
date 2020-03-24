import gp
import gp_c
import numpy as np
from importlib import reload
import system
from colorama import Fore, Back, Style,init

init()

potential=r'potential(psi_real,psi_imag,x){ 100. * (psi_real*psi_real + psi_imag*psi_imag) + 0.5 *( x[0]*x[0] + x[1]*x[1])}'

shape=(128,128)
domain=( (0.,5.),(-5,5))
geo=gp.geometry( shape,domain)
#print(geo.positions(0))
y=np.exp(-geo.positions(0)**2 - geo.positions(1)**2)
#print(y.shape)

print( Fore.YELLOW + "Compiling..." + Fore.RESET)
system.compileLocalPotential(potential,"cylindrical")

reload(gp_c)

print( Fore.YELLOW + "Running..." + Fore.RESET)
gp_c.run(y,geo)

