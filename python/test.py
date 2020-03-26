import gp
import gp_c
import numpy as np
from importlib import reload
import system
from colorama import Fore, Back, Style,init
from math import *

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
sim=gp.gp_simulation_cylindrical("out/initial_real","out/initial_imag",geo)
# excites the dipole mode
delta_psi=np.exp(- 1j * pi  * np.tanh(sim.z/0.8))

psi0=(sim.phi_real + 1j * sim.phi_imag)*delta_psi

sim.phi_real=np.real(psi0)
sim.phi_imag=np.imag(psi0)

print( Fore.YELLOW + "Running..." + Fore.RESET)
gp_c.run(sim.phi_real,sim.phi_imag,geo)