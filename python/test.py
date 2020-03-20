import gp
import gp_c
import numpy as np


shape=(128,128)
domain=( (0.,5.),(-5,5))
geo=gp.geometry( shape,domain)

#print(geo.positions(0))
y=np.exp(-geo.positions(0)**2 - geo.positions(1)**2)
#print(y.shape)

gp_c.run(y,geo)