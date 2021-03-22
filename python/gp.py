import numpy as np
from math import *
import h5py


class geometry:
	def __init__(self,shape,domain,symmetry="none",coordinates="cartesian",bc=None,grown_shape=None):
		'''
		Creates the geometry of a cell centered grid. Does not store fields.
		'''

		self.shape=shape
		self.lower_edges= np.array([float(bound[0]) for bound in domain])
		self.higher_edges= np.array([float(bound[1]) for bound in domain])

		self.step = [ (h - l)/n for l,h,n in zip(self.lower_edges,self.higher_edges,self.shape) ]
		self.dimensions=len(shape)
		self.symmetry=symmetry
		self.coordinates=coordinates
		self.bc=bc
		if grown_shape is None:
			grown_shape=shape
		self.grown_shape=grown_shape

	def positions(self,axis):

		x=np.arange(self.lower_edges[axis],self.higher_edges[axis],self.step[axis])
		x+=0.5*self.step[axis]

		if self.dimensions == 1 :
			return x

		if self.dimensions == 2:
			if axis==0:
				return np.outer(x,np.ones(shape=(self.shape[1])))
			else :
				if axis==1:
					return np.outer(np.ones(shape=self.shape[0]),x)
		
		if self.dimensions == 3 :
			if axis == 0:
				return np.outer(x,np.ones(shape=(self.shape[1],self.shape[2]))).reshape(self.shape)
			else:
				if axis == 1:
					tmp=np.outer(x,np.ones(shape=(self.shape[2])))
					return np.outer(np.ones(shape=(self.shape[0])),tmp).reshape(self.shape)
				else:
					if axis == 2:
						return  np.outer(np.ones(shape=(self.shape[0],self.shape[1])),x).reshape(self.shape)



class wavefunction:

	def __init__(self,geo=None,phi=None):
		'''
		Stores an array of grids(geometries) and associated wavefunction 
		value arrays.
		'''

		self.geometry=geo
		self.value=phi

		
	def saveHDF5(self,filename):
		f= h5py.File(filename,'w')
		geo=self.geometry
		phi=self.value
		g=f
		g.create_dataset( name="phi",data=phi,dtype=np.float64)
		#g.attrs["lower_edges"]=geo.lower_edges
		g.attrs.create(name="lower_edges", data=geo.lower_edges, shape=geo.lower_edges.shape, dtype=np.float64)
		g.attrs.create(name="higher_edges", data=geo.higher_edges, shape=geo.higher_edges.shape, dtype=np.float64)
		g.attrs.create("shape", data=geo.shape, shape=(geo.dimensions,), dtype=np.int32)

		f.close()
	
	def loadHDF5(self,filename):
		f= h5py.File("gaussian.hdf5",'r')
		g=f
		domain =[
				[left, right] for left,right in zip(g.attrs["lower_edges"],g.attrs["higher_edges"]) 
				]
		shape=g.attrs["shape"]
		self.phi=g["phi"][()]
		self.geometry=gp.geometry(shape=shape,domain=domain)

		
		f.close()



			
			









