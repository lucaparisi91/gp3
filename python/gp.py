import gp_c
import numpy as np

class geometry(gp_c.geometry):
	def __init__(self,shape,domain):

		super().__init__(shape)
		if len(shape) != 2:
			raise NotImplementedError()

		self.lower_edges= [bound[0] for bound in domain]
		self.higher_edges= [bound[1] for bound in domain]
		self.step = [ (h - l)/n for l,h,n in zip(self.lower_edges,self.higher_edges,self.shape) ]
		self.dimensions=len(shape)
	def positions(self,axis):

		x=np.arange(self.lower_edges[axis],self.higher_edges[axis],self.step[axis])
		x+=0.5*self.step[axis]
		if axis==0:
			return np.outer(x,np.ones(shape=(self.shape[1])))
		else :
			if axis==1:
				return np.outer(np.ones(shape=self.shape[0]),x)

