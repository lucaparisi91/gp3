import numpy as np
import yt
from math import *
import glob
import re
from contextlib import redirect_stdout,redirect_stderr
import io
import os
import sys

class geometry:
	def __init__(self,shape,domain,symmetry="none",coordinates="cartesian",bc=None,grown_shape=None):
		self.shape=shape
		self.lower_edges= [bound[0] for bound in domain]
		self.higher_edges= [bound[1] for bound in domain]
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



