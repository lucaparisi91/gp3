import gp_c
import numpy as np
import yt
from math import *
import glob
import re
from contextlib import redirect_stdout,redirect_stderr
import io
import os
import sys

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



class gp_simulation_cylindrical:
	def __init__(self,folder_name_real,folder_name_imag,geo):
		out_stream = io.StringIO()
		
		self.yt_input_real=yt.load(folder_name_real);
		self.yt_input_imag=yt.load(folder_name_imag);
		

		data_real=self.yt_input_real.all_data()
		data_imag=self.yt_input_imag.all_data()
		self.geo=geo
		self.r=np.array(data_real.fcoords)[:,0]
		self.z=np.array(data_real.fcoords)[:,1]
		self.phi_real=self.yt_input_real.all_data()["phi"]
		self.phi_imag=self.yt_input_imag.all_data()["phi"]


		self.r=self.r.reshape(geo.shape)
		self.z=self.z.reshape(geo.shape)
		self.phi_real=self.phi_real.reshape(geo.shape)
		self.phi_imag=self.phi_imag.reshape(geo.shape)

		self.time=float(self.yt_input_imag.current_time)

	def _integrate(self,y):
		return 2*pi*np.sum(y*self.r)*self.geo.step[0] *self.geo.step[1]
	def average(self,f):
		y=f(self.r,self.z)*(self.phi_real**2 + self.phi_imag**2)
		return self._integrate(y)
	def norm(self):
		return self.average(lambda r,z: 1.)
	def density(self):
		return np.array(self.phi_real**2 + self.phi_imag**2)




class gp_simulations:

	def __init__(self,folder,geo):
		self.folder=folder
		self.geo=geo
		self.phi_reals=glob.glob(self.folder + '/phi_real*')
		
		indices=[]
		pattern=r'.*phi_real([0-9]+)'
		for filename in self.phi_reals:
			result = re.match(pattern, filename)
			if result:
				indices.append(int(result.group(1)))
		self.indices=np.array(indices)

		self.phi_imags= [ re.sub("phi_real","phi_imag",phi_real)  for phi_real in self.phi_reals]

		sorted_data=np.array(sorted(zip(self.indices,self.phi_reals,self.phi_imags) ))

		self.indices=sorted_data[:,0]
		self.phi_reals=sorted_data[:,1]
		self.phi_imags=sorted_data[:,2]
	def __len__(self):
		return len(self.phi_reals)

	def __iter__(self):
		
		for phi_real,phi_imag in zip(self.phi_reals,self.phi_imags):
			yield gp_simulation_cylindrical(phi_real,phi_imag,self.geo)
	def __repr__(self):

		return "<folder=" + str(self.folder) +", len=" + str(len(self))+ "= >"
	def __getitem__(self,index):
			return gp_simulation_cylindrical(self.phi_reals[index],self.phi_imags[index],self.geo)