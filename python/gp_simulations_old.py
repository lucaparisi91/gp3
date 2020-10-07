class gp_simulation:
	def __init__(self,folder_name_real,folder_name_imag,geo):
		out_stream = io.StringIO()
		
		self.yt_input_real=yt.load(folder_name_real)
		self.yt_input_imag=yt.load(folder_name_imag)

		data_real=self.yt_input_real.all_data()
		data_imag=self.yt_input_imag.all_data()
		self.geo=geo

		if geo.coordinates == "cylindrical" :
			self.r=np.array(data_real.fcoords)[:,0]
			self.z=np.array(data_real.fcoords)[:,1]
			self.r=self.r.reshape(geo.shape)
			self.z=self.z.reshape(geo.shape)


		else:

			if geo.coordinates == "spherical":
				self.r=np.array(data_real.fcoords)[:,0]

			else:

				if (  geo.dimensions >= 1 ):
					self.x=np.array(data_real.fcoords)[:,0].reshape(geo.shape)

					if (geo.dimensions >= 2):
						self.y=np.array(data_real.fcoords)[:,1].reshape(geo.shape)

						if (geo.dimensions >= 3):
							self.z=np.array(data_real.fcoords)[:,2].reshape(geo.shape)
					
		
		self.phi_real=self.yt_input_real.all_data()["phi"]
		self.phi_imag=self.yt_input_imag.all_data()["phi"]

		self.phi_real=self.phi_real.reshape(geo.shape)
		self.phi_imag=self.phi_imag.reshape(geo.shape)

		self.time=float(self.yt_input_imag.current_time)

	def _coords(self):

		if self.geo.dimensions == 3:
			return (self.x, self.y,self.z)
		if self.geo.coordinates == "cylindrical":
			return (self.r,self.z)
		if self.geo.coordinates == "spherical":
			return (self.r)
		

	def _integrate(self,y):
		if self.geo.coordinates == "cylindrical":
			return 2*pi*np.sum(y*self.r)*self.geo.step[0] *self.geo.step[1]
		else:

			if self.geo.coordinates == "spherical":
				return 4*pi*np.sum(y*self.r**2)*self.geo.step[0]
			else:
				res = np.sum(y)
				dv=1
				for i in range(self.geo.dimensions):
					dv*=self.geo.step[i]
				return res*dv

	
	def average(self,f):
		if hasattr(f, '__call__'):

			coords = self._coords()
			integrand=f(*coords)
		else:
			integrand=f


		y=integrand*(self.phi_real**2 + self.phi_imag**2)
		return self._integrate(y)
	def norm(self):
		return self.average(1.)
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