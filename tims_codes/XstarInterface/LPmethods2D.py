import os
import sys
import timeit
import numpy as np
import pickle
from scipy.optimize import brentq,newton,bisect

Fund = {'mp':1.6726231e-24, 'me':9.1093898e-28, 'c':3.0e10, 'h':6.6260755e-27, \
				'e':4.8032068e-10, 'kb':1.380658e-16, 'G':6.67259e-8}
# calculate Thomson cross section 
Fund['sig_th'] = (8.*np.pi/3.)*(Fund['e']**2/(Fund['me']*Fund['c']**2))**2

def get_default_paths(dump,run_name,ion):
	# create path names to pickled objects
	dump = dump
	opacity_folder = 'ion_opacities/' 
	out_folder = 'output/' 
	output_folder = out_folder + run_name + '/' + ion
	pickled_sigmas = output_folder + '/sigmas/'
	pickled_profiles = output_folder + '/profiles/' 	
	paths = [pickled_sigmas,pickled_profiles]
	for path in paths:
		if not os.path.exists(path):
			os.makedirs(path)
			
	default_paths = {}
	default_paths['sigmas'] =  pickled_sigmas + 'sigmas_' + str(dump) + '.p'
	default_paths['profiles'] = pickled_profiles
	# paths to stored calculations
	default_paths['tau_vs_dist'] = pickled_profiles + 'tau_xorys_' + str(dump)  + '.p'
	default_paths['taus'] = pickled_profiles + 'taus_' + str(dump)  + '.p'	
	default_paths['fluxes'] = pickled_profiles + 'fluxes_' + str(dump) + '.p'
	default_paths['cog'] = pickled_profiles + 'cog_' + str(dump) + '.p'
	default_paths['cog_normed'] = pickled_profiles + 'cog_normed_' + str(dump) + '.p'
	
	return default_paths
	
def bilinear_interp(q11, q21, q12, q22, 
					x1, x2, y1, y2, x, y):
	x2x1 = x2 - x1
	y2y1 = y2 - y1
	x2x  = x2 - x
	y2y  = y2 - y
	yy1  = y - y1
	xx1  = x - x1
	return (q11*x2x*y2y + q21*xx1*y2y + q12*x2x*yy1 + q22*xx1*yy1)/(x2x1 * y2y1)

class opacity_table:
	def __init__(self,data_file):
		self.data = []
		"""
		Read in the ion data.
		"""
		print("[opacity_table]:\nReading the data file {}.".format(data_file))
		with open(data_file) as f:
			for line in f:
				self.data.append(line.split())
		print("[opacity_table]:\nDone!")
		"""
		break it down into smaller chunks
		"""
		self.data   = np.array(self.data[1:])

		self.xis    = self.data[:,1].astype(float)
		self.temps  = np.log10(self.data[:,0].astype(float)*1e4)

		self.data      = np.array(self.data[:,2:])
		self.num_lines = len(self.data[0])/3
		self.num_temps = len(np.unique(self.temps))
		self.num_xis   = len(np.unique(self.xis))

		line_arr = []

		for i in range(int(self.num_lines)):
			line_arr.append(self.data[0][1+i*3])
		self.lines = line_arr
		self.opact_dict = dict()

		print("[opacity_table]:\nCreating interpolation functions for the following lines:")
		print(self.lines)
		
		for i,line in enumerate(line_arr):
			self.opact_dict.update({line :
				self.data[:,2+i*3].astype(float)})
	"""
	Function to get the opacities after we've 
	initialized the tables and such.

	we get opacities by passing log(xi) and log(T)
	"""
	def get_opacity(self,xi,T,line=None,verbose=False):
		if not line:
			line = self.lines[0]
# 		print("[opacity_table]:\n(xi,T) = ({},{})".format(xi,T))
		l1 = np.where(np.unique(self.temps) < T)[0][-1]*self.num_xis
		l2 = np.where(np.unique(self.xis) < xi)[0][-1]
		l3 = l1 + self.num_xis
		Nmax = self.opact_dict[line].shape[0]
		
# 		l3 = min(l3,Nmax-1)
# 		print("Nmax = {}".format(Nmax))
# 		print("num_xis = {}, (l1,l2,l3) = ({},{},{})".format(self.num_xis,l1,l2,l3))
# 		q11 = self.opact_dict[line][l1+l2]
# 		q12 = self.opact_dict[line][l1+l2+1]
# 		q21 = self.opact_dict[line][l1+l2]
# 		q22 = self.opact_dict[line][l3+l2+1]
# 		x1  = self.temps[l1+l2]
# 		x2  = self.temps[l3+l2]
# 		y1  = self.xis[l1+l2]
# 		y2  = self.xis[l1+l2+1]
		
		# Shalini's version (March 2020)
		q11 = self.opact_dict[line][l1+l2]
		q12 = self.opact_dict[line][l1+l2+1]
		q21 = self.opact_dict[line][l3+l2]
		q22 = self.opact_dict[line][l3+l2+1]
		x1 = self.temps[l1+l2]
		x2 = self.temps[l3+l2]
		y1 = self.xis[l1+l2]
		y2 = self.xis[l3+l2+1]
		
		if verbose:
			print(q11,q12,q21,q22,x1,x2,y1,y2)

		ans = bilinear_interp(q11, q21, q12, q22, 
							  x1, x2, y1, y2, T, xi)
		if verbose:
			print(ans)
		return ans
		
					
class LineProfiles:
	def __init__(self,dump,run_name,simdata,units,line_data,Nwidths=3.,Nv=1000,v_wind_cgs=0.,d_fac=1.,tau_cutoff=0.5,y_observer=0,soln_is_1D=0,recalc=0):
		"""
			Inputs:
			dump		: integer marking the simulation output dump
			simdata		: dictionary containing all simulation variables
			units  	: dictionary containing cgs unit conversions
			line_data	: dictionary containing atomic data
			Nwidths		: number of thermal widths around line center
			Nv			: number of points to use for velocity axis
			v_wind		: bulk velocity offset (cloud may live in a wind)
	
			Note: it is assumed that lines wavelengths are in Angstroms
		""" 
		self.soln_is_1D = soln_is_1D
		self.Fund = Fund
		
		print("[LineProfiles constructor]:\nCreating LineProfiles instance.")
		
		# get ion, desired lines, and opacity data 
		self.ion = line_data['ion']
		self.A_Z = line_data['A_Z']
		self.lines = line_data['lines']
		self.table_path = line_data['path'] + line_data['SED'] + '/'
		pickled_ion = self.table_path + self.ion + '.p'
		try:
			print("[LineProfiles constructor]:\nAttempting to open {}".format(pickled_ion))
			self.XSTAR_LUT = pickle.load(open(pickled_ion, "rb"))
			print("[LineProfiles constructor]:\nSuccess!")
		except:
			print("[LineProfiles constructor]:\nBuilding opacity LUT using table\n{}".format(self.table_path + line_data['table']))
			self.XSTAR_LUT = opacity_table(self.table_path + line_data['table'])
			pickle.dump(self.XSTAR_LUT,open(pickled_ion,"wb"),protocol=2)
			
		# create path names to pickled objects
		default_paths = get_default_paths(dump,run_name,self.ion)
		self.pickled_sigmas =  default_paths['sigmas']
		self.pickled_profiles = default_paths['profiles']
		# paths to stored calculations (optional)
		self.pickled_tau_xorys = default_paths['tau_vs_dist']
		self.pickled_taus = default_paths['taus']
		self.pickled_fluxes = default_paths['fluxes']
		self.pickled_cog = default_paths['cog']
		self.pickled_cog_normed = default_paths['cog_normed']
			
		# initialize new line_data dictionary
		self.line_data = {}
		for line in self.lines:
			self.line_data[line] = {}
		
		# record rest frequencies and oscillator strengths
		for line in self.lines:
			lambda0 = float(line)*1e-8  # converts Angstroms to cgs
			nu0 = Fund['c']/lambda0
			self.line_data[line]['nu0'] = nu0
# 			self.line_data[line]['f'] = line_data['f_b']
# 			self.line_data[line]['f'] = line_data['f_r']
		
		# construct the freq shift (velocity) axis
		vc = simdata['vc']
		self.v0_cgs = v_wind_cgs + vc*units['v'] # average bulk velocity
		v0 = self.v0_cgs/Fund['c']
		y_max = Nwidths*units['v']/Fund['c']
		if Nv > 1:
			self.y = np.linspace(v0-y_max,v0+y_max,Nv) # y = (nu-nu_0)/nu_0 = v/c
		else:
			self.y = np.array([v0])
			
		if soln_is_1D:
			y_min = -5e7/Fund['c']
			y_max = 5e6/Fund['c']
			self.y = np.linspace(y_min,y_max,Nv)
			
		self.vel = self.y*Fund['c']/1e5 # vel in units of km/s
		
		# calculate the flow velocity
		if y_observer:
			v_los = simdata['velocity_y'] # LOS velocity is the x-component only
		else:
			v_los = simdata['velocity_x']
		vlocal_cgs = v_los*units['v']
		self.v = -(v0 + vlocal_cgs/Fund['c']) # minus sign for blue shift
		
		# calculate thermal line width
# 		T_cgs = units['T']*self.T
		T_cgs = simdata['temperature']
		m_i = line_data['m_i']*Fund['mp']
		vt2 = 2*Fund['kb']*T_cgs/m_i
		self.vt_cgs = np.sqrt(vt2)
		self.vt = self.vt_cgs/Fund['c']
				
		# Lengths of arrays
		res = self.vt.shape
		if len(res) == 1:
			self.Nx = res[0]
			self.Ny = 1
		else: #2D data
			self.Nx = res[1]
			self.Ny = res[0]
		self.Nv = Nv 		  # our frequency resolution
		
		# log the (xi,T) values
		xi_cgs = simdata['xi']
		self.log_xi = np.log10(xi_cgs)
		self.log_T  = np.log10(T_cgs)
		
		# calculate dx from grid spacing
		if self.Ny == 1:
			x = simdata['X']
			self.dx = units['L']*(x[1]-x[0])
		else: #2D calculation
			x = simdata['X']
			y = simdata['Y']
			self.dx = units['L']*(x[0,1]-x[0,0])
			self.dy = units['L']*(y[1,0]-y[0,0])
			self.xarr,self.yarr = x,y
		
		# record the number density
		self.n_cgs = d_fac*units['n']*simdata['density']
		self.L0 = units['L']
		
		# calculate the column density
		dx = x[1:] - x[:-1]
		self.N = np.cumsum(0.5*(self.n_cgs[1:] + self.n_cgs[:-1])*dx)
		
		self.norm = 1./np.sqrt(np.pi)
		self.units = units
		
		# observer is in x-direction by default
		self.observer_at_90deg = False
		if y_observer:
			self.observer_at_90deg = True
			
		# by default, all rows will be summed over
		self.i_min = 0
		self.i_max = self.Nx - 1
		self.j_min = 0
		self.j_max = self.Ny - 1
		
		# cutoff between optically thick/thin
		self.DOWN_AVERAGING = False
		self.tau_cutoff = tau_cutoff
		
		# density assumed by xstar
		self.n_xstar = 1e8 #1e10
		
		# load line-center opacity LUT or calculate from xstar table
		self.recalc = recalc
		if recalc:
			self.line_sigmas = self.linecenter_opacities(self.pickled_sigmas)
		else:
			try:
				self.line_sigmas = pickle.load(open(self.pickled_sigmas,"rb"))
			except:
				print("[LineProfiles constructor]:\nCalculating opacities using XSTAR data...")
				self.line_sigmas = self.linecenter_opacities(self.pickled_sigmas)
				
			
	def linecenter_opacities(self,pickled_sigmas):
		""" returns 2D line-center opacity array (frequency independent) """
		line_sigmas = []
		for line in self.lines:
			sigma = np.zeros_like(self.vt)
			line_sigmas.append(sigma)
		func = self.XSTAR_LUT.get_opacity
		print("[linecenter_opacities]:\nlog(xi).shape = {}, log(T).shape = {}".format(self.log_xi.shape,self.log_T.shape))
		
		start = timeit.default_timer()
		for i in range(self.Nx):
			sys.stdout.write("\rLooking up line-center opacity for index %d " % i)
			sys.stdout.flush()
			if self.Ny > 1:
				for j in range(self.Ny):
					for sigma,line in zip(line_sigmas,self.lines):
# 						print line
						sigma[j][i] = func(self.log_xi[j][i],self.log_T[j][i],line=line)
						sigma[j][i] *= (self.n_cgs[j][i]/self.n_xstar)
			else:
				for sigma,line in zip(line_sigmas,self.lines):
# 					print line
					sigma[i] = func(self.log_xi[i],self.log_T[i],line=line)
					sigma[i] *= (self.n_cgs[i]/self.n_xstar)
			
		stop = timeit.default_timer()
		runtime = round((stop - start)/60.,2)
		print("[linecenter_opacities]:\nOpacity table construction took {} minutes.".format(runtime))
		
		print("[linecenter_opacities]:\nPickling line-center opacities to {}".format(self.pickled_sigmas))
		pickle.dump(line_sigmas, open(pickled_sigmas, "wb"), protocol=2)
		
		return line_sigmas
	
	def phi_nu(self,y):
		""" returns 2D line shape array at a single frequency """
		arg = (y - self.v)/self.vt
		return self.norm*np.exp(-arg**2)
	
	def alpha_nu(self,phi_nu,alpha0):
		""" returns 2D opacity array at a single frequency """
# 		alpha = phi_nu*alpha0/(1.+self.v)
		alpha = phi_nu*alpha0  # 3/21/2020: could it just be this??
		return alpha
	
	# def tau_nu_1(self,alphaofx):
	# 	""" returns 1D x-integrated optical depth array at a single frequency
	# 		note: this calc assumes a uniform grid (i.e. dx is pulled 
	# 		outside the sum)
	# 	"""
		
	# 	# first calculate tau(x,y)
	# 	tau_nu_xy = np.array_like(alphaofx)
	# 	if self.Ny == 1:
	# 		tau_nu_xy[0] = 0.
	# 		for n in range(1,self.Nx):
	# 			tau_nu_xy[n] = np.sum(0.5*(alphaofx[0:n] + alphaofx[1:n+1]))*self.dx
	# 	else: 
	# 		tau_nu_xy[:,0] = 0.
	# 		for n in range(1,self.Nx): #fix this index
	# 			tau_nu_xy[:,i] = np.sum(0.5*(alphaofx[:,0:n] + alphaofx[:,1:n+1]),axis=1)*self.dx
		
	# 	# now integrate flux over y to arrive at tau_eff(x)
	# 	if self.Ny == 1: # nothing to integrate
	# 		tau_nu_x = tau_nu_xy
	# 	else: 
	# 		flux_xy = np.exp(-tau_nu_xy)
	# 		flux_x = np.sum(0.5*(flux_xy[self.j_min:self.j_max,:] + flux_xy[self.j_min+1:self.j_max+1,:]),axis=0)
	# 		flux_x /= (1+self.j_max-self.j_min)
	# 		tau_nu_x = -np.log(flux_x)
			
	# 	return tau_nu_x
	
	def tau_nu_x(self,alphas):
		""" calculates tau_nu(x) at a single frequency
			note: this calc assumes a uniform grid (i.e. dx is pulled 
			outside the sum)
		"""
		
		if self.Ny == 1:
			tau_nu_x = np.sum(0.5*(alphas[0:-1] + alphas[1:]))*self.dy
		else: 
			tau_nu_x = np.sum(0.5*(alphas[0:-1,:] + alphas[1:,:]),axis=0)*self.dy
		
		return tau_nu_x
	

	def tau_nu_y(self,alphas):
		""" calculates tau_nu(y) at a single frequency
			note: this calc assumes a uniform grid (i.e. dx is pulled 
			outside the sum)
		"""
		
		if self.Ny == 1:
			tau_nu_y = np.sum(0.5*(alphas[0:-1] + alphas[1:]))*self.dx
		else: 
			tau_nu_y = np.sum(0.5*(alphas[:,0:-1] + alphas[:,1:]),axis=1)*self.dx
		
		return tau_nu_y
		
	def cumtau_nu(self,alphas):
		""" calculates cumulative tau_nu(y) at a single frequency
			note: this calc assumes a uniform grid (i.e. dx is pulled 
			outside the sum)
		"""
		
# 		return np.cumsum(0.5*(alphas[0:-1] + alphas[1:]))*self.dx
		return np.cumsum(alphas)*self.dx

	def N_H_x(self,n_H):
		""" calculates column density N_H
			note: this calc assumes a uniform grid (i.e. dx is pulled 
			outside the sum)
		"""
		
		if self.Ny == 1:
			N_H_x = np.sum(0.5*(n_H[0:-1] + n_H[1:]))*self.dy
		else: 
			N_H_x = np.sum(0.5*(n_H[0:-1,:] + n_H[1:,:]),axis=0)*self.dy
		
		return N_H_x

	def N_H_y(self,n_H):
		""" calculates column density N_H
			note: this calc assumes a uniform grid (i.e. dx is pulled 
			outside the sum)
		"""
		
		if self.Ny == 1:
			N_H_y = np.sum(0.5*(n_H[0:-1] + n_H[1:]))*self.dx
		else:
			N_H_y = np.sum(0.5*(n_H[:,0:-1] + n_H[:,1:]),axis=1)*self.dx
		
		return N_H_y	
			
	def flux_nu(self,tau_nu_xory):
		""" Returns a scalar. 
			Integrates flux(y) over y or flux(x) over x.
			Note: the x or y range by default is the entire 
			simulation box, but the user can specify a narrower 
			range by setting (i_min,i_max) or (j_min,j_max).
		"""
		if self.Ny == 1:
			return np.exp(-tau_nu_xory)
		else: 
			if self.observer_at_90deg:
				m = self.i_min
				n = self.i_max 
			else:
				m = self.j_min
				n = self.j_max
			tau_nu_xory = tau_nu_xory[m:n+1]
			Nvals = 1+n-m
			sorted = np.argsort(tau_nu_xory)
			flux_xory = np.exp(-tau_nu_xory[sorted])
			flux = np.sum(0.5*(flux_xory[0:Nvals-1] + flux_xory[1:Nvals]))
			#flux = np.sum(0.5*(flux_xory[m:n] + flux_xory[m+1:n+1]))
			
			flux /= Nvals
			return flux
			
	def tau_nu_old(self,tau_nu_xory):
		""" Returns a scalar.
			Integrates tau_nu(x) over x or tau_nu(y) over y.
			Note: the x or y-range by default is the entire simulation
			box, but the user can specify a narrower range by
			setting (i_min,i_max) or (j_min,j_max).
		"""
		
		if self.Ny == 1:
			return tau_nu_xory
		else: 
			if self.observer_at_90deg:
				m = self.i_min
				n = self.i_max 
			else:
				m = self.j_min
				n = self.j_max
			tau_nu_xory = tau_nu_xory[m:n+1]
			Nvals = 1+n-m
			sorted = np.argsort(tau_nu_xory)
			tau_nu_xory = tau_nu_xory[sorted]
			tau_nu_avg = np.sum(0.5*(tau_nu_xory[0:Nvals-1] + tau_nu_xory[1:Nvals]))
			#tau_nu_avg = np.sum(0.5*(tau_nu_xory[m:n] + tau_nu_xory[m+1:n+1]))
			if 0:
				print("max(taus) = {}, sum(taus) = {}, sum(taus)/(1+n-m) = {}".format(np.max(tau_nu_xory[0:Nvals]),tau_nu_avg,tau_nu_avg/Nvals))
			tau_nu_avg /= Nvals
			return tau_nu_avg
			
	def tau_nu(self,tau_nu_xory,flux=None):
		""" Returns a scalar.
			Integrates tau_nu(x) over x or tau_nu(y) over y.
			Note: the x or y-range by default is the entire simulation
			box, but the user can specify a narrower range by
			setting (i_min,i_max) or (j_min,j_max).
		"""
		
		if self.Ny == 1:
			return tau_nu_xory
		else: 
			if self.observer_at_90deg:
				m = self.i_min
				n = self.i_max 
			else:
				m = self.j_min
				n = self.j_max
			tau_nu_xory = tau_nu_xory[m:n+1]
			if self.DOWN_AVERAGING:
				#tau_cutoff = np.mean(tau_nu_xory[self.covered_lc])
				tau_nu_xory = tau_nu_xory[self.covered_lc]
				Nvals = len(tau_nu_xory)
				sorted = np.argsort(tau_nu_xory)
				tau_nu_xory = tau_nu_xory[sorted]
				tau_cutoff = np.sum(0.5*(tau_nu_xory[0:Nvals-1] + tau_nu_xory[1:Nvals]))
				tau_cutoff /= Nvals
				if tau_cutoff < self.tau_min:
					tau_cutoff = self.tau_min
			else:
				if flux==None:
					tau_cutoff = self.tau_cutoff
				else:
					tau_cutoff = max(self.tau_cutoff,-np.log(flux))
					print("[tau_nu]: tau_cutoff = {}".format(tau_cutoff))
			covered = np.where(tau_nu_xory > tau_cutoff)
			covered = covered[0]
			tau_nu_xory = tau_nu_xory[covered]
			Nvals = len(tau_nu_xory)
			sorted = np.argsort(tau_nu_xory)
			tau_nu_xory = tau_nu_xory[sorted]
			tau_nu_avg = np.sum(0.5*(tau_nu_xory[0:Nvals-1] + tau_nu_xory[1:Nvals]))
			#tau_nu_avg = np.sum(0.5*(tau_nu_xory[m:n] + tau_nu_xory[m+1:n+1]))
			if 0:
				print("max(taus) = {}, sum(taus) = {}, sum(taus)/(1+n-m) = {}".format(np.max(tau_nu_xory[0:Nvals]),tau_nu_avg,tau_nu_avg/Nvals))
			tau_nu_avg /= Nvals
			return tau_nu_avg
			
	def covering_fraction(self,tau_xory,i_nu,flux=None):
		""" covering fraction at a single frequency """
		if self.Ny == 1: #this assumes a uniform grid
			print("ERROR: covering_fractions not defined in 1D!")
			sys.exit()
		else: 
			if self.observer_at_90deg:
				m = self.i_min
				n = self.i_max 
			else:
				m = self.j_min
				n = self.j_max
			if self.DOWN_AVERAGING:
				tau_cutoff = np.mean(tau_xory[self.covered_lc])
				if tau_cutoff < self.tau_min:
					tau_cutoff = self.tau_min
				print("i_nu = {}, tau_cutoff = {}".format(i_nu,tau_cutoff))
			else:
				if flux == None:
					tau_cutoff = self.tau_cutoff
				else:
					tau_cutoff = max(self.tau_cutoff,-np.log(flux))
			covered = (tau_xory[m:n+1] > tau_cutoff)
			C = np.sum(covered)/np.float(1+n-m)

			return C
			
	def calc_profiles(self,line_sigmas,do_pickle,use_flux):
		""" optical depth profile """
		tau_nus = {}; flux_nus = {}; tau_nu_xorys = {}
		keys = self.lines + ['total']
		for key in keys:
			flux_nus[key] = np.zeros(self.Nv)	
			tau_nus[key] = np.zeros(self.Nv)
			if self.observer_at_90deg:
				tau_nu_xorys[key] = np.zeros((self.Nv,self.Nx))
			else:
				Npts = self.Ny
				tau_nu_xorys[key] = np.zeros((self.Nv,Npts))
				
		if self.observer_at_90deg:
			tau_total = np.zeros((self.Nv,self.Nx))
		else:
			tau_total = np.zeros((self.Nv,self.Ny))
		
		if self.DOWN_AVERAGING: # calculate line center coverage
			phi_arr = self.phi_nu(self.y[self.Nv/2-1])
			alpha_arr = self.alpha_nu(phi_arr,line_sigmas[0])
			if self.observer_at_90deg:
				tau_nu_xory = self.tau_nu_x(alpha_arr)
			else:
				tau_nu_xory = self.tau_nu_y(alpha_arr)
			covered = np.where(tau_nu_xory > self.tau_cutoff)
			self.covered_lc = covered[0]
		
		indx = 0
		for sigma_arr in line_sigmas:
			line = self.lines[indx]
			print("[calc_profiles]:\ncalculating tau-profile for line {}".format(line))
			indx += 1
			for i in range(self.Nv):
				phi_arr = self.phi_nu(self.y[i])
				alpha_arr = self.alpha_nu(phi_arr,sigma_arr)
				if self.observer_at_90deg:
					tau_nu_x = self.tau_nu_x(alpha_arr)
					tau_nu_xorys[line][i] = tau_nu_x
				else:
					tau_nu_y = self.tau_nu_y(alpha_arr)
					tau_nu_xorys[line][i] = tau_nu_y

				flux_nus[line][i] = self.flux_nu(tau_nu_xorys[line][i])
				tau_total[i] += tau_nu_xorys[line][i]
				if indx == len(self.lines):
					tau_nu_xorys['total'][i] = tau_total[i]
					flux_nus['total'][i] = self.flux_nu(tau_total[i])
				if use_flux:
					tau_nus[line][i] = self.tau_nu(tau_nu_xorys[line][i],flux_nus[line][i])
				else:
					tau_nus[line][i] = self.tau_nu(tau_nu_xorys[line][i])
			flux_nus['vel_kms'] = self.vel
			print("Done.")
		
		if do_pickle:
			print("[calc_profiles]:\nPickling tau profiles to {}".format(self.pickled_taus))
			pickle.dump(tau_nus, open(self.pickled_taus, "wb"), protocol=2)
			print("[calc_profiles]:\nPickling flux profiles to {}".format(self.pickled_fluxes))
			pickle.dump(flux_nus, open(self.pickled_fluxes, "wb"), protocol=2)
			if self.Ny > 1 : # don't pickle for 1D runs (b/c tau_nu_xorys = tau_nus)
				print("[calc_profiles]:\nPickling tau_xory profiles to {}".format(self.pickled_tau_xorys))
				pickle.dump(tau_nu_xorys, open(self.pickled_tau_xorys, "wb"), protocol=2)
		
		return tau_nus,flux_nus,tau_nu_xorys
			
	# main calc - previously named make_slab_profiles
	def calc_final_profiles(self,do_pickle=False,recalc=False,use_flux=False):
		# if use_flux=True, tau_cutoff is set to -np.log(flux)
		
		self.tau_profiles = {}
		self.line_profiles = {}
		self.tau_xory_profiles = {}

		if recalc or self.recalc:
			print("[calc_final_profiles]:\nre-calculating optical depth profiles...")
			tau_profiles, flux_profiles, tau_nu_xorys = self.calc_profiles(self.line_sigmas,do_pickle,use_flux)
		elif do_pickle:
			try:
				print("[calc_final_profiles]:\nAttempting to open pickle files...")
				tau_profiles = pickle.load(open(self.pickled_taus,"rb"))
				flux_profiles = pickle.load(open(self.pickled_fluxes,"rb"))
				tau_nu_xorys = pickle.load(open(self.pickled_tau_xorys,"rb"))
			except:
				print("[calc_final_profiles]:\ncalculating optical depth profiles...")
				tau_profiles, flux_profiles, tau_nu_xorys = self.calc_profiles(self.line_sigmas,do_pickle,use_flux)
		
		keys = self.lines + ['total']
		for key in keys:
			self.tau_profiles[key] = tau_profiles[key]
			self.line_profiles[key] = flux_profiles[key]
			self.tau_xory_profiles[key] = tau_nu_xorys[key]

		# final profiles
		self.tau_profiles['total_tau'] = np.zeros_like(self.y)
		for line in self.lines:
			self.tau_profiles['total_tau'] += tau_profiles[line]
		self.line_profiles['total_tau'] = np.exp(-self.tau_profiles['total_tau'])
	
	def equiv_widths(self,flux_vs_nu,v_t):
		width_vs_nu = 1. - flux_vs_nu
		dys = self.y[1:] - self.y[:-1]
		eq_width = np.sum(0.5*(width_vs_nu[0:-1] + width_vs_nu[1:])*dys/v_t) # eq_width in units of thermal width
# 		eq_width = np.sum(0.5*(width_vs_nu[0:-1] + width_vs_nu[1:])*dys)
		return eq_width
	
	def calc_curveofgrowth(self,do_pickle):
		""" optical depth profile """
		cog = {}; cog_normed = {}
		cumtau_nu = {}; tau_nu_vs_r = {}; flux_nu_vs_r = {}
		keys = self.lines + ['total']
		
		Nr = self.Nx
		for key in keys:
			cog[key] = np.zeros(Nr)
			cog_normed[key] = np.zeros(Nr)	
			cumtau_nu[key] = np.zeros(self.Nv)
			tau_nu_vs_r[key] = np.zeros((self.Nv,Nr))
			flux_nu_vs_r[key] = np.zeros((self.Nv,Nr))
			tau_total = np.zeros((self.Nv,Nr))
		
		indx = 0
		for sigma_arr in self.line_sigmas:
			line = self.lines[indx]
			print("[calc_curveofgrowth]:\ncalculating cumulative taus,fluxes for line {}".format(line))
			indx += 1
			for i in range(self.Nv):
				phi_arr = self.phi_nu(self.y[i])
				alpha_arr = self.alpha_nu(phi_arr,sigma_arr)

				cumtau_nu = self.cumtau_nu(alpha_arr)
				tau_nu_vs_r[line][i] = cumtau_nu
				flux_nu_vs_r[line][i] = np.exp(-cumtau_nu)
				tau_total[i] += tau_nu_vs_r[line][i]
				
				if indx == len(self.lines):
					tau_nu_vs_r['total'][i] = tau_total[i]
					flux_nu_vs_r['total'][i] = np.exp(-tau_total[i])
		indx = 0	
		for line in self.lines:
			print("[calc_curveofgrowth]:\ncalculating equivalent widths for line {}".format(line))
			indx += 1
			for i in range(Nr):
				i_min = flux_nu_vs_r[line][:,i].argmin()
				cog[line][i] = self.equiv_widths(flux_nu_vs_r[line][:,i],self.vt[0])
				cog_normed[line][i] = self.equiv_widths(flux_nu_vs_r[line][:,i],self.vt[i])
				print(self.vt[i],tau_nu_vs_r[line][:,i][i_min],flux_nu_vs_r[line][:,i][i_min],cog[line][i])
				
				if indx == len(self.lines):
					cog['total'][i] = self.equiv_widths(flux_nu_vs_r['total'][:,i],self.vt[0])
					cog_normed['total'][i] = self.equiv_widths(flux_nu_vs_r['total'][:,i],self.vt[i])
		
		if do_pickle:
			print("[calc_curveofgrowth]:\nPickling equivalent width profiles to {}".format(self.pickled_cog))
			pickle.dump(cog, open(self.pickled_cog, "wb"), protocol=2)
			pickle.dump(cog_normed, open(self.pickled_cog_normed, "wb"), protocol=2)
		
		return cog,cog_normed
		
	def get_LOS_opacities(self,indices,y_dir=False):
		
		los_opacities_red = []
		los_opacities_blue = []
		ctr = 0
		for sigmas in self.line_sigmas:
			if ctr == 0:
				los_opacities = los_opacities_blue
			else:
				los_opacities = los_opacities_red
			ctr += 1
			for idx in indices:
				print("idx = {}".format(idx))
				print(sigmas[0,:].shape)
				print(sigmas[:,0].shape)
				if self.observer_at_90deg and not(y_dir):
					print("y = {}".format(self.yarr[0][idx]))
					los_opacities.append(sigmas[:,idx])
				else:
					print("x = {}".format(self.xarr[idx,0]))
					los_opacities.append(sigmas[idx,:])
		
		return los_opacities_blue,los_opacities_red
	
	def get_LOS_opacities1D(self):
		
		los_opacities_red = []
		los_opacities_blue = []
		ctr = 0
		for sigmas in self.line_sigmas:
			if ctr == 0:
				los_opacities = los_opacities_blue
			else:
				los_opacities = los_opacities_red
			ctr += 1

			los_opacities.append(sigmas)
		
		return los_opacities_blue,los_opacities_red
			
	def calc_covering_fractions_exact(self,Intensities=None):
		""" covering factor profile """
		Cs = {}
		for line in self.lines:
			Cs[line] = np.zeros(self.Nv)
		
		if self.DOWN_AVERAGING: # calculate line center coverage
			phi_arr = self.phi_nu(self.y[self.Nv/2-1])
			alpha_arr = self.alpha_nu(phi_arr,line_sigmas[0])
			if self.observer_at_90deg:
				tau_nu_xory = self.tau_nu_x(alpha_arr)
			else:
				tau_nu_xory = self.tau_nu_y(alpha_arr)
			covered = np.where(tau_nu_xory > self.tau_cutoff)
			self.covered_lc = covered[0]
			
		indx = 0
		for sigmaofx in self.line_sigmas:
			line = self.lines[indx]
			print("calculating exact covering fraction for line {}".format(line))
			if Intensities != None:
				I = Intensities[indx]
			indx += 1
			for i in range(self.Nv):
				tau_xy = self.tau_xory_profiles[line][i]
				if Intensities == None:
					C = self.covering_fraction(tau_xy,i)
				else:
					C = self.covering_fraction(tau_xy,i,flux = I[i])
				Cs[line][i] = C
			print("Done.")
		
		return Cs
		
	def calc_ionfraction_maps(self,fs,i_y=None):
		""" calculates f_ion """
		Maps = {}
		for line in self.lines:
			Maps[line] = np.zeros((self.Ny,self.Nx))
		
		# choose line center as default
		if i_y == None:
			i_y = self.Nv/2-1
			
		indx = 0
		for sigma in self.line_sigmas:
			print("sigma_max = {}, sigma_min = {}".format(np.max(sigma),np.min(sigma)))
			line = self.lines[indx]
			a12 = (np.pi*Fund['e']**2/Fund['me']/Fund['c'])*fs[indx]
			print("calculating abundance map for line {}".format(line))
			indx += 1
			Maps[line] = sigma*self.line_data[line]['nu0']*self.vt/(a12*self.n_xstar)/self.A_Z
			print("Done.")
		
		return Maps
		
	def calc_AMD(self):
		n = self.n_cgs
		n_ip1 = np.zeros((self.Ny,self.Nx))
		n_ip1[:,:-1] = self.n_cgs[:,1:]
		n_ip1[:,-1] = self.n_cgs[:,0]
		dn_by_dx = (n_ip1 - n)/self.dx
		n_sqd = n**2
		print("calc_AMD: np.min(n^2), np.max(n^2) = {},{}".format(np.min(n_sqd),np.max(n_sqd)))
		print("calc_AMD: np.min(dn/dx), np.max(dn/dx) = {},{}".format(np.min(dn_by_dx),np.max(dn_by_dx)))
		
		AMD = -n_sqd/dn_by_dx
		# set places where dn/dx = 0 to their neighbors
		ii,jj = np.where(dn_by_dx == 0.)
		for i in range(len(ii)):
			AMD[ii[i],jj[i]] = AMD[ii[i],jj[i]+1] 
		return AMD,dn_by_dx
		
	def calc_covering_fraction_PPC(self,R):
		""" covering factor profile """
		Cs = np.zeros(self.Nv)
		
		Ib = self.line_profiles[self.lines[0]]
		Ir = self.line_profiles[self.lines[1]]
		for i in range(self.Nv):
			if R == 0.5: # analytic case
				N = Ir[i]**2 - 2.*Ir[i] + 1.
				D = Ib[i]    - 2.*Ir[i] + 1.
				C = N/D
				print("PPC calc: i = {}, 1/R = {}, N = {}, D = {}, C = {}".format(i,1./R,N,D,C))
			else:
				f = lambda C: (Ib[i] + C - 1.)/C - ((Ir[i] + C - 1.)/C)**(1./R)
				C = brentq(f,0.,1.,rtol=1e-20,maxiter=500,full_output=False,disp=True)
				#C = newton(f, 15., fprime=None, args=(), tol=1e-13, maxiter=50,fprime2=None)
				print("PPC calc: i = {}, 1/R = {}, C = {}".format(i,1./R,C))
			Cs[i] = C
	
		return Cs
		
		
				
		
		
