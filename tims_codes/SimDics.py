import athena_read   # Athena++'s reader
import my_athena_reader	 # Randall's custom Athena++ tab reader
import numpy as np
import sys
import code_units as units

def my_var(grid_dic,d,p,v_r,v_t,v_p,B_r=None,B_t=None,B_p=None):
	if 0:
		Bp = np.sqrt(B_r**2 + B_t**2)
		v_Ap = Bp/np.sqrt(d)
		vp = np.sqrt(v_r**2 + v_t**2)
		qty = v_Ap/vp
	if 0:
		B = np.sqrt(B_r**2 + B_t**2 + B_p**2)
		v_A = B/np.sqrt(d)
		qty = v_A/v_p
	if 0:
		qty = grid_dic['theta_c_grid']
	if 0:
		B_x = B_r*grid_dic['sin_thetas'] + B_t*grid_dic['cos_thetas']
		B_z = B_r*grid_dic['cos_thetas'] - B_t*grid_dic['sin_thetas']
		theta_field = np.arctan(B_z/B_x)*180./np.pi
		qty = 90. - theta_field
		qty = theta_field
	if 0:
		thetas = np.arctan(grid_dic['sin_thetas']/grid_dic['cos_thetas'])
		theta_field = (thetas + np.arctan(B_t/B_r))*180./np.pi
		qty = 90 - theta_field*grid_dic['r_c']
		qty = np.arctan(B_t/B_r)*180./np.pi
	if 0:
		B_x = B_r*grid_dic['sin_thetas'] + B_t*grid_dic['cos_thetas']
		B_z = B_r*grid_dic['cos_thetas'] - B_t*grid_dic['sin_thetas']
		Bx_hat = B_x/np.sqrt(B_x**2 + B_z**2)
		Bz_hat = B_z/np.sqrt(B_x**2 + B_z**2)
		x_hat = grid_dic['x_c']/np.sqrt(grid_dic['x_c']**2 + grid_dic['y_c']**2)
		z_hat = grid_dic['y_c']/np.sqrt(grid_dic['x_c']**2 + grid_dic['y_c']**2)
# 		qty = np.arccos(Bx_hat*x_hat + Bz_hat*x_hat)*180./np.pi
		qty = np.arccos(Bx_hat*x_hat)*180./np.pi
	if 0:
		Br_hat = B_r/np.sqrt(B_r**2 + B_t**2)
		Bt_hat = B_t/np.sqrt(B_r**2 + B_t**2)
		x_hat = grid_dic['x_c']/np.sqrt(grid_dic['x_c']**2 + grid_dic['y_c']**2)
# 		z_hat = grid_dic['y_c']/np.sqrt(grid_dic['x_c']**2 + grid_dic['y_c']**2)
		qty = np.arccos(Br_hat*x_hat + Bt_hat*x_hat)*180./np.pi
	if 0: # used for the WP12 streamline_angles.py benchmark
		v_x = v_r*grid_dic['sin_thetas'] + v_t*grid_dic['cos_thetas']
		v_z = v_r*grid_dic['cos_thetas'] - v_t*grid_dic['sin_thetas']
		nhat_dot_xhat = 1./np.sqrt(1. + (v_z/v_x)**2)
		qty = np.arccos(nhat_dot_xhat)*180./np.pi
	if 1:
		B_x = B_r*grid_dic['sin_thetas'] + B_t*grid_dic['cos_thetas']
		B_z = B_r*grid_dic['cos_thetas'] - B_t*grid_dic['sin_thetas']
		nhat_dot_xhat = 1./np.sqrt(1. + (B_z/B_x)**2)
		qty = np.arccos(nhat_dot_xhat)*180./np.pi
	if 0:
		vp = np.sqrt(v_r**2 + v_t**2)
		mfd = d*vp
		grad_mfd_r,grad_mfd_t = calc_grad_r(mfd,grid_dic['r_c_grid']), calc_grad_t(mfd,grid_dic['theta_c_grid'])
		n_r,n_t = v_r/vp,v_t/vp
		dlnA = -(n_r*grad_mfd_r + n_t*grad_mfd_t)/mfd
		qty = dlnA
	return qty

def my_var_label(logval):
	if logval:
# 		return r'$\log{(v_A/v_\phi)}$'
		return r'$\log{(\theta_B)}$'
	else:
# 		return r'$v_A/v_\phi$'
		return r'$\theta_B$'


def uov_scalar(qty,grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p):
	""" 
	This function is to be defined by the user to 
	perform any manipulations of scalar user output 
	variables. The function uov_vector() calls this 
	function upon loading each uov.

	Below is an example of performing some 
	manipulations assuming index=0 and qty=div(B).  
	"""

	if 0: # div(B) calc - change to 0 to ignore
	# For analyzing divergence of B, we likely want
	# to take a log, so here we set all zero values to 
	# a small number.  We also normalize div(B) by B.
		filter = np.where(qty == 0.)
		qty[filter] = 1e-25

		# normalize div B by B
		B = np.sqrt(B_r**2 + B_t**2 + B_p**2)
		qty = qty #/B
	
	else: # user can add modify as needed 
		qty = qty

	return qty


def uov_vector(uovs,grid_dic,d,p,v_r,v_t,v_p,B_r=None,B_t=None,B_p=None):
	""" 
	This function is to be defined by the user to 
	perform any manipulations of vector user output 
	variables. The function load_vars() calls this 
	function after adding all vector components 
	(which are individual uov's) to an array, which 
	is passed here as the input "uovs".

	A scalar uov will be passed along to uov_scalar()

	Currently this function is set to expect the
	components of current density, j = curl(B), in
	spherical coordinates.  The (r,theta,phi) 
	components of j correspond to indices (0,1,2).
	Below we further take the cross-product with
	B to arrive at the Lorentz force,f_L = jxB.  
	"""

	N = len(uovs)
	print("[uov_vector]: N = {}".format(N))
	j = {}
	for n in range(N):
		key = np.str(n+1)
		j[key] = np.array(uovs[n])

	if N is 1: # not a vector; call scalar routine instead
		qty = uov_scalar(uovs[0],grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p)
		return [qty]
	
	# now we compute jxB
	elif N is 2:
		j_r,j_t = j['1'],j['2']
		jcrossB_p = j_r*B_t - j_t*B_r
		return [jcrossB_p]

	elif N is 3:
		j_r,j_t,j_p = j['1'],j['2'],j['3']
		jcrossB_r = j_t*B_p - j_p*B_t
		jcrossB_t = j_r*B_p - j_p*B_r
		jcrossB_p = j_r*B_t - j_t*B_r
		return [jcrossB_r, jcrossB_t, jcrossB_p]
# 		return np.sqrt(j_r**2 + j_t**2 + j_p**2)
	
	# user can add modify as needed 
	else: # if N > 3 
		uov_array = []
		for key in j.keys():
			uov_array.append(j[key])
		return uov_array



def user_out_var_label(logval):
	if logval:
		return r'$j$' #r'$\log{((\nabla \cdot B)/B)}$'
	else:
		return r'$div(B)/B$'


def dump_dir_dic(run_names,run_paths,dump_tags,uov_tags=None):
	"""
	inuts:
	run_names - always an array of strings specifying the dump folder
	run_paths - either a string or a dictionary of strings
			  - the strings are names of the paths to the dump folder
	dump_tags - either a string or a dictionary of strings
			  - the strings are names of the actual output dumps
	uov_tags  - (optional) either a string or a dictionary of strings
			  - dump_tags of "user output variables" (see Athena++ docs)

	output:
	dump_dirs - a dictionary containing the complete paths to individual dumps
				for all specified simulations 
	"""
	dump_dirs = {}
	
	Nsims = len(run_names)
	
	# paths to all the sim directories 
	if isinstance(run_paths, dict):
		for n in range(1,Nsims+1):
			key = np.str(n)
			if key in run_paths.keys():
				dump_dirs[key+'path'] = run_paths[key] + run_names[n-1] + '/'
			else:
				dump_dirs[key+'path'] = run_paths['1'] + run_names[n-1] + '/'
				print("[dump_dir_dic]: n = {}, path = {}".format(n,dump_dirs[key+'path']))
	else: # run_paths is a string
		for n in range(1,Nsims+1):
			dump_dirs[np.str(n)+'path'] = run_paths + run_names[n-1] + '/'
	
	# path to the main data files
	if isinstance(dump_tags, dict):
		for n in range(1,Nsims+1):
			key = np.str(n)
			if key in dump_tags.keys():
				dump_dirs[key] = dump_dirs[key+'path'] + dump_tags[key]
			else:
				dump_dirs[key] = dump_dirs[key+'path'] + dump_tags['1']
	else: # dump_tags is a string not a dic
		for n in range(1,Nsims+1):
			dump_dirs[np.str(n)] = dump_dirs[np.str(n)+'path'] + dump_tags

	# path to user output data files
	if uov_tags is not None: 
		if isinstance(uov_tags, dict):
			for n in range(1,Nsims+1):
				key = np.str(n) + 'uov'
				if np.str(n) in dump_tags.keys():
					dump_dirs[key] = dump_dirs[np.str(n)+'path'] + uov_tags[np.str(n)]
				else:
					dump_dirs[key] = dump_dirs[np.str(n)+'path'] + uov_tags['1']
		else: # uov_tags is a string not a dic
			for n in range(1,Nsims+1):
				dump_dirs[np.str(n)+'uov'] = dump_dirs[np.str(n)+'path'] + uov_tags
			
	return dump_dirs

def load_labels1D(logval=False):
	""" Returns a dictionary of latex plot labels for all variables 
		currently implemented
	
	Inputs:
	logval - if True, returns r'$\log(qty)$' instead of just r'$qty$'

	Example usage:
	vars = ['B']
	SimDic = load_vars(i_dump,vars,dump_dirs,grid_dics)
	LabelDic = load_labels(logval)
	plot(r_c,qty,label=LabelDic['B'])
	"""
	
	LabelDic = {}	
	hydro_vars = {}
	hc_vars = {}
	misc_vars = {}
	mhd_vars = {}
	
	if logval:
		# labels for primitive hydro variables 
		hydro_vars['d'] = r'$\log(\rho\:[\rm{g\,cm^{-3}}])$'
		hydro_vars['p'] = r'$\log(p)$'
		hydro_vars['v'] = r'$\log(v)$'
		
		# labels for derived hydro variables 
		hydro_vars['v'] = r'$\log(v)$'
		hydro_vars['vp'] = r'$\log(v_{poloidal})$'
		hydro_vars['mfd'] = r'$\log(\rho v_{poloidal})$'
		hydro_vars['cs'] = r'$\log(c_s)$'
		hydro_vars['c_iso'] = r'$\log(c_{iso})$'
		hydro_vars['s'] = r'$s=\log(s)$' 
		hydro_vars['T'] = r'$\log(T \:[\rm{K}])$'
		hydro_vars['M'] = r'$\log(M)$'
		hydro_vars['Be'] = r'$\log(Be)$'
		hydro_vars['v_sqd'] = r'$\log(v^2)$'
		hydro_vars['KE'] = r'$\log[(1/2)\rho v^2]$'
		hydro_vars['TE'] = r'$\log(E)$'

		hc_vars['rhoL'] = r'$\log(\mathcal{\rho L})$'
		hc_vars['netL'] = r'$\log(\mathcal{L})$'
		hc_vars['netC'] = r'$\log(\mathcal{C})$'
		hc_vars['netH'] = r'$\log(\mathcal{H})$'
		hc_vars['netLoverT'] = r'$\log(\mathcal{L}/T)$'
		hc_vars['netLoverv'] = r'$\log(\mathcal{L}/v)$'
		hc_vars['div_q'] = r'$\log(\nabla\cdot q)$'
		hc_vars['gradBe'] = r'$\log(dBe/dx)$'
		hc_vars['grad_s'] = r'$\log(ds/dx)$'


		# labels for MHD quantities
		mhd_vars['B_r'] = r'$\log(B_r)$'
		mhd_vars['B_t'] = r'$\log(B_\theta)$'
		mhd_vars['B_p'] = r'$\log(B_\phi)$'
		mhd_vars['B'] = r'$\log(B)$'
		mhd_vars['Bp'] = r'$\log(B_{poloidal})$'
		mhd_vars['E_r'] = r'$\log(E_r)$'
		mhd_vars['E_t'] = r'$\log(E_\theta)$'
		mhd_vars['E_p'] = r'$\log(E_\phi)$'
		mhd_vars['E'] = r'$\log(E)$'
		mhd_vars['Ep'] = r'$\log(E_{poloidal})$'
		mhd_vars['E_p_over_Ep'] = r'$\log(E_\phi/E_{poloidal})$'
		mhd_vars['beta'] = r'$\log(\beta)$'
		mhd_vars['v_A'] = r'$\log(v_{A})$'
		mhd_vars['v_Ap'] = r'$\log(v_{A})$'
		mhd_vars['v_fast'] = r'$\log(v_{fast})$'
		mhd_vars['v_slow'] = r'$\log(v_{slow})$'
		mhd_vars['Bphi_over_Bp'] = r'$\log(B_\phi/B_{poloidal})$'
		mhd_vars['B2'] = r'$\log{(B^2)}$'
		mhd_vars['A'] = r'$\log{A}$'
		mhd_vars['j_r'] = r'$\log[j_r = (\nabla\times B)_r]$'
		mhd_vars['j_t'] = r'$\log[j_\theta = (\nabla\times B)_t]$'
		mhd_vars['j_p'] = r'$\log[j_\phi = (\nabla\times B)_p]$'
		mhd_vars['jp'] = r'$\log[j_p = (\nabla\times B)_{poloidal}]$'
		mhd_vars['j'] = r'$\log[j = \nabla\times B]$'
		mhd_vars['jcrossB_r'] = r'$\log[j_r = (j\times B)_r]$'
		mhd_vars['jcrossB_t'] = r'$\log[j_\theta = (j\times B)_t]$'
		mhd_vars['jcrossB_p'] = r'$\log[j_\phi = (j\times B)_p]$'
		mhd_vars['jcrossBp'] = r'$\log[j_p = (j\times B)_{poloidal}]$'
		mhd_vars['jcrossB'] = r'$\log[j = j\times B]$'
		mhd_vars['omega_A'] = r'$\log[omega_A]$'
		mhd_vars['Fmag_par'] = r'$\log(F_\parallel)$'
		mhd_vars['Fmag_phi'] = r'$\log(F_\phi)$'
		mhd_vars['gradP_par'] = r'$\log[(\nabla p)_\parallel]$'
		mhd_vars['Be_mag'] = r'$\log(Be_{mag})$'


		# label for my_var, user_out_var, & floor maps
		misc_vars['tau'] = r'$\log(\tau)$'
		misc_vars['xi'] = r'$\log(\xi)$'
		misc_vars['dlnA'] = r'$\log[-\hat{n}\cdot\nabla(\rho v_{poloidal})/\rho v_{poloidal}]$'
		misc_vars['div_v'] = r'$\log(\nabla \cdot v)]$'
		misc_vars['my_var'] = my_var_label(logval)
		misc_vars['user_out_var'] = user_out_var_label(logval)
		misc_vars['d_floor'] = r'$\rho_{floor}$' + ' triggers'
		misc_vars['p_floor'] = r'$p_{floor}$' + ' triggers'
		misc_vars['dt_cs_r'] = r'$\log(dt_{cs,r})$'
		misc_vars['dt_cs_t'] = r'$\log(dt_{cs,\theta})$'
		misc_vars['dt_v_r'] = r'$\log(dt_{v,r})$'
		misc_vars['dt_v_t'] = r'$\log(dt_{v,\theta})$'
		misc_vars['dt_vA_r'] = r'$\log(dt_{v_A,r})$'
		misc_vars['dt_vA_t'] = r'$\log(dt_{v_A,\theta})$'

	else:
		# labels for primitive hydro variables 
		hydro_vars['d'] = r'$\rho$'
		hydro_vars['p'] = r'$p$'
		hydro_vars['v_r'] = r'$v_r$'
		hydro_vars['v_t'] = r'$v_t$'
		hydro_vars['v_p'] = r'$v_\phi/v_{kep}$'
		
		# labels for derived hydro variables 
		hydro_vars['v'] = r'$v$'
		hydro_vars['vp'] = r'$v_{poloidal}$'
		hydro_vars['mfd'] = r'$\rho v_{poloidal}$'
		hydro_vars['cs'] = r'$c_s$'
		hydro_vars['c_iso'] = r'$c_{iso}$'
		hydro_vars['s'] = r'$s=\log(p/\rho^\gamma)$' 
		hydro_vars['T'] = r'$T/T_C$'
		hydro_vars['T_vir'] = r'$T/T_{vir}$'
		hydro_vars['M'] = r'$M$'
		hydro_vars['Be'] = r'$Be$'
		hydro_vars['v_sqd'] = r'$v^2$'
		hydro_vars['KE'] = r'$(1/2)\rho v^2$'
		hydro_vars['TE'] = r'$E$'

		hc_vars['rhoL'] = r'$\mathcal{\rho L}$'
		hc_vars['netL'] = r'$\mathcal{L}$'
		hc_vars['netC'] = r'$\mathcal{C}$'
		hc_vars['netH'] = r'$\mathcal{H}$'
		hc_vars['netLoverT'] = r'$\mathcal{L}/T$'
		hc_vars['netLoverv'] = r'$\mathcal{L}/v$'
		hc_vars['div_q'] = r'$\nabla\cdot q$'
		hc_vars['gradBe'] = r'$dBe/dx$'
		hc_vars['grad_s'] = r'$ds/dx$'

		# labels for MHD quantities
		mhd_vars['B_r'] = r'$B_r$'
		mhd_vars['B_t'] = r'$B_\theta$'
		mhd_vars['B_p'] = r'$B_\phi$'
		mhd_vars['B'] = r'$B$'
		mhd_vars['Bp'] = r'$B_{poloidal}$'
		mhd_vars['E_r'] = r'$E_r$'
		mhd_vars['E_t'] = r'$E_\theta$'
		mhd_vars['E_p'] = r'$E_\phi$'
		mhd_vars['E'] = r'$E$'
		mhd_vars['Ep'] = r'$E_{poloidal}$'
		mhd_vars['E_p_over_Ep'] = r'$E_\phi/E_{poloidal}$'
		mhd_vars['beta'] = r'$\beta$'
		mhd_vars['v_A'] = r'$v_{A}$'
		mhd_vars['v_Ap'] = r'$v_{A}$'
		mhd_vars['v_fast'] = r'$v_{fast}$'
		mhd_vars['v_slow'] = r'$v_{slow}$'
		mhd_vars['Bphi_over_Bp'] = r'$B_\phi/B_{poloidal}$'
		mhd_vars['B2'] = r'$B^2$'
		mhd_vars['A'] = r'$A$'
		mhd_vars['j_r'] = r'$j_r = (\nabla\times B)_r$'
		mhd_vars['j_t'] = r'$j_\theta = (\nabla\times B)_t$'
		mhd_vars['j_p'] = r'$j_\phi = (\nabla\times B)_p$'
		mhd_vars['jp'] = r'$j_p = (\nabla\times B)_{poloidal}$'
		mhd_vars['j'] = r'$j = \nabla\times B$'
		mhd_vars['jcrossB_r'] = r'$j_r = (j\times B)_r$'
		mhd_vars['jcrossB_t'] = r'$j_\theta = (j\times B)_t$'
		mhd_vars['jcrossB_p'] = r'$j_\phi = (j\times B)_p$'
		mhd_vars['jcrossBp'] = r'$j_p = (j\times B)_{poloidal}$'
		mhd_vars['jcrossB'] = r'$j = j\times B$'
		mhd_vars['omega_A'] = r'$omega_A$'
		mhd_vars['Fmag_par'] = r'$F_\parallel$'
		mhd_vars['Fmag_phi'] = r'$F_\phi$'
		mhd_vars['gradP_par'] = r'$(\nabla p)_\parallel$'
		mhd_vars['Be_mag'] = r'$Be_{mag}$'

		# label for my_var, user_out_var, floor maps, & dt maps
		misc_vars['tau'] = r'$\tau$'
		misc_vars['xi'] = r'$\xi$'
		misc_vars['dlnA'] = r'$-\hat{n}\cdot\nabla(\rho v_{poloidal})/\rho v_{poloidal}$'
		misc_vars['div_v'] = r'$\nabla \cdot v$'
		misc_vars['my_var'] = my_var_label(0)
		misc_vars['user_out_var'] = user_out_var_label(0)
		misc_vars['d_floor'] = r'$\rho_{floor}$' + ' triggers'
		misc_vars['p_floor'] = r'$p_{floor}$' + ' triggers'
		misc_vars['dt_cs_r'] = r'$dt_{cs,r}$'
		misc_vars['dt_cs_t'] = r'$dt_{cs,\theta}$'
		misc_vars['dt_v_r'] = r'$dt_{v,r}$'
		misc_vars['dt_v_t'] = r'$dt_{v,\theta}$'
		misc_vars['dt_vA_r'] = r'$dt_{v_A,r}$'
		misc_vars['dt_vA_t'] = r'$dt_{v_A,\theta}$'


	# add all the labels to LabelDic
	LabelDic['hydro_keys'] = hydro_vars.keys()
	LabelDic['hc_keys'] = hc_vars.keys()
	LabelDic['vorticity_keys'] = ['w_r','w_t','w_p','wp','w']
	LabelDic['misc_keys'] = misc_vars.keys()
	LabelDic['mhd_keys'] = mhd_vars.keys()
	for key in LabelDic['hydro_keys']:
		LabelDic[key] = hydro_vars[key]
	for key in LabelDic['hc_keys']:
		LabelDic[key] = hc_vars[key]
	for key in LabelDic['misc_keys']:
		LabelDic[key] = misc_vars[key]
	for key in LabelDic['mhd_keys']:
		LabelDic[key] = mhd_vars[key]
	
	return LabelDic

def load_labels(logval=False):
	""" Returns a dictionary of latex plot labels for all variables 
		currently implemented
	
	Inputs:
	logval - if True, returns r'$\log(qty)$' instead of just r'$qty$'

	Example usage:
	vars = ['B']
	SimDic = load_vars(i_dump,vars,dump_dirs,grid_dics)
	LabelDic = load_labels(logval)
	plot(r_c,qty,label=LabelDic['B'])
	"""
	
	LabelDic = {}	
	hydro_vars = {}
	misc_vars = {}
	mhd_vars = {}
	
	if logval:
		# labels for primitive hydro variables 
		hydro_vars['d'] = r'$\log(\rho\:[\rm{g\,cm^{-3}}])$'
		hydro_vars['p'] = r'$\log(p)$'
		hydro_vars['v_r'] = r'$\log(v_r)$'
		hydro_vars['v_t'] = r'$\log(v_t)$'
		hydro_vars['v_p'] = r'$\log(v_\phi/v_{kep})$'
		
		# labels for derived hydro variables 
		hydro_vars['v'] = r'$\log(v)$'
		hydro_vars['vp'] = r'$\log(v_{poloidal})$'
		hydro_vars['mfd'] = r'$\log(\rho v_{poloidal})$'
		hydro_vars['cs'] = r'$\log(c_s)$'
		hydro_vars['c_iso'] = r'$\log(c_{iso})$'
		hydro_vars['s'] = r'$s=\log(s)$' 
		hydro_vars['T'] = r'$\log(T \:[\rm{K}])$'
		hydro_vars['T_vir'] = r'$\log(T/T_{vir})$'
		hydro_vars['M'] = r'$\log(M)$'
		hydro_vars['Mp'] = r'$\log(M_{poloidal})$'
		hydro_vars['l'] = r'$\log(l)$'
		hydro_vars['omega'] = r'$\log(\Omega)$'
		hydro_vars['Be'] = r'$\log(Be)$'
		hydro_vars['w_r'] = r'$\log[(\nabla\times v)_r]$'
		hydro_vars['w_t'] = r'$\log[(\nabla\times v)_t]$'
		hydro_vars['w_p'] = r'$\log[(\nabla\times v)_p]$'
		hydro_vars['wp'] = r'$\log[(\nabla\times v)_{poloidal}]$'
		hydro_vars['w'] = r'$\log[\nabla\times v]$'

		# labels for MHD quantities
		mhd_vars['B_r'] = r'$\log(B_r)$'
		mhd_vars['B_t'] = r'$\log(B_\theta)$'
		mhd_vars['B_p'] = r'$\log(B_\phi)$'
		mhd_vars['B'] = r'$\log(B)$'
		mhd_vars['Bp'] = r'$\log(B_{poloidal})$'
		mhd_vars['E_r'] = r'$\log(E_r)$'
		mhd_vars['E_t'] = r'$\log(E_\theta)$'
		mhd_vars['E_p'] = r'$\log(E_\phi)$'
		mhd_vars['E'] = r'$\log(E)$'
		mhd_vars['Ep'] = r'$\log(E_{poloidal})$'
		mhd_vars['E_p_over_Ep'] = r'$\log(E_\phi/E_{poloidal})$'
		mhd_vars['beta'] = r'$\log(\beta)$'
		mhd_vars['v_A'] = r'$\log(v_{A})$'
		mhd_vars['v_Ap'] = r'$\log(v_{A})$'
		mhd_vars['v_fast'] = r'$\log(v_{fast})$'
		mhd_vars['v_slow'] = r'$\log(v_{slow})$'
		mhd_vars['Bphi_over_Bp'] = r'$\log(B_\phi/B_{poloidal})$'
		mhd_vars['B2'] = r'$\log{(B^2)}$'
		mhd_vars['A'] = r'$\log{A}$'
		mhd_vars['j_r'] = r'$\log[j_r = (\nabla\times B)_r]$'
		mhd_vars['j_t'] = r'$\log[j_\theta = (\nabla\times B)_t]$'
		mhd_vars['j_p'] = r'$\log[j_\phi = (\nabla\times B)_p]$'
		mhd_vars['jp'] = r'$\log[j_p = (\nabla\times B)_{poloidal}]$'
		mhd_vars['j'] = r'$\log[j = \nabla\times B]$'
		mhd_vars['jcrossB_r'] = r'$\log[j_r = (j\times B)_r]$'
		mhd_vars['jcrossB_t'] = r'$\log[j_\theta = (j\times B)_t]$'
		mhd_vars['jcrossB_p'] = r'$\log[j_\phi = (j\times B)_p]$'
		mhd_vars['jcrossBp'] = r'$\log[j_p = (j\times B)_{poloidal}]$'
		mhd_vars['jcrossB'] = r'$\log[j = j\times B]$'
		mhd_vars['omega_A'] = r'$\log[omega_A]$'
		mhd_vars['Fmag_par'] = r'$\log(F_\parallel)$'
		mhd_vars['Fmag_phi'] = r'$\log(F_\phi)$'
		mhd_vars['gradP_par'] = r'$\log[(\nabla p)_\parallel]$'
		mhd_vars['Be_mag'] = r'$\log(Be_{mag})$'


		# label for my_var, user_out_var, & floor maps
		misc_vars['tau'] = r'$\log(\tau)$'
		misc_vars['xi'] = r'$\log(\xi)$'
		misc_vars['dlnA'] = r'$\log[-\hat{n}\cdot\nabla(\rho v_{poloidal})/\rho v_{poloidal}]$'
		misc_vars['div_v'] = r'$\log(\nabla \cdot v)]$'
		misc_vars['my_var'] = my_var_label(logval)
		misc_vars['user_out_var'] = user_out_var_label(logval)
		misc_vars['d_floor'] = r'$\rho_{floor}$' + ' triggers'
		misc_vars['p_floor'] = r'$p_{floor}$' + ' triggers'
		misc_vars['dt_cs_r'] = r'$\log(dt_{cs,r})$'
		misc_vars['dt_cs_t'] = r'$\log(dt_{cs,\theta})$'
		misc_vars['dt_v_r'] = r'$\log(dt_{v,r})$'
		misc_vars['dt_v_t'] = r'$\log(dt_{v,\theta})$'
		misc_vars['dt_vA_r'] = r'$\log(dt_{v_A,r})$'
		misc_vars['dt_vA_t'] = r'$\log(dt_{v_A,\theta})$'

	else:
		# labels for primitive hydro variables 
		hydro_vars['d'] = r'$\rho$'
		hydro_vars['p'] = r'$p$'
		hydro_vars['v_r'] = r'$v_r$'
		hydro_vars['v_t'] = r'$v_t$'
		hydro_vars['v_p'] = r'$v_\phi/v_{kep}$'
		
		# labels for derived hydro variables 
		hydro_vars['v'] = r'$v$'
		hydro_vars['vp'] = r'$v_{poloidal}$'
		hydro_vars['mfd'] = r'$\rho v_{poloidal}$'
		hydro_vars['cs'] = r'$c_s$'
		hydro_vars['c_iso'] = r'$c_{iso}$'
		hydro_vars['s'] = r'$s=\log(p/\rho^\gamma)$' 
		hydro_vars['T'] = r'$T/T_C$'
		hydro_vars['T_vir'] = r'$T/T_{vir}$'
		hydro_vars['M'] = r'$M$'
		hydro_vars['Mp'] = r'$M_{poloidal}$'
		hydro_vars['l'] = r'$l$'
		hydro_vars['omega'] = r'$\Omega$'
		hydro_vars['Be'] = r'$Be$'
		hydro_vars['w_r'] = r'$(\nabla\times v)_r$'
		hydro_vars['w_t'] = r'$(\nabla\times v)_t$'
		hydro_vars['w_p'] = r'$(\nabla\times v)_p$'
		hydro_vars['wp'] = r'$(\nabla\times v)_{poloidal}$'
		hydro_vars['w'] = r'$\nabla\times v$'

		# labels for MHD quantities
		mhd_vars['B_r'] = r'$B_r$'
		mhd_vars['B_t'] = r'$B_\theta$'
		mhd_vars['B_p'] = r'$B_\phi$'
		mhd_vars['B'] = r'$B$'
		mhd_vars['Bp'] = r'$B_{poloidal}$'
		mhd_vars['E_r'] = r'$E_r$'
		mhd_vars['E_t'] = r'$E_\theta$'
		mhd_vars['E_p'] = r'$E_\phi$'
		mhd_vars['E'] = r'$E$'
		mhd_vars['Ep'] = r'$E_{poloidal}$'
		mhd_vars['E_p_over_Ep'] = r'$E_\phi/E_{poloidal}$'
		mhd_vars['beta'] = r'$\beta$'
		mhd_vars['v_A'] = r'$v_{A}$'
		mhd_vars['v_Ap'] = r'$v_{A}$'
		mhd_vars['v_fast'] = r'$v_{fast}$'
		mhd_vars['v_slow'] = r'$v_{slow}$'
		mhd_vars['Bphi_over_Bp'] = r'$B_\phi/B_{poloidal}$'
		mhd_vars['B2'] = r'$B^2$'
		mhd_vars['A'] = r'$A$'
		mhd_vars['j_r'] = r'$j_r = (\nabla\times B)_r$'
		mhd_vars['j_t'] = r'$j_\theta = (\nabla\times B)_t$'
		mhd_vars['j_p'] = r'$j_\phi = (\nabla\times B)_p$'
		mhd_vars['jp'] = r'$j_p = (\nabla\times B)_{poloidal}$'
		mhd_vars['j'] = r'$j = \nabla\times B$'
		mhd_vars['jcrossB_r'] = r'$j_r = (j\times B)_r$'
		mhd_vars['jcrossB_t'] = r'$j_\theta = (j\times B)_t$'
		mhd_vars['jcrossB_p'] = r'$j_\phi = (j\times B)_p$'
		mhd_vars['jcrossBp'] = r'$j_p = (j\times B)_{poloidal}$'
		mhd_vars['jcrossB'] = r'$j = j\times B$'
		mhd_vars['omega_A'] = r'$omega_A$'
		mhd_vars['Fmag_par'] = r'$F_\parallel$'
		mhd_vars['Fmag_phi'] = r'$F_\phi$'
		mhd_vars['gradP_par'] = r'$(\nabla p)_\parallel$'
		mhd_vars['Be_mag'] = r'$Be_{mag}$'

		# label for my_var, user_out_var, floor maps, & dt maps
		misc_vars['tau'] = r'$\tau$'
		misc_vars['xi'] = r'$\xi$'
		misc_vars['dlnA'] = r'$-\hat{n}\cdot\nabla(\rho v_{poloidal})/\rho v_{poloidal}$'
		misc_vars['div_v'] = r'$\nabla \cdot v$'
		misc_vars['my_var'] = my_var_label(0)
		misc_vars['user_out_var'] = user_out_var_label(0)
		misc_vars['d_floor'] = r'$\rho_{floor}$' + ' triggers'
		misc_vars['p_floor'] = r'$p_{floor}$' + ' triggers'
		misc_vars['dt_cs_r'] = r'$dt_{cs,r}$'
		misc_vars['dt_cs_t'] = r'$dt_{cs,\theta}$'
		misc_vars['dt_v_r'] = r'$dt_{v,r}$'
		misc_vars['dt_v_t'] = r'$dt_{v,\theta}$'
		misc_vars['dt_vA_r'] = r'$dt_{v_A,r}$'
		misc_vars['dt_vA_t'] = r'$dt_{v_A,\theta}$'


	# add all the labels to LabelDic
	LabelDic['hydro_keys'] = hydro_vars.keys()
	LabelDic['vorticity_keys'] = ['w_r','w_t','w_p','wp','w']
	LabelDic['misc_keys'] = misc_vars.keys()
	LabelDic['mhd_keys'] = mhd_vars.keys()
	for key in LabelDic['hydro_keys']:
		LabelDic[key] = hydro_vars[key]
	for key in LabelDic['misc_keys']:
		LabelDic[key] = misc_vars[key]
	for key in LabelDic['mhd_keys']:
		LabelDic[key] = mhd_vars[key]
	
	return LabelDic

def generate_grid_data1D(dump_dirs,i_dump=0,i_min=0,i_max=None,ext='tab'):
	""" Loads multiple simulation runs into a single dictionary
	
	Inputs:
	dump_dirs - paths to all sims (the output of dump_dir_dic())
	i_dump - integer specifying dump #
	"""
	print("[generate_grid_data]:")

	all_keys = dump_dirs.keys()
	keys_to_use = []
	for key in all_keys:
		keys_to_use.append(key[0])
	keys_to_use = np.unique(keys_to_use)
	
	Nsims = len(keys_to_use)
	
	GridDics = {}
	
	for n in range(Nsims):
		
		idx = np.str(n+1)
		
		# open the vtk file using athena_read
		if i_dump >= 0:
			dump_file = dump_dirs[idx] + str(i_dump).zfill(5) + '.' + ext
		else:
			dump_file = dump_dirs[idx] + ext
		if ext == 'tab':
			athena_dic = my_athena_reader.read_tab(dump_file)
		elif ext == 'vtk':
			athena_dic = athena_read.vtk(dump_file)
		else:
			athena_dic = athena_read.athdf(dump_file)
		
		# store output of each sim to a dictionary
# 		x_c = athena_dic[0][0][:,0][i_min:i_max]
		x_c = athena_dic['x1v'][i_min:i_max]
		grid_dic = {'x_c':x_c}  # position data

		grid_dic['i_min'],grid_dic['i_max'] = i_min,i_max
		
		# store grid_dic in GridDics
		print("\tDone.")
		print("\tGridDic[{}] has the following keys:\n{}".format(n+1,grid_dic.keys()))
		GridDics[idx] = grid_dic

	
	return GridDics


def generate_grid_data(dump_dirs,i_dump=0,zoom_region=None,R_zoom=None,iskip=1,jskip=1,kskip=1,ext='vtk'):
	""" Loads multiple simulation runs into a single dictionary
	
	Inputs:
	dump_dirs - paths to all sims (the output of dump_dir_dic())
	i_dump - integer specifying dump #

	Example usage:
	vars = ['d', 'v_r']
	SimDic = load_vars(i_dump,vars,dump_dirs,grid_dics)
	Sim['1']['d'] is an array of density for the 1st sim
	Sim['2']['v_r'] is an array of r-velocity for the 2nd sim
	etc...
	"""
	print("[generate_grid_data]:")

	all_keys = dump_dirs.keys()
	keys_to_use = []
	for key in all_keys:
		keys_to_use.append(key[0])
	keys_to_use = np.unique(keys_to_use)
	
	Nsims = len(keys_to_use)
	
	GridDics = {}
	
	for n in range(Nsims):
		
		idx = np.str(n+1)
		
		# open the vtk file using athena_read
		dump_file = dump_dirs[idx] + str(i_dump).zfill(5) + '.' + ext
		if ext == 'vtk':
			athena_dic = athena_read.vtk(dump_file)
		else:
			athena_dic = athena_read.athdf(dump_file)


		# store output of each sim to a dictionary
		grid_dic = {}
		
		# Store zoom_region as indices
		if zoom_region is None:
			grid_dic['i_min'] = 0
			grid_dic['i_max'] = None
			grid_dic['j_min'] = 0
			grid_dic['j_max'] = None
			grid_dic['k_min'] = 0
			grid_dic['k_max'] = None
		else:
			grid_dic['i_min'] = zoom_region['i_min']
			grid_dic['i_max'] = zoom_region['i_max']
			grid_dic['j_min'] = zoom_region['j_min']
			grid_dic['j_max'] = zoom_region['j_max']

			if len(athena_dic['x3f']) > 2:
				grid_dic['k_min'] = zoom_region['k_min']
				grid_dic['k_max'] = zoom_region['k_max']
			else:
				grid_dic['k_min'] = 0
				grid_dic['k_max'] = None

		i_min = grid_dic['i_min']
		i_max = grid_dic['i_max']
		j_min = grid_dic['j_min']
		j_max = grid_dic['j_max']
		k_min = grid_dic['k_min']
		k_max = grid_dic['k_max']

		print(i_min,i_max)
		print(j_min,j_max)
		

		# calculate cell centered grid'
		is3D = False
		if ext == 'athdf':
			if i_max == None:
				x = athena_dic['x1f'][i_min:]
			else:
				x = athena_dic['x1f'][i_min:i_max+1]
			if j_max == None:
				y = athena_dic['x2f'][j_min:]
			else:
				y = athena_dic['x2f'][j_min:j_max+1]

			if len(athena_dic['x3f']) > 2:
				is3D = True
				if k_max == None:
					z = athena_dic['x3f'][k_min:]
				else:
					z = athena_dic['x3f'][k_min:k_max+1]
		else:
			if i_max == None:
				x = athena_dic[0][i_min:]
			else:
				x = athena_dic[0][i_min:i_max+1]
			if j_max == None:
				y = athena_dic[1][j_min:]
			else:
				y = athena_dic[1][j_min:j_max+1]

			if athena_dic[2].shape[0] > 1:
				is3D = True
				if k_max == None:
					z = athena_dic[2][k_min:]
				else:
					z = athena_dic[2][k_min:k_max+1]

		x_c = x[:-1] + 0.5 * np.ediff1d(x)
		y_c = y[:-1] + 0.5 * np.ediff1d(y)
		grid_dic['x'] =x
		grid_dic['y'] = y
		grid_dic['x_c'] = x_c
		grid_dic['y_c'] = y_c
		if is3D:
			z_c = z[:-1] + 0.5 * np.ediff1d(z)
			grid_dic['z'] = z
			grid_dic['z_c'] = z_c
	
		# Determine zoom-in index, i_max, corresponding to x_zoom
		Nx = int(len(x_c)) 
		if R_zoom == None or R_zoom < x_c[0] or R_zoom > x_c[Nx-1]:
			print("\t [generate_grid_data]: ")
			print("\t R_zoom is outside r-range = [{},{}]".format(x_c[0],x_c[Nx-1]))
			print("\t Setting x_zoom to max(x)")
			i_max2 = Nx
			i_zoom = Nx-1
		else:
			R_diag = 1.1*np.sqrt(2.)*R_zoom
			dxs_cut = np.abs(r_c - R_diag)
			dxs_zoom = np.abs(r_c - R_zoom)
			i_max2 = min(Nr-1,dRs_cut.argmin())
			i_zoom = min(Nr-1,dRs_zoom.argmin())
		grid_dic['i_zoom'] = i_zoom
		# reset i_max
		# grid_dic['i_max'] = i_max
		
		print("\tGenerating various grids for simulation {}...\n".format(n+1))
	
		# Generate Cartesian cell edge grid
		# grid_dic['r_grid'], grid_dic['theta_grid'] = np.meshgrid(r[0:i_max+1],theta)
		# grid_dic['x'] = grid_dic['r_grid']*np.sin(grid_dic['theta_grid'])
		# grid_dic['y'] = grid_dic['r_grid']*np.cos(grid_dic['theta_grid'])

		# # Generate Cartesian cell-centered grid
		# sin_thetas, cos_thetas = np.sin(theta_c_grid),np.cos(theta_c_grid)
		# grid_dic['r_c_grid'], grid_dic['theta_c_grid'] = r_c_grid, theta_c_grid
		# grid_dic['x_c'], grid_dic['y_c']  = r_c_grid*sin_thetas, r_c_grid*cos_thetas
		x_c_grid, y_c_grid = np.meshgrid(x_c,y_c)
		grid_dic['x_c_grid'], grid_dic['y_c_grid'] = x_c_grid, y_c_grid
		grid_dic['x_c_avg'] = 0.5*(x_c_grid[1:,1:] + x_c_grid[:-1,:-1])
		grid_dic['y_c_avg'] = 0.5*(y_c_grid[1:,1:] + y_c_grid[:-1,:-1])
		
		# # Save the grid of sin_thetas/cos_thetas 
		# grid_dic['sin_thetas'], grid_dic['cos_thetas'] = sin_thetas, cos_thetas

		# Calculate index grids
		i_grid = np.zeros_like(x_c_grid)
		for j in range(len(y_c[j_min:j_max])):
			i_grid[j,:] = j
			for i in range(len(x_c[i_min:i_max])):
				i_grid[j,i] = i
			
		j_grid = np.zeros_like(y_c_grid)
		for i in range(len(x_c[i_min:i_max])):
			j_grid[:,i] = i
			for j in range(len(y_c[j_min:j_max])):
				j_grid[j,i] = j

		grid_dic['i_grid'], grid_dic['j_grid'] = i_grid, j_grid

		# Calculate sparse grids for displaying vector fields - 0 to j_min indices
		x_sparse = x_c_grid[j_min::jskip,i_min::iskip]
		y_sparse = y_c_grid[j_min::jskip,i_min::iskip]
		i_sparse = grid_dic['i_grid'][j_min::jskip,i_min::iskip]
		j_sparse = grid_dic['j_grid'][j_min::jskip,i_min::iskip]
		grid_dic['x_sparse0'],grid_dic['y_sparse0'] = x_sparse, y_sparse 
		grid_dic['i_sparse0'],grid_dic['j_sparse0'] = i_sparse, j_sparse 
	
		# Calculate sparse grids for displaying vector fields - j_min to j_max indices
		x_sparse = x_c_grid[j_min::jskip,i_min::iskip]
		y_sparse = y_c_grid[j_min::jskip,i_min::iskip]
		i_sparse = i_grid[j_min::jskip,i_min::iskip]
		j_sparse = j_grid[j_min::jskip,i_min::iskip]
		grid_dic['x_sparse'],grid_dic['y_sparse'] = x_sparse, y_sparse 
		grid_dic['i_sparse'],grid_dic['j_sparse'] = i_sparse, j_sparse 


		# Grid of cartesian cell sizes used for streamline/B-field line calc
		# lengths=np.ones([len(theta_c),len(r_c)])
		# #Loop over the r and theta arrays
		# for j in range (len(theta_c)):
		# 	for i in range(len(r_c)):
		# 		if j>0:	
		# 			dl1=np.abs(r_c[i]*np.cos(theta_c[j])-r_c[i]*np.cos(theta_c[j-1]))
		# 		else:
		# 			dl1=np.abs(r_c[i]*np.cos(theta_c[j])-r_c[i]*np.cos(theta_c[j+1]))
		# 		if i>0:
		# 			dl2=np.abs(r_c[i]-r_c[i-1])
		# 		else:
		# 			dl2=np.abs(r_c[i]-r_c[i+1])

		# 		lengths[j][i]=np.min([dl1,dl2])
		# grid_dic['lengths'] = lengths
		
		# store grid_dic in GridDics
		print("\tDone.")
		print("\tGridDic[{}] has the following keys:\n{}".format(n+1,grid_dic.keys()))
		GridDics[idx] = grid_dic
	
	return GridDics


def calc_grad_r(qty,r_c_grid):
	"""
	Take the gradient with respect to r.
	r_c_grid is the actual grid used in Athena,
	so we allow for a non-uniform mesh here.
	"""

	grad_qty = np.zeros_like(qty)

	#for j in range(qty.shape[0]):
	grad_qty[:,:-1] = (qty[:,1:]-qty[:,:-1])/(r_c_grid[:,1:]-r_c_grid[:,:-1])

	# now copy the gradients at r_out to the 
	# last r_c zone in qty so the shapes match 
	grad_qty[:,-1] = grad_qty[:,-2]

	return grad_qty

def calc_grad_t(qty,theta_c_grid):
	"""
	Take the gradient with respect to theta.
	theta_c_grid is the actual grid used in Athena,
	so we allow for a non-uniform mesh here.
	"""
	
	grad_qty = np.zeros_like(qty)

	grad_qty[:-1,:] = (qty[1:,:]-qty[:-1,:])/(theta_c_grid[1:,:]-theta_c_grid[:-1,:])

	# now copy the gradients at r_out to the 
	# last r_c zone in qty so the shapes match 
	grad_qty[-1,:] = grad_qty[-2,:]

	return grad_qty

def add_gradients(var_keys,sim_dics,grid_dics):
	
	Nsims = len(grid_dics.keys()) 
	
	for n in range(Nsims):
		
		idx = np.srt(n+1)
		
		# load grid for this sim
		grid_dic = grid_dics[idx]

		# we will add the gradients to sim_dic
		sim_dic = sim_dics[idx]

		for key in var_keys:
			grad_r_key = 'grad_' + key + '_r'
			grad_t_key = 'grad_' + key + '_t'
			sim_dic[grad_r_key] = calc_grad_r(sim_dic[key],grid_dic['r_c_grid'])
			sim_dic[grad_t_key] = calc_grad_t(sim_dic[key],grid_dic['theta_c_grid'])

def calc_div(q_r,q_t,grid_dic):
	
	r_c_grid,theta_c_grid = grid_dic['r_c_grid'],grid_dic['theta_c_grid']
	cos_theta_avg = 0.5*(grid_dic['cos_thetas'][1:,1:] + grid_dic['cos_thetas'][:-1,:-1])
	sin_theta_avg = 0.5*(grid_dic['sin_thetas'][1:,1:] + grid_dic['sin_thetas'][:-1,:-1])
	
	dq_dr = (q_r[:,1:]-q_r[:,:-1])/(r_c_grid[:,1:]-r_c_grid[:,:-1])
	dq_dtheta = (q_t[1:,:]-q_t[:-1,:])/(theta_c_grid[1:,:]-theta_c_grid[:-1,:])
	
	r_c_avg = 0.5*(r_c_grid[1:,1:] + r_c_grid[:-1,:-1])
	theta_c_avg = 0.5*(theta_c_grid[1:,1:] + theta_c_grid[:-1,:-1])
	q_r_avg = 0.5*(q_r[1:,1:] + q_r[:-1,:-1])
	q_t_avg = 0.5*(q_t[1:,1:] + q_t[:-1,:-1])
	div_q = dq_dr[0:-1,:] + 2.*q_r_avg/r_c_avg + \
			dq_dtheta[:,0:-1]/r_c_avg + \
			(q_t_avg/r_c_avg)*(cos_theta_avg/sin_theta_avg)
	
	return div_q

def calc_optical_depth(grid_dic,dens):
	print("\n\n [calc_optical_depth]: WTF! \n\n")

	# r_grid, theta_grid = np.meshgrid(r,theta_c)
	r_grid,theta_grid = grid_dic['r_grid'],grid_dic['theta_c_grid']
	tau_grid = np.zeros_like(dens)
	
	r,r_c,theta_c = grid_dic['r'],grid_dic['r_c'],grid_dic['theta_c']
	drs = r[1:]-r[:-1]
	for j in range(len(theta_c)):
		for i in range(1,len(r_c)+1):
			tau_grid[j,i-1] = np.sum(dens[j,0:i]*drs[0:i])
	tau_grid *= units.tau

	return tau_grid

def PBCshift(Npts,i):
	if i < 0:
		return (Npts + i)
	elif i >= Npts:
		return (i - Npts)		
	else:
		return i

def pbc_gradient(var,x):
	# calculate gradient assuming periodic BCs
	dx_inv = 1./(x[1]-x[0])
	N = len(var)
	var_prime = np.zeros(N)
	for i in range(N):
		vL = 0.5*(var[i]+var[PBCshift(N,i-1)])
		vR = 0.5*(var[i]+var[PBCshift(N,i+1)])
		var_prime[i] = dx_inv*(vR-vL)
	return var_prime

def load_vars1D(i_dump,vars,dump_dirs,grid_dics,gamma=5./3.,MHD_check=0,ext='tab'):
	""" Loads multiple simulation runs into a single dictionary
	
	Inputs:
	i_dump - integer specifying dump # 
	vars - array of variables to collect
	dump_dirs - paths to all sims (the output of dump_dir_dic())
	grid_dics - grid data for all sims (the output of generate_grid_data())
	i_maxs - if set to True, the size of all output arrays will be limited
			 to i_max, which is contained in grid_dics
	gamma - adiabatic index

	Returns: A dictionary containing data of the variables listed in vars for
	each simulation specified in dump_dirs for dump i_dump

	Example usage:
	vars = ['d', 'v_r']
	SimDic = load_vars(i_dump,vars,dump_dirs,grid_dics)
	Sim['1']['d'] is an array of density for the 1st sim
	Sim['2']['v_r'] is an array of r-velocity for the 2nd sim
	etc...
	"""
	
	Nsims = len(grid_dics.keys()) 
	
	# load all keys
	label_dic = load_labels1D()

	print("\n\n [load_vars]: vars = {} \n\n".format(vars))
	
	global VarDics
	VarDics = {}
	for n in range(Nsims):
		
		idx = np.str(n+1)
		
		# load grid for this sim
		grid_dic = grid_dics[idx]
		
		# check if i_min
		if grid_dic['i_min'] == 0:
			i_min = 0
		else:
			i_min = grid_dic['i_min']

		# check if i_max
		if grid_dic['i_max'] is None:
			i_max = None
		else:
			i_max = grid_dic['i_max']
		
		# open the vtk file using athena_read
		if i_dump >= 0:
			dump_file = dump_dirs[idx] + str(i_dump).zfill(5) + '.' + ext
		else:
			dump_file = dump_dirs[idx] + ext
		if ext == 'tab':
			athena_dic = my_athena_reader.read_tab(dump_file)
		elif ext == 'vtk':
			athena_dic = athena_read.vtk(dump_file)
		else:
			athena_dic = athena_read.athdf(dump_file)
		
		# store output of each sim to 1 dictionary
		VarDics[idx] = {}
		
		# these are set to 0 if vars in the same group have been 
		# calculated to prevent duplicate computations
		global group1, group2, group3, group4, group5
		group1, group2, group3, group4, group5 = 1, 1, 1, 1, 1
		
		# unpack primitive variables
		global x,d,p,v,T
		x = athena_dic['x1v'][i_min:i_max]
		d = athena_dic['rho'][i_min:i_max]*units.d
		p = athena_dic['press'][i_min:i_max]*units.p
		v_x = athena_dic['vel1'][i_min:i_max]*units.v
		v_y = athena_dic['vel2'][i_min:i_max]*units.v
		v = np.sqrt(v_x**2 + v_y**2)
		T = units.T*p/d
		xi = units.xi/d/x**2
		cs = units.cs*np.sqrt(T)
		Be = 0.5*v**2 + cs**2/(gamma - 1.) 
		
		# user out variables flagged with this key
		uov_key = idx + 'uov'

		global MHDrun 
		try:
			B_r = athena_dic[3]['Bcc'][0,:,0:i_max,0]
			B_t = athena_dic[3]['Bcc'][0,:,0:i_max,1]
			B_p = athena_dic[3]['Bcc'][0,:,0:i_max,2]
			MHDrun = True	
		except:
			MHDrun = False	

		for var in vars:
			try:
				if var == 'd':
					qty = d
				elif var == 'p':
					qty = p
				elif var == 'v': #total velocity
					qty = v	
				elif var == 'v_x': #total velocity
					qty = v_x	
				elif var == 'v_y': #total velocity
					qty = v_y	
				elif var == 'T': #temperature
					qty = T
				elif var == 'mfd': # mass flux density
					qty = d*v
				elif var == 'cs': # adiabatic sound speed
					qty = cs
				elif var == 'c_iso': # adiabatic sound speed
					qty = np.sqrt(p/d)
				elif var == 's': # entropy
					qty = np.log(p/d**gamma)
				elif var == 'M': # poloidal Mach number
					qty = v/cs	
				elif var == 'Be': # Bernoulli function 
					qty = Be #0.5*v**2 + cs**2/(gamma - 1.)
				elif var == 'v_sqd': #total velocity
					qty = v*v
				elif var == 'KE': # Bernoulli function 
					qty = 0.5*d*v**2 
				elif var == 'TE': # Bernoulli function 
					qty = 0.5*d*v**2 + p/(gamma - 1.) 
				 

				# other derived hydro quantities				
				elif  var in label_dic['hc_keys']: 
					
					if group1:
						group1 = 0

						xi_eq = 190.
						hcModel = blondin.HCcurves(gamma,nT=1e13)
						eqvals = hcModel.getEquilibriumValues(xi_eq)
						# athena = hcModel.getCodeUnits(eqvals,1.)
						# L_T,L_rho = hcModel.getL_TandL_rho(xi_eq,eqvals)
					xi = xi_eq/d
					gfac = 1./(gamma*(gamma-1.))
					hcModel.setCodeUnits(eqvals)
					netL = gfac*d*hcModel.netLosses(xi,T)
					netC = gfac*d*(hcModel.Lb(T) + hcModel.Ll(xi,T))
					netH = gfac*d*(hcModel.Gc(xi,T) + hcModel.Gx(xi,T))
					# netL = hcModel.netLosses(xi,T*eqvals['T0'])*eqvals['norm']
					hcModel.resetCgsUnits()
					if var == 'netL': 
						qty = netL
					if var == 'rhoL': 
						qty = d*netL
					elif var == 'netC': 
						qty = netC
					elif var == 'netH': 
						qty = netH
					elif var == 'netLoverT':
						qty = netL/T
					elif var == 'netLoverv': 
						qty = netL/v
					elif var == 'gradBe':
						Be = 0.5*v**2 + cs**2/(gamma - 1.) 
						qty = pbc_gradient(Be,x)
					elif var == 'div_q':
						qty = -units.kappa*pbc_gradient(pbc_gradient(T,x),x)
					elif var == 'grad_s':
						s = np.log(p/d**gamma)
						qty = pbc_gradient(s,x)
					
				# MHD quantities	
				elif var in label_dic['mhd_keys']:
		
					if group2:
						try:
							B = np.sqrt(B_r**2 + B_t**2 + B_p**2)
							Bp = np.sqrt(B_r**2 + B_t**2) # poloidal field
							v_A = B/np.sqrt(d)
							v_Ap = Bp/np.sqrt(d)
							# compute fast and slow wave speeds, e.g. Shu Vol 2 page 307
							cos_phi = np.cos(Bp/B)
							vA2 = v_A**2
							cs2 = cs**2
							discrim = (vA2 + cs2)**2 - 4.*vA2*cs2*cos_phi**2
							sqrt_check = np.min(discrim)
							if sqrt_check < 0.:
								print("[SimDics.load_vars(): undefined sqrt! -> min(discrim) = {}".format(sqrt_check))
								print("[SimDics.load_vars(): taking an abs() meaning v_fast/v_slow are not defined!!")
								sqrt_term = np.sqrt(np.abs(discrim))
							else:
								sqrt_term = np.sqrt(discrim)
							v_fast = np.sqrt(0.5*(vA2 + cs2 + sqrt_term))
							v_slow = np.sqrt(0.5*(vA2 + cs2 - sqrt_term))
							group2 = 0
						except:
							if MHD_check: 
								continue
							else:
								raise ValueError(var)
					if var == 'B_r':
						qty = B_r
					elif var == 'B_t':
						qty = B_t
					elif var == 'B_p':
						qty = B_p
					elif var == 'v_A': # total Alven speed
						qty = v_A
					elif var == 'v_Ap': # poloidal Alven speed
						qty = v_Ap
					elif var == 'B': # magntude of B-field
						qty = B
					elif var == 'Bp':  # poloidal B-field
						qty = Bp
					elif var == 'E_r':
						qty = (v_p*B_t - v_t*B_p)/units.c
					elif var == 'E_t':
						qty = (v_r*B_p - v_p*B_r)/units.c
					elif var == 'E_p':
						qty = (v_t*B_r - v_r*B_t)/units.c
					elif var == 'Ep':
						qty = np.sqrt(E_r**2 + E_t**2)
					elif var == 'E':
						qty = np.sqrt(E_r**2 + E_t**2 + E_p**2)
					elif var == 'E_p_over_Ep':
						E_r = (v_p*B_t - v_t*B_p)/units.c
						E_t = (v_r*B_p - v_p*B_r)/units.c
						E_p = (v_t*B_r - v_r*B_t)/units.c
						qty = E_p/np.sqrt(E_r**2 + E_t**2)
					elif var == 'beta': # plasma beta
						qty = 2.*p/Bp**2
					elif var == 'B2': # plasma beta
						qty = B**2
					elif var == 'Bphi_over_Bp': 
						qty = B_p/Bp
					elif var == 'v_fast': 
						qty = v_fast
					elif var == 'v_slow': 
						qty = v_slow
					elif var == 'A': # See Contopoulos 1996 eq 14
						v_x = v_r*grid_dic['sin_thetas'] + v_t*grid_dic['cos_thetas']
						v_z = v_r*grid_dic['cos_thetas'] - v_t*grid_dic['sin_thetas']
						B_x = B_r*grid_dic['sin_thetas'] + B_t*grid_dic['cos_thetas']
						B_z = B_r*grid_dic['cos_thetas'] - B_t*grid_dic['sin_thetas']
						x_c = grid_dic['x_c']
						qty = x_c*(v_z*B_x - v_x*B_z)
					elif var == 'omega_A': # fastness parameter, e.g. Tzeferacos+13 eq 23
						vp = np.sqrt(v_r**2 + v_t**2)
						kappa = d*vp/Bp
						Omega = v_p/grid_dic['x_c'] - kappa*B_p/d/grid_dic['x_c']
# 						M_Ap = vp/v_Ap
# 						filter = np.where(np.abs(M_Ap - 1.)) < 1e-3
						qty = Omega*grid_dic['x_c']/v_Ap
					elif var == 'Fmag_par':
						grad_Bphi_r = calc_grad_r(grid_dic['x_c']*B_p,grid_dic['r_c_grid'])
						grad_Bphi_t = calc_grad_t(grid_dic['x_c']*B_p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_Bphi_par = n_r*grad_Bphi_r + n_t*grad_Bphi_t
						qty = -(B_p/grid_dic['x_c'])*grad_Bphi_par
					elif var == 'Fmag_phi':
						grad_Bphi_r = calc_grad_r(grid_dic['x_c']*B_p,grid_dic['r_c_grid'])
						grad_Bphi_t = calc_grad_t(grid_dic['x_c']*B_p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_Bphi_par = n_r*grad_Bphi_r + n_t*grad_Bphi_t
						qty = (Bp/grid_dic['x_c'])*grad_Bphi_par
					elif var == 'gradP_par':
						grad_p_r = calc_grad_r(p,grid_dic['r_c_grid'])
						grad_p_t = calc_grad_t(p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_p_par = n_r*grad_p_r + n_t*grad_p_t
						qty = -grad_p_par
					elif var == 'Be_mag':
						vp = np.sqrt(v_r**2 + v_t**2)
						kappa = d*vp/Bp
						Omega = v_p/grid_dic['x_c']
						Omega_mag = -kappa*B_p/d/grid_dic['x_c']
						Omega_B = Omega + Omega_mag
						magnetic = -(Omega_B/kappa)*grid_dic['x_c']*B_p
						kinetic_poloidal = 0.5*vp**2 
						kinetic_rotational = 0.5*v_p**2 
						enthalpy = cs**2/(gamma-1.)
						gravity = -units.GM/grid_dic['r_c_grid']
						qty = kinetic_poloidal + kinetic_rotational + enthalpy + gravity + magnetic
						
						
					if uov_key in dump_dirs.keys():
						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
						print("[load_vars]: uov_file = {}".format(uov_file))
						uov_dic = athena_read.vtk(uov_file)
						
						# unpack current density, j
						j_r = uov_dic[3]['user_out_var3'][0,:,0:i_max]
						j_t = uov_dic[3]['user_out_var4'][0,:,0:i_max]
						j_p = uov_dic[3]['user_out_var5'][0,:,0:i_max]

						# get rid of any nan's along north pole
						if np.abs(j_t[0,0]) == np.inf:
							j_t[0,:] = j_t[-1,0] # use south pole values on north pole
						
						jp = np.sqrt(j_r**2 + j_t**2)
						j = np.sqrt(j_r**2 + j_t**2 + j_p**2)
						

						
						# now compute the Lorentz force, jxB
						jcrossB_r = j_t*B_p - j_p*B_t
						jcrossB_t = j_r*B_p - j_p*B_r
						jcrossB_p = j_r*B_t - j_t*B_r
						jcrossBp = np.sqrt(jcrossB_r**2 + jcrossB_t**2)
						jcrossB = np.sqrt(jcrossB_r**2 + jcrossB_t**2 + jcrossB_p**2)
						
						if var == 'j_r':
							qty = j_r
						elif var == 'j_t':
							qty = j_t
						elif var == 'j_p':
							qty = j_p
						elif var == 'jp':
							qty = jp
						elif var == 'j':
							qty = j
						elif var == 'jcrossB_r':
							qty = jcrossB_r
						elif var == 'jcrossB_t':
							qty = jcrossB_t
						elif var == 'jcrossB_p':
							qty = jcrossB_p
						elif var == 'jcrossBp':
							qty = jcrossBp
						elif var == 'jcrossB':
							qty = jcrossB
					else:
						print("Post-processing current density calculation not yet implemented!")

				# quantities derived from optical depth 
				elif var in ['tau','xi']:
					if group3:
						tau = calc_optical_depth(grid_dic,d)
						group3 = 0
					if var == 'tau':
						qty = tau
					elif var == 'xi':  # photoionization parameter
						qty = units.xi*np.exp(-tau)/(d*grid_dic['r_c_grid']**2)	
						
				# flow Area diagnostics
				elif var in ['dlnA','div_v']:
					if var == 'dlnA':
						vp = np.sqrt(v_r**2 + v_t**2)
						mfd = d*vp
						grad_mfd_r = calc_grad_r(mfd,grid_dic['r_c_grid'])
						grad_mfd_t = calc_grad_t(mfd,grid_dic['theta_c_grid'])
						n_r,n_t = v_r/vp,v_t/vp
						dlnA = -(n_r*grad_mfd_r + n_t*grad_mfd_t)/mfd
						qty = dlnA
					elif var == 'div_v':
# 						p_avg = 0.5*(p[1:,1:] + p[:-1,:-1])
						qty = calc_div(v_r,v_t,grid_dic)
					group5 = 0
						

				# numerics: cell crossing times and floors - e.g., dt_cs_t = r(dtheta)/|cs|
				elif var in ['dt_cs_r','dt_cs_t','dt_v_r','dt_v_t','dt_vA_r','dt_vA_t','d_floor','p_floor']:
					if group4: # need all speeds
						cs = units.cs*np.sqrt(gamma*p/d)
						v = np.sqrt(v_r**2 + v_t**2 + v_p**2)

						dtheta_grid = grid_dic['theta_grid'][1:,:] - grid_dic['theta_grid'][0:-1,:] 
						rdtheta_grid = grid_dic['r_c_grid']*dtheta_grid[:,0:-1]
						dr_grid = grid_dic['r_grid'][:,1:] - grid_dic['r_grid'][:,0:-1]
						dr_grid = dr_grid[0:-1,:]

						group4 = 0

					if var == 'dt_cs_r': # sound crossing time in radial dir
						qty = dr_grid/cs
					elif var == 'dt_cs_t': # sound crossing time in polar dir
						qty = rdtheta_grid/cs
					if var == 'dt_v_r': # flow crossing time in radial dir
						qty = dr_grid/v
					elif var == 'dt_v_t': # flow crossing time in polar dir
						qty = rdtheta_grid/v
					elif var == 'dt_vA_r': # Alven crossing time in radial dir
						qty = dr_grid/v_A
					elif var == 'dt_vA_t': # Alven crossing time in polar dir
						qty = rdtheta_grid/v_A
					elif var == 'd_floor': # map of zones hitting dens floor
						qty = np.copy(d)
						qty[max_filter] = 1e20 # np.nan
					elif var == 'p_floor': # map of zones hitting pres floor
						qty = np.copy(p)
						qty[min_filter] = 1e20 #np.nan
					
				
				# scratch variable		
				elif var == 'my_var':
					if MHDrun:
						qty = my_var(grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p)
					else:
						qty = my_var(grid_dic,d,p,v_r,v_t,v_p)
				
				# Athena++ user output scalar variable
# 				elif var.startswith('user_out'): #and not(uov_vector):
# 					# open the vtk file using athena_read
# 					uov_key = idx + 'uov'
# 					if uov_key in dump_dirs.keys():
# 						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
# 						uov_dic = athena_read.vtk(uov_file)
# 						qty = uov_dic[3][var][0,:,0:i_max]
# 						qty = uov_scalar(qty,int(var[-1]))

				# Athena++ user output vector variable
				elif var.startswith('user_out'):
					print("\n\n USER_OUT_VAR!! \n\n")
					
					index = int(var[-1])
					var_base_name = var[:-1]
					uovs = []

					# open the vtk file using athena_read
					uov_key = idx + 'uov'
					if uov_key in dump_dirs.keys():
						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
						print("[load_vars]: uov_file = {}".format(uov_file))
						uov_dic = athena_read.vtk(uov_file)
						for n in range(index+1):
							uov_var = var_base_name + np.str(n)
							print(uov_var)
							uovs.append(uov_dic[3][uov_var][0,:,0:i_max])
						if MHDrun:
							uovs = uov_vector(uovs,grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p)
						else:
							uovs = uov_vector(uovs,grid_dic,d,p,v_r,v_t,v_p)
						
						# if vector uov, add each component to VarDics
						if len(uovs) > 1:
							for n in range(index+1):
								uov_var = var_base_name + np.str(n)
								VarDics[idx][uov_var] = uovs[n]
							qty = VarDics[idx][var] 
						else:
							qty = uovs[0]

				else: # variable not recognized - throw exception 
					raise ValueError(var)	
			
				VarDics[idx][var] = qty

			except ValueError as err:
				label_dic = load_labels()
# 				print "[loadvars()]: label_dic['mhd_keys'] = {}".format(label_dic['mhd_keys'])
				if var in label_dic['mhd_keys']:
					print("\n******\nERROR!\n******")
					print("MHD variable {} was recognized but not found for simulation:\n{}".format(err.args[0],dump_file))
					print("This must not be an MHD simulation.")
					print("Choose from these hydro quantities:")
					print(label_dic['hydro_keys'])
				else:
					print("\n******\nERROR!\n******")
					print("Variable {} not recognized. Choose from:".format(err.args[0]))
					print("hydro vars:\n{}".format(label_dic['hydro_keys']))
					print("mhd vars:\n{}".format(label_dic['mhd_keys']))
					print("misc vars:\n{}".format(label_dic['misc_keys']))
				sys.exit()

	
	return VarDics

def load_vars(i_dump,vars,dump_dirs,grid_dics,gamma=5./3.,MHD_check=0,ext='vtk',k=0):
	""" Loads multiple simulation runs into a single dictionary
	
	Inputs:
	i_dump - integer specifying dump # 
	vars - array of variables to collect
	dump_dirs - paths to all sims (the output of dump_dir_dic())
	grid_dics - grid data for all sims (the output of generate_grid_data())
	i_maxs - if set to True, the size of all output arrays will be limited
			 to i_max, which is contained in grid_dics
	gamma - adiabatic index

	Returns: A dictionary containing data of the variables listed in vars for
	each simulation specified in dump_dirs for dump i_dump

	Example usage:
	vars = ['d', 'v_r']
	SimDic = load_vars(i_dump,vars,dump_dirs,grid_dics)
	Sim['1']['d'] is an array of density for the 1st sim
	Sim['2']['v_r'] is an array of r-velocity for the 2nd sim
	etc...
	"""
	
	Nsims = len(grid_dics.keys()) 
	
	# load all keys
	label_dic = load_labels()

	print("\n\n [load_vars]: vars = {} \n\n".format(vars))
	
	global VarDics
	VarDics = {}
	for n in range(Nsims):
		
		idx = np.str(n+1)
		
		# load grid for this sim
		grid_dic = grid_dics[idx]
		
		# unpack region indices
		i_min = grid_dic['i_min']
		i_max = grid_dic['i_max']
		j_min = grid_dic['j_min']
		j_max = grid_dic['j_max']
		k_min = grid_dic['k_min']
		k_max = grid_dic['k_max']
		
		# open the vtk file using athena_read
		dump_file = dump_dirs[idx] + str(i_dump).zfill(5) + '.' + ext
		if ext == 'vtk':
			athena_dic = athena_read.vtk(dump_file)
		else:
			athena_dic = athena_read.athdf(dump_file)

		
		# store output of each sim to 1 dictionary
		VarDics[idx] = {}
		
		# these are set to 0 if vars in the same group have been 
		# calculated to prevent duplicate computations
		global group1, group2, group3, group4, group5
		group1, group2, group3, group4, group5 = 1, 1, 1, 1, 1
		
		# unpack primitive variables
		global d,p,v_r,v_t,v_p
		if ext == 'athdf':
			d = athena_dic['rho'][k,j_min:j_max,i_min:i_max]
			p = athena_dic['press'][k,j_min:j_max,i_min:i_max]
			v_x = athena_dic['vel1'][k,j_min:j_max,i_min:i_max]
			v_y = athena_dic['vel2'][k,j_min:j_max,i_min:i_max]
			v_z = athena_dic['vel3'][k,j_min:j_max,i_min:i_max]	
		else:
			d = athena_dic[3]['rho'][k,j_min:j_max,i_min:i_max]
			p = athena_dic[3]['press'][k,j_min:j_max,i_min:i_max]
			v_x = athena_dic[3]['vel'][k,j_min:j_max,i_min:i_max,0]
			v_y = athena_dic[3]['vel'][k,j_min:j_max,i_min:i_max,1]
			v_z = athena_dic[3]['vel'][k,j_min:j_max,i_min:i_max,2]
		cs = units.cs*np.sqrt(gamma*p/d)
		
		# user out variables flagged with this key
		uov_key = idx + 'uov'

		global MHDrun 
		try:
			B_r = athena_dic['Bcc1'][k,j_min:j_max,i_min:i_max]
			B_t = athena_dic['Bcc2'][k,j_min:j_max,i_min:i_max]
			B_p = athena_dic['Bcc3'][k,j_min:j_max,i_min:i_max]
			MHDrun = True	
		except:
			MHDrun = False	

		for var in vars:
			try:
				if var == 'd':
					qty = d*units.d
				elif var == 'p':
					qty = p
				elif var == 'v_x':
					qty = v_x
				elif var == 'v_y':
					qty = v_y
				elif var == 'v_z':
					qty = v_z
			
				# vorticity quantities
				elif var in label_dic['vorticity_keys']: 
					if uov_key in dump_dirs.keys():
						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
						print("[load_vars]: uov_file = {}".format(uov_file))
						uov_dic = athena_read.vtk(uov_file)
						w_r = uov_dic[3]['user_out_var0'][0,:,0:i_max]
						w_t = uov_dic[3]['user_out_var1'][0,:,0:i_max]
						w_p = uov_dic[3]['user_out_var2'][0,:,0:i_max]
						wp = np.sqrt(w_r**2 + w_t**2)
						w = np.sqrt(w_r**2 + w_t**2 + w_p**2)
						if var == 'w_r':
							qty = w_r
						elif var == 'w_t':
							qty = w_t
						elif var == 'w_p':
							qty = w_p
						elif var == 'wp':
							qty = wp
						elif var == 'w':
							qty = w
					else:
						print("Post-processing vorticity calculation not yet implemented!")
						# sys.exit()
						
				# other derived hydro quantities				
				elif  var in label_dic['hydro_keys']: 
# 					print "[SimDics::loadVars()]: group1 = {}, var = {}".format(group1,var)
					if group1:
						vp = np.sqrt(v_x**2 + v_y**2)
						v = np.sqrt(v_x**2 + v_y**2 + v_z**2)
						group1 = 0
					
					if var == 'v': # poloidal velocity
						qty = v
					elif var == 'vp': # poloidal velocity
						qty = vp
					elif var == 'mfd': # mass flux density
						qty = d*vp
					elif var == 'cs': # adiabatic sound speed
						qty = cs
					elif var == 'c_iso': # adiabatic sound speed
						qty = units.cs*np.sqrt(p/d)
					elif var == 's': # entropy
						qty = np.log(p/d**gamma)
					elif var == 'T': # temperature
						qty = units.T*p/d
						# filter = qty > 1e-20
						# qty[filter] = 1.
					elif var == 'T_vir': # virial temperature
						qty = 3.*grid_dic['r_c_grid']*p/d
					elif var == 'Mp': # poloidal Mach number
						qty = vp/cs	
					elif var == 'l': # magnitude of specific angular momentum
						qty = grid_dic['r_c_grid']*v_p*grid_dic['sin_thetas']
					elif var == 'omega': # angular velocity
						qty = v_p/(grid_dic['r_c_grid']*grid_dic['sin_thetas'])
					elif var == 'Be': # Bernoulli function 
						qty = 0.5*v**2 + cs**2/(gamma - 1.) - \
						units.GM/grid_dic['r_c_grid'] 
					
				# MHD quantities	
				elif var in label_dic['mhd_keys']:
		
					if group2:
						try:
							B = np.sqrt(B_r**2 + B_t**2 + B_p**2)
							Bp = np.sqrt(B_r**2 + B_t**2) # poloidal field
							v_A = B/np.sqrt(d)
							v_Ap = Bp/np.sqrt(d)
							# compute fast and slow wave speeds, e.g. Shu Vol 2 page 307
							cos_phi = np.cos(Bp/B)
							vA2 = v_A**2
							cs2 = cs**2
							discrim = (vA2 + cs2)**2 - 4.*vA2*cs2*cos_phi**2
							sqrt_check = np.min(discrim)
							if sqrt_check < 0.:
								print("[SimDics.load_vars(): undefined sqrt! -> min(discrim) = {}".format(sqrt_check))
								print("[SimDics.load_vars(): taking an abs() meaning v_fast/v_slow are not defined!!")
								sqrt_term = np.sqrt(np.abs(discrim))
							else:
								sqrt_term = np.sqrt(discrim)
							v_fast = np.sqrt(0.5*(vA2 + cs2 + sqrt_term))
							v_slow = np.sqrt(0.5*(vA2 + cs2 - sqrt_term))
							group2 = 0
						except:
							if MHD_check: 
								continue
							else:
								raise ValueError(var)
					if var == 'B_r':
						qty = B_r
					elif var == 'B_t':
						qty = B_t
					elif var == 'B_p':
						qty = B_p
					elif var == 'v_A': # total Alven speed
						qty = v_A
					elif var == 'v_Ap': # poloidal Alven speed
						qty = v_Ap
					elif var == 'B': # magntude of B-field
						qty = B
					elif var == 'Bp':  # poloidal B-field
						qty = Bp
					elif var == 'E_r':
						qty = (v_p*B_t - v_t*B_p)/units.c
					elif var == 'E_t':
						qty = (v_r*B_p - v_p*B_r)/units.c
					elif var == 'E_p':
						qty = (v_t*B_r - v_r*B_t)/units.c
					elif var == 'Ep':
						qty = np.sqrt(E_r**2 + E_t**2)
					elif var == 'E':
						qty = np.sqrt(E_r**2 + E_t**2 + E_p**2)
					elif var == 'E_p_over_Ep':
						E_r = (v_p*B_t - v_t*B_p)/units.c
						E_t = (v_r*B_p - v_p*B_r)/units.c
						E_p = (v_t*B_r - v_r*B_t)/units.c
						qty = E_p/np.sqrt(E_r**2 + E_t**2)
					elif var == 'beta': # plasma beta
						qty = 2.*p/Bp**2
					elif var == 'B2': # plasma beta
						qty = B**2
					elif var == 'Bphi_over_Bp': 
						qty = B_p/Bp
					elif var == 'v_fast': 
						qty = v_fast
					elif var == 'v_slow': 
						qty = v_slow
					elif var == 'A': # See Contopoulos 1996 eq 14
						v_x = v_r*grid_dic['sin_thetas'] + v_t*grid_dic['cos_thetas']
						v_z = v_r*grid_dic['cos_thetas'] - v_t*grid_dic['sin_thetas']
						B_x = B_r*grid_dic['sin_thetas'] + B_t*grid_dic['cos_thetas']
						B_z = B_r*grid_dic['cos_thetas'] - B_t*grid_dic['sin_thetas']
						x_c = grid_dic['x_c']
						qty = x_c*(v_z*B_x - v_x*B_z)
					elif var == 'omega_A': # fastness parameter, e.g. Tzeferacos+13 eq 23
						vp = np.sqrt(v_r**2 + v_t**2)
						kappa = d*vp/Bp
						Omega = v_p/grid_dic['x_c'] - kappa*B_p/d/grid_dic['x_c']
# 						M_Ap = vp/v_Ap
# 						filter = np.where(np.abs(M_Ap - 1.)) < 1e-3
						qty = Omega*grid_dic['x_c']/v_Ap
					elif var == 'Fmag_par':
						grad_Bphi_r = calc_grad_r(grid_dic['x_c']*B_p,grid_dic['r_c_grid'])
						grad_Bphi_t = calc_grad_t(grid_dic['x_c']*B_p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_Bphi_par = n_r*grad_Bphi_r + n_t*grad_Bphi_t
						qty = -(B_p/grid_dic['x_c'])*grad_Bphi_par
					elif var == 'Fmag_phi':
						grad_Bphi_r = calc_grad_r(grid_dic['x_c']*B_p,grid_dic['r_c_grid'])
						grad_Bphi_t = calc_grad_t(grid_dic['x_c']*B_p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_Bphi_par = n_r*grad_Bphi_r + n_t*grad_Bphi_t
						qty = (Bp/grid_dic['x_c'])*grad_Bphi_par
					elif var == 'gradP_par':
						grad_p_r = calc_grad_r(p,grid_dic['r_c_grid'])
						grad_p_t = calc_grad_t(p,grid_dic['theta_c_grid'])
						n_r,n_t = B_r/Bp,B_t/Bp
						grad_p_par = n_r*grad_p_r + n_t*grad_p_t
						qty = -grad_p_par
					elif var == 'Be_mag':
						vp = np.sqrt(v_r**2 + v_t**2)
						kappa = d*vp/Bp
						Omega = v_p/grid_dic['x_c']
						Omega_mag = -kappa*B_p/d/grid_dic['x_c']
						Omega_B = Omega + Omega_mag
						magnetic = -(Omega_B/kappa)*grid_dic['x_c']*B_p
						kinetic_poloidal = 0.5*vp**2 
						kinetic_rotational = 0.5*v_p**2 
						enthalpy = cs**2/(gamma-1.)
						gravity = -units.GM/grid_dic['r_c_grid']
						qty = kinetic_poloidal + kinetic_rotational + enthalpy + gravity + magnetic
						
						
					if uov_key in dump_dirs.keys():
						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
						print("[load_vars]: uov_file = {}".format(uov_file))
						uov_dic = athena_read.vtk(uov_file)
						
						# unpack current density, j
						j_r = uov_dic[3]['user_out_var3'][0,:,0:i_max]
						j_t = uov_dic[3]['user_out_var4'][0,:,0:i_max]
						j_p = uov_dic[3]['user_out_var5'][0,:,0:i_max]

						# get rid of any nan's along north pole
						if np.abs(j_t[0,0]) == np.inf:
							j_t[0,:] = j_t[-1,0] # use south pole values on north pole
						
						jp = np.sqrt(j_r**2 + j_t**2)
						j = np.sqrt(j_r**2 + j_t**2 + j_p**2)
						

						
						# now compute the Lorentz force, jxB
						jcrossB_r = j_t*B_p - j_p*B_t
						jcrossB_t = j_r*B_p - j_p*B_r
						jcrossB_p = j_r*B_t - j_t*B_r
						jcrossBp = np.sqrt(jcrossB_r**2 + jcrossB_t**2)
						jcrossB = np.sqrt(jcrossB_r**2 + jcrossB_t**2 + jcrossB_p**2)
						
						if var == 'j_r':
							qty = j_r
						elif var == 'j_t':
							qty = j_t
						elif var == 'j_p':
							qty = j_p
						elif var == 'jp':
							qty = jp
						elif var == 'j':
							qty = j
						elif var == 'jcrossB_r':
							qty = jcrossB_r
						elif var == 'jcrossB_t':
							qty = jcrossB_t
						elif var == 'jcrossB_p':
							qty = jcrossB_p
						elif var == 'jcrossBp':
							qty = jcrossBp
						elif var == 'jcrossB':
							qty = jcrossB
					else:
						print("Post-processing current density calculation not yet implemented!")

				# quantities derived from optical depth 
				elif var in ['tau','xi']:
					if group3:
						tau = calc_optical_depth(grid_dic,d)
						group3 = 0
					if var == 'tau':
						qty = tau
					elif var == 'xi':  # photoionization parameter
						qty = units.xi*np.exp(-tau)/(d*grid_dic['r_c_grid']**2)	
						
				# flow Area diagnostics
				elif var in ['dlnA','div_v']:
					if var == 'dlnA':
						vp = np.sqrt(v_r**2 + v_t**2)
						mfd = d*vp
						grad_mfd_r = calc_grad_r(mfd,grid_dic['r_c_grid'])
						grad_mfd_t = calc_grad_t(mfd,grid_dic['theta_c_grid'])
						n_r,n_t = v_r/vp,v_t/vp
						dlnA = -(n_r*grad_mfd_r + n_t*grad_mfd_t)/mfd
						qty = dlnA
					elif var == 'div_v':
# 						p_avg = 0.5*(p[1:,1:] + p[:-1,:-1])
						qty = calc_div(v_r,v_t,grid_dic)
					group5 = 0
						

				# numerics: cell crossing times and floors - e.g., dt_cs_t = r(dtheta)/|cs|
				elif var in ['dt_cs_r','dt_cs_t','dt_v_r','dt_v_t','dt_vA_r','dt_vA_t','d_floor','p_floor']:
					if group4: # need all speeds
						cs = units.cs*np.sqrt(gamma*p/d)
						v = np.sqrt(v_r**2 + v_t**2 + v_p**2)

						dtheta_grid = grid_dic['theta_grid'][1:,:] - grid_dic['theta_grid'][0:-1,:] 
						rdtheta_grid = grid_dic['r_c_grid']*dtheta_grid[:,0:-1]
						dr_grid = grid_dic['r_grid'][:,1:] - grid_dic['r_grid'][:,0:-1]
						dr_grid = dr_grid[0:-1,:]

						group4 = 0

					if var == 'dt_cs_r': # sound crossing time in radial dir
						qty = dr_grid/cs
					elif var == 'dt_cs_t': # sound crossing time in polar dir
						qty = rdtheta_grid/cs
					if var == 'dt_v_r': # flow crossing time in radial dir
						qty = dr_grid/v
					elif var == 'dt_v_t': # flow crossing time in polar dir
						qty = rdtheta_grid/v
					elif var == 'dt_vA_r': # Alven crossing time in radial dir
						qty = dr_grid/v_A
					elif var == 'dt_vA_t': # Alven crossing time in polar dir
						qty = rdtheta_grid/v_A
					elif var == 'd_floor': # map of zones hitting dens floor
						qty = np.copy(d)
						qty[max_filter] = 1e20 # np.nan
					elif var == 'p_floor': # map of zones hitting pres floor
						qty = np.copy(p)
						qty[min_filter] = 1e20 #np.nan
					
				
				# scratch variable		
				elif var == 'my_var':
					if MHDrun:
						qty = my_var(grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p)
					else:
						qty = my_var(grid_dic,d,p,v_r,v_t,v_p)
				
				# Athena++ user output scalar variable
# 				elif var.startswith('user_out'): #and not(uov_vector):
# 					# open the vtk file using athena_read
# 					uov_key = idx + 'uov'
# 					if uov_key in dump_dirs.keys():
# 						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
# 						uov_dic = athena_read.vtk(uov_file)
# 						qty = uov_dic[3][var][0,:,0:i_max]
# 						qty = uov_scalar(qty,int(var[-1]))

				# Athena++ user output vector variable
				elif var.startswith('user_out'):
					print("\n\n USER_OUT_VAR!! \n\n")
					index = int(var[-1])
					var_base_name = var[:-1]
					uovs = []

					# open the vtk file using athena_read
					uov_key = idx + 'uov'
					if uov_key in dump_dirs.keys():
						uov_file = dump_dirs[uov_key] + str(i_dump).zfill(5) + '.vtk'
						print("[load_vars]: uov_file = {}".format(uov_file))
						uov_dic = athena_read.vtk(uov_file)
						for n in range(index+1):
							uov_var = var_base_name + np.str(n)
							print(uov_var)
							uovs.append(uov_dic[3][uov_var][0,:,0:i_max])
						if MHDrun:
							uovs = uov_vector(uovs,grid_dic,d,p,v_r,v_t,v_p,B_r,B_t,B_p)
						else:
							uovs = uov_vector(uovs,grid_dic,d,p,v_r,v_t,v_p)
						
						# if vector uov, add each component to VarDics
						if len(uovs) > 1:
							for n in range(index+1):
								uov_var = var_base_name + np.str(n)
								VarDics[idx][uov_var] = uovs[n]
							qty = VarDics[idx][var] 
						else:
							qty = uovs[0]

				else: # variable not recognized - throw exception 
					raise ValueError(var)	
			
				VarDics[idx][var] = qty

			except ValueError as err:
				label_dic = load_labels()
# 				print "[loadvars()]: label_dic['mhd_keys'] = {}".format(label_dic['mhd_keys'])
				if var in label_dic['mhd_keys']:
					print("\n******\nERROR!\n******")
					print("MHD variable {} was recognized but not found for simulation:\n{}".format(err.args[0],dump_file))
					print("This must not be an MHD simulation.")
					print("Choose from these hydro quantities:")
					print(label_dic['hydro_keys'])
				else:
					print("\n******\nERROR!\n******")
					print("Variable {} not recognized. Choose from:".format(err.args[0]))
					print("hydro vars:\n{}".format(label_dic['hydro_keys']))
					print("mhd vars:\n{}".format(label_dic['mhd_keys']))
					print("misc vars:\n{}".format(label_dic['misc_keys']))
				sys.exit()

	
	return VarDics

def load1DSimulationData(GridDic,SimDic,units,defaults=True):

	r = GridDic['x_c']
	simdata = {}
	simdata['X'] = r
	

	simdata['density'] = SimDic['d']
	simdata['velocity_x'] = SimDic['v_x']
	
	if defaults:
		gamma = 5./3.
		simdata['pressure'] = SimDic['p']
		simdata['temperature'] = SimDic['p']/SimDic['d']*units['T']
# 		simdata['temperature'] =  gamma*simdata['pressure']/simdata['density'] # local cloud sims
		simdata['xi'] = units['xi']/SimDic['d']/r**2
		
	print("[load1DSimulationData]:\nDone loading data.")
	return simdata
	
def load2DSimulationData(GridDic,SimDic,defaults=True):

	X = GridDic['x_c_grid']
	Y = GridDic['y_c_grid']
	simdata = {}
	simdata['X'] = X
	simdata['Y'] = Y
	

	simdata['density'] = SimDic['d']
	simdata['velocity_x'] = SimDic['v_x']
	simdata['velocity_y'] = SimDic['v_y']
	
	if defaults:
		gamma = 5./3.
		simdata['pressure'] = SimDic['p']
		simdata['temperature'] = gamma*SimDic['p']/SimDic['d']
	
	# copy SimDic data to simdata
# 	for key in SimDic.keys():
# 		simdata[key] = SimDic[key]
		
	print("[load2DSimulationData]:\nDone loading data.")
	return simdata