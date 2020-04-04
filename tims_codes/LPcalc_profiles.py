# import python libraries
import os
import sys
import array as ar
import matplotlib.animation as animation
import pickle
import importlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *

# import user libraries
import SimDics
import XstarInterface.LPmethods2D as XI

# load companion scripts
# import inputs; reload(inputs)
import imp
import inputs1D as inputs; imp.reload(inputs)
import fig_setup; imp.reload(fig_setup)

Fund = {'mp':1.6726231e-24, 'me':9.1093898e-28, 'c':3.0e10, 'h':6.6260755e-27, \
				'e':4.8032068e-10, 'kb':1.380658e-16, 'G':6.67259e-8}
				
# data -----------------------------------------------------------------------
# load Athena data

n_dump = inputs.n_dump # index of dump to use in array inputs.dumps

# get paths to sim output
dump_dirs = SimDics.dump_dir_dic(run_names=inputs.run_names,run_paths=inputs.run_paths,dump_tags=inputs.dump_tags,uov_tags=inputs.uov_tags)

# load grid dictionaries from sim(s)
try:
	grid_dics = pickle.load(open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "rb"))
except:
	if inputs.ext == 'tab':
		grid_dics = SimDics.generate_grid_data1D(dump_dirs,i_dump=inputs.dumps[n_dump])
	else:
		grid_dics = SimDics.generate_grid_data(dump_dirs,i_dump=inputs.dumps[n_dump],zoom_region=inputs.zoom_region,ext=inputs.ext)
	pickle.dump(grid_dics, open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "wb"), protocol=2)
finally:
	if inputs.RELOAD_GRID:
		if inputs.ext == 'tab':
			grid_dics = SimDics.generate_grid_data1D(dump_dirs,i_dump=inputs.dumps[n_dump])
		else:
			grid_dics = SimDics.generate_grid_data(dump_dirs,i_dump=inputs.dumps[n_dump],zoom_region=inputs.zoom_region,ext=inputs.ext)
		pickle.dump(grid_dics, open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "wb"), protocol=2)

# load variables
vars = ['d','T','p','v_x','v_y']

# loop over all outputs
for n_dump in range(len(inputs.dumps)):

	if inputs.ext == 'tab':
		sim_dics = SimDics.load_vars1D(inputs.dumps[n_dump],vars,dump_dirs,grid_dics)
		simdata = SimDics.load1DSimulationData(grid_dics['1'],sim_dics['1'],inputs.units)
	else:
		sim_dics = SimDics.load_vars(inputs.dumps[n_dump],vars,dump_dirs,grid_dics,gamma=5./3.,ext=inputs.ext)
		simdata = SimDics.load2DSimulationData(grid_dics['1'],sim_dics['1'])
	
	simdata['simfolder'] = inputs.run_names[0]
	simdata['gamma'] = 5./3.

	# add any net velocity (used with cloud tracking runs)
	vc = 0. #vc_data[dump,1]
	simdata['vc'] = vc

	X2D = simdata['X']
	d2D = simdata['density']
	T2D = simdata['temperature']

	line_data = inputs.line_data	# atomic line data


	# figure setup ---------------------------------------------------------------
	# make the figure template
	fig = figure(figsize=(6,7.))
	fig,axes = fig_setup.vel_space_plot(fig,simdata,line_data)

	time = inputs.dt*inputs.dumps[n_dump]
	time_years = round(time/(3600*24*365),1)
# 	axes['1'].set_title('time = ' + str(time) +  ' (' + str(time_years) + ' years)')
	axes['1'].set_title(inputs.run_names[0] + ' (dump = ' + str(time) + ')')
	print("time = {}, vc = {} km/s".format(time,vc*inputs.units['v']/1e5))

	# calculation -----------------------------------------------------------------
	# Line profile calcs

	# call to Xstar interface - LPcalc is a new instance of LineProfile class
	LPcalc = XI.LineProfiles(inputs.dumps[n_dump],inputs.run_names[0],simdata,inputs.units,line_data,Nwidths=1,Nv=1000,soln_is_1D=inputs.flag1D,recalc=inputs.RECALC_LUT) 

	# sub calc ---------------------------

	# calculate both optical depth and line profiles
	LPcalc.calc_final_profiles(do_pickle = inputs.DO_PICKLE) 
	Ib = LPcalc.line_profiles[line_data['lines'][0]] # LP of blue line in doublet	
	It = LPcalc.line_profiles['total']

	v = LPcalc.vel # LOS velocity in km/s

	# Alternatively, use wavelength units (untested)
	if 0: # wavelength axis in Angstroms
		lambda1 = float(line_data['lines'][0])/(1.+LPcalc.y) 
		lambda2 = float(line_data['lines'][0])/(1.+LPcalc.y)
		v1,v2 = lambda1,lambda2
		v = 0.5*(v1 + v2)
		Npts = len(v)
	else:
		v1 = v2 = v  # LOS velocity is the same for each line

	# sub calc ---------------------------
	# final optical depth profile as a function of LOS velocity
	# i.e. this is tau(nu) = (1/L) \int tau(x) dx, where L = \int dx
	tau_b = LPcalc.tau_profiles[line_data['lines'][0]]
	
	if len(line_data['lines']) > 1:
		Ir = LPcalc.line_profiles[line_data['lines'][1]] # LP of red line in doublet
		tau_r = LPcalc.tau_profiles[line_data['lines'][1]]



	# plots  ----------------------------------------------------------------------
	# -----------------------------------------------------------------------------

	# 1st panel ---------------------------
	axes['1'].semilogy(v,tau_b,'b-')
	if len(line_data['lines']) > 1:
		axes['1'].semilogy(v,tau_r,'r-')
	axes['1'].set_ylim(5e-3,1e3)
	axes['1'].set_ylabel(r'$\tau_\nu$')

	# 2nd panel ---------------------------
	axes['2'].plot(v,Ib,'b-')
	if len(line_data['lines']) > 1:
		axes['2'].plot(v,Ir,'r-')
		axes['2'].plot(v,It,color='purple',ls='-')
		axes['2'].plot(v,Ir-Ib,'k-')
		axes['1'].plot(v,Ir-Ib,'k-')
	axes['2'].set_ylim(0.,1.)
	axes['2'].set_ylabel(r'$I_\nu$')
		
	# 2D calcs only
# 	if inputs.SHOW_IMAGE:	# PLOT density image
# 		d_min,d_max = 0.12,8. #np.min(d2D),np.max(d2D)
# 		if inputs.Y_OBSERVER:
# 			im = axes['img'].imshow(d2D.transpose(),vmin=d_min, vmax=d_max, cmap='Spectral',extent=axes['img_box'])	
# 		else:
# 			im = axes['img'].imshow(d2D,vmin=d_min, vmax=d_max, cmap='Spectral',extent=axes['img_box'])
# 		cbar = plt.colorbar(im, cax=axes['cax'])
# 		fig_setup.fix_cbar_ticks(axes['cax'])


	# save fig --------------------------------------------------------------------
	# -----------------------------------------------------------------------------
	save_path = 'output/' + inputs.run_names[0] + '/' + LPcalc.ion + '/pngs/' 
	if not os.path.exists(save_path):
		os.makedirs(save_path)
	fig_title =  save_path + inputs.run_names[0] + '_profiles_' + str(inputs.dumps[n_dump]) + '.png'

	# fig.savefig('figs/LPcalc_profiles.png', dpi=160)
	fig.savefig(fig_title, dpi=160)
	plt.show()
