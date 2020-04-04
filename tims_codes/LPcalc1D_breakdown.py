# import python libraries
import os
import sys
import array as ar
import pickle
import matplotlib
import importlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *

# import user libraries
import SimDics # Athena interface
import XstarInterface.LPmethods2D as XI

# load companion scripts
import imp
import inputs1D as inputs; imp.reload(inputs)
import fig_setup; imp.reload(fig_setup)


# data -----------------------------------------------------------------------
# load Athena data

n_dump = inputs.n_dump # index of dump to use in array inputs.dumps

# get paths to both sims
dump_dirs = SimDics.dump_dir_dic(run_names=inputs.run_names,run_paths=inputs.run_paths,dump_tags=inputs.dump_tags,uov_tags=inputs.uov_tags)

# load grid dictionaries from sim(s)
try:
	grid_dics = pickle.load(open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "rb"))
except:
	grid_dics = SimDics.generate_grid_data1D(dump_dirs,i_dump=inputs.dumps[n_dump],ext=inputs.ext)
	pickle.dump(grid_dics, open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "wb"), protocol=2)
finally:
	if inputs.RELOAD_GRID:
		grid_dics = SimDics.generate_grid_data1D(dump_dirs,i_dump=inputs.dumps[n_dump],ext=inputs.ext)
		pickle.dump(grid_dics, open(inputs.grid_path + inputs.default_grid_name + '_GridDics.p', "wb"), protocol=2)

# load variables
vars = ['d','T','p','v_x']
if inputs.ext == 'tab':
	sim_dics = SimDics.load_vars1D(inputs.dumps[n_dump],vars,dump_dirs,grid_dics)
	simdata = SimDics.load1DSimulationData(grid_dics['1'],sim_dics['1'],inputs.units)
else:
	sim_dics = SimDics.load_vars(inputs.dumps[n_dump],vars,dump_dirs,grid_dics,gamma=5./3.,ext=inputs.ext)
	simdata = SimDics.load2DSimulationData(grid_dics['1'],sim_dics['1'])


simdata['simfolder'] = inputs.run_names[0]
simdata['output_folder'] = simdata['simfolder'] + '/' + inputs.line_data['ion']
simdata['gamma'] = 5./3.
simdata['vc'] = inputs.vc

line_data = inputs.line_data

# figure setup ---------------------------------------------------------------
# make the figure template
fig = figure(figsize=(6,10))
fig,axes = fig_setup.pos_space_plot(fig,simdata,line_data)

time = inputs.dt*inputs.dumps[n_dump]
time_years = round(time/(3600*24*365),0)
# axes['1'].set_title('time = ' + str(time) +  ' (' + str(time_years) + 'years)')
axes['1'].set_title(inputs.run_names[0] + ' ' + line_data['ion'] + ' (dump = ' + str(time) + ')')
# print("time = {}, vc = {} km/s").format(time,vc*eqvals['c0']/1e5)

# calculation -----------------------------------------------------------------
# Line profile calcs

# call to Xstar interface - LPcalc is a new instance of LineProfile class
LPcalc = XI.LineProfiles(inputs.dumps[n_dump],inputs.run_names[0],simdata,inputs.units,line_data,Nwidths=1,Nv=1000,soln_is_1D=inputs.flag1D,recalc=inputs.RECALC_LUT) 

# sub calc ---------------------------
# calculate both optical depth and line profiles
LPcalc.calc_final_profiles(do_pickle = inputs.DO_PICKLE) # includes LP calc
Ib = LPcalc.line_profiles[line_data['lines'][0]] # this is the LP 


# sub calc ---------------------------
# optical depth profile at line center vs radius, i.e. tau(r) = \int alpha0 dr
tau_y_lc = LPcalc.tau_xory_profiles[line_data['lines'][0]][int(LPcalc.Nv/2)]
tau_lcs = np.zeros(LPcalc.Nv)
for i in range(LPcalc.Nv):
	tau_lcs[i] = LPcalc.tau_xory_profiles[line_data['lines'][0]][i][0]



# calculation -----------------------------------------------------------------
# unpack T/xi from hydro profiles 

y1 = simdata['X']
Ny = simdata['X'].shape

d1 = simdata['density']
d1_cgs = inputs.units['rho']*d1

T1 = simdata['temperature']
T1_cgs = T1
T5 = T1_cgs/1e5

xi1 = simdata['xi']
xmin,xmax = np.min(simdata['X']),np.max(simdata['X'])

xmin,xmax = np.round(xmin),np.round(xmax)

# get the opacties
los_opacities_blue,los_opacities_red = LPcalc.get_LOS_opacities1D()
# rb_ratio0 = los_opacities_blue[0]/los_opacities_red[0]
Nvals = len(los_opacities_blue[0])

if 1:  # plot phi and alpha(nu) (test for Shalini)
	fig0, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
	phi0 = LPcalc.phi_nu(0.)
	phi_p = LPcalc.phi_nu(1e-4)
	phi_n = LPcalc.phi_nu(-1e-4)
	
	ax1.semilogx(y1,phi0,label='$y=0$')
	ax1.semilogx(y1,phi_p,label='$y=10^{-4}$')
	ax1.semilogx(y1,phi_n,ls='--',label='$y=-10^{-4}$')
	ax1.set_ylabel(r'$\phi(y)$')
	
	phi_dic = {}
	for key,val in zip(['r','0','p','n'],[y1,phi0,phi_p,phi_n]):
		phi_dic[key] = val
	pickle.dump(phi_dic, open('phi_dic.p', "wb"), protocol=2)
	
	alpha0_b = LPcalc.alpha_nu(phi0,los_opacities_blue)[0]
	alphap_b = LPcalc.alpha_nu(phi_p,los_opacities_blue)[0]
	alphan_b = LPcalc.alpha_nu(phi_n,los_opacities_blue)[0]
# 	alpha0_r = LPcalc.alpha_nu(phi0,los_opacities_red)[0]
# 	alphap_r = LPcalc.alpha_nu(phi_p,los_opacities_red)[0]
# 	alphan_r = LPcalc.alpha_nu(phi_n,los_opacities_red)[0]

	alpha_dic = {}
	for key,val in zip(['r','0','p','n'],[y1,alpha0_b,alphap_b,alphan_b]):
		alpha_dic[key] = val
	pickle.dump(alpha_dic, open('alpha_dic.p', "wb"), protocol=2)
	
	alpha0 = 1e-15
	ax2.loglog(y1,alpha0_b/alpha0)
	ax2.loglog(y1,alphap_b/alpha0)
	ax2.loglog(y1,alphan_b/alpha0)
	ax2.set_ylim(1e-10,1e1)
	xlabel('r')
	ylabel(r'$\alpha(y)/10^{-15}$')
	
	ax1.legend()
	show()

if 0: # plot of line-center optical depth vs. LOS velocity (test for Shalini)
	fig0, (ax1) = plt.subplots(nrows=1, ncols=1)
	ax1.plot(LPcalc.vel,tau_lcs)
	ax1.set_xlabel(r'$v_{los}$' + ' ' + '[km/s]')
	ax1.set_ylabel(r'$\tau_{\nu_0}(r_{out})$')


# plots  ----------------------------------------------------------------------
# -----------------------------------------------------------------------------
axes['1'].set_xlim(xmin,xmax)

# 1st panel ---------------------------
axes['1'].semilogy(y1,d1,'.k-')
# axes['1'].set_ylim(0.1,15.)
axes['1'].set_ylabel(r'$n$')

# 2nd panel ---------------------------
axes['2'].semilogy(y1,xi1,'k-')
axes['2b'].semilogy(y1,T1_cgs,'r--')
axes['2'].set_ylim(1e1,1e4) 
axes['2b'].set_ylim(1e4,5e6)
axes['2'].set_ylabel(r'$\xi$')
axes['2b'].set_ylabel(r'$T$')

# 3rd panel ---------------------------
axes['3'].semilogy(y1,los_opacities_blue[0],'b-',label=r'$\kappa_{blue}$')
# axes['3'].semilogy(y1,los_opacities_red[0],'r-',label=r'$\kappa_{red}$')
# axes['3'].set_xlim(0.,1.)
axes['3'].set_ylim(1e-23,1e-14)
axes['3'].set_ylabel(r'$\alpha_{\nu_0}$')
axes['3'].set_xlabel(r'$r$')

# 4th panel ---------------------------
# velocity vs. distance 
axes['4'].plot(y1,simdata['velocity_x']/1e5,'k-')
axes['4'].set_ylim(1e0,6.5e2)

# add legends
legend = axes['3'].legend(loc='best', shadow=False, fancybox=True, framealpha = 0.1)


# save fig --------------------------------------------------------------------
# -----------------------------------------------------------------------------
save_path = 'output/' + simdata['output_folder'] + '/pngs/' 
if not os.path.exists(save_path):
	os.makedirs(save_path)
fig_title =  save_path + inputs.run_names[0] + '_breakdown_' + str(inputs.dumps[n_dump]) + '.png'

fig.savefig(fig_title, dpi=160)
plt.show()
