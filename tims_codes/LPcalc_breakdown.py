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
import AthenaInterface.AthenaUtils as AI
import XstarInterface.LPmethods2D as XI

# load companion scripts
import inputs; reload(inputs)
import fig_setup; reload(fig_setup)


# data -----------------------------------------------------------------------
# load Athena data

n_dump = inputs.n_dump # index of dump to use in array inputs.dumps

simdata = AI.load2DSimulationData(inputs.dumps[n_dump],inputs.Sim1)

# add any net velocity (used with cloud tracking runs)
vc = 0. #vc_data[dump,1]
if inputs.Y_OBSERVER:
	vc = 0.
simdata['vc'] = vc

X2D = simdata['X']
d2D = simdata['density']
T2D = simdata['temperature']

# abbreviate commonly used data in inputs
eqvals = inputs.eqvals			# physical units
line_data = inputs.line_data	# atomic line data
xi0 = inputs.xi0				# ionization parameter



# figure setup ---------------------------------------------------------------
# make the figure template
inputs.SHOW_IMAGE = 1
fig,axes = fig_setup.breakdown_plot(inputs.SHOW_IMAGE)

# set x-limits 
# (applies to ax1, ax2, & ax3 since ax2/ax3 share y-axis of ax1)
axes['1'].set_xlim(0.,1.)

time = inputs.Sim1['dt']*inputs.dumps[n_dump]
time_hours = round(time*eqvals['t_th']/(3600*24),1)
axes['1'].set_title('time = ' + str(time) +  ' (' + str(time_hours) + ' days)')
print("time = {}, vc = {} km/s").format(time,vc*eqvals['c0']/1e5)

# calculation -----------------------------------------------------------------
# Line profile calcs

# call to Xstar interface - LPcalc is a new instance of LineProfile class
LPcalc = XI.LineProfiles(inputs.dumps[n_dump],simdata,eqvals,line_data,Nwidths=1,Nv=100,y_observer=inputs.Y_OBSERVER) 

# sub calc ---------------------------
# calculate both optical depth and line profiles
LPcalc.calc_final_profiles(recalc=0) # includes LP calc
Ib = LPcalc.line_profiles[line_data['lines'][0]] # this is the LP

# sub calc ---------------------------
# final optical depth profile as a function of LOS velocity
# i.e. this is tau(nu) = (1/L) \int tau(x) dx, where L = \int dx
tau_b = LPcalc.tau_profiles[line_data['lines'][0]]

# gather some indices of special optical depth values
i_taumax = tau_b.argmax() # index of tau_max
if line_data['ion'] == 'OVIII': 
	tau_filter = np.where(tau_b > 1.)
	tau_filter = tau_filter[0]
	i_tau1 = tau_filter[0]	# index of tau = 1
	
	tau_filter = np.where(tau_b > 0.5)
	tau_filter = tau_filter[0]
	i_taupt5 = tau_filter[0] # index of tau = 0.5

# sub calc ---------------------------
# optical depth profiles as a function of position for 
# several LOS velocities, i.e. each is tau(x) = \int alpha dx
# 1. optical depth profile at line center
tau_y_lc = LPcalc.tau_xory_profiles[line_data['lines'][0]][LPcalc.Nv/2]
# 2. optical depth profile at LOS velocity where tau_b = 1
tau_y_1 = LPcalc.tau_xory_profiles[line_data['lines'][0]][i_tau1]
# 3. optical depth profile at LOS velocity where tau_b = 0.5
tau_y_pt5 = LPcalc.tau_xory_profiles[line_data['lines'][0]][i_taupt5]



# calculation -----------------------------------------------------------------
# hydro profiles - slices of den/temp maps through center of domain

y1 = X2D[0,:]
y1 = y1[::-1]
Ny,Nx = X2D.shape

if inputs.Y_OBSERVER:
	d1 = d2D[:,Nx/2-1]
	d2D_8 = eqvals['n0']*d2D/1e8
	d1_cgs = eqvals['n0']*d1
	
	T1 = T2D[:,Nx/2-1]
	T1_cgs = eqvals['T0']*T1
	T5 = T1_cgs/1e5
	
	xi1 = xi0/d1
else:
	d1 = d2D[Ny/2-1,:]
	d2D_8 = eqvals['n0']*d2D/1e8
	d1_cgs = eqvals['n0']*d1
	
	T1 = T2D[Ny/2-1,:]
	T1_cgs = eqvals['T0']*T1
	T5 = T1_cgs/1e5
	
	xi1 = xi0/d1



# calculation -----------------------------------------------------------------
# covering fraction profiles 
# plot from ipython using
# >> plot(v,C1); plot(v,C2)
	
Cs = LPcalc.calc_covering_fractions_exact()
C1 = Cs[line_data['lines'][0]] # 1st line in doublet
C2 = Cs[line_data['lines'][1]] # 2nd line in doublet


# plots  ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 1st panel ---------------------------
axes['1'].plot(y1,d1,'.k-')
axes['1'].set_ylim(0.1,3.)

# 2nd panel ---------------------------
axes['2'].plot(y1,xi1,'k-')
axes['2'].plot(y1,T5,'r-')
axes['2'].set_ylim(5e1,4.5e2) 
axes['2b'].set_ylim(np.min(T5),np.max(T5))

# 3rd panel ---------------------------
if inputs.Y_OBSERVER:
	los_opacities_blue,los_opacities_red = LPcalc.get_LOS_opacities([Ny/2-1],y_dir=True)
else:
	los_opacities_blue,los_opacities_red = LPcalc.get_LOS_opacities([Nx/2-1],y_dir=False)
rb_ratio0 = los_opacities_blue[0]/los_opacities_red[0]
Nvals = len(los_opacities_red[0])

axes['3'].plot(y1,los_opacities_blue[0],'b-',label=r'$\kappa_{blue}$')
axes['3'].plot(y1,los_opacities_red[0],'r-',label=r'$\kappa_{red}$')
axes['3'].set_xlim(0.,1.)
axes['3'].set_ylim(1e-13,5e-10)

# add legend
legend = axes['3'].legend(loc='best', shadow=False, fancybox=True, framealpha = 0.1)

# 4th panel --------------------------- 
# optical depths vs. distance along edge of domain

axes['4'].plot(y1,tau_y_lc,'.k-')
axes['4'].plot(y1,np.ones_like(y1)*tau_b[LPcalc.Nv/2-1],'b-',label=r'$\tau(v=0)$')
axes['4'].set_ylim(0,20.)

# add legend
legend = axes['4'].legend(loc='best', shadow=False, fancybox=True, framealpha = 0.1)

		
if 0: # overplot vertical lines marking covering region
	covered = np.where(tau_y_lc > 1.)
	covered = covered[0]
	axes['4'].axvline(y1[covered[0]], ymin=1e-12, ymax=1e6, color='k', ls='--', lw=2)
	axes['4'].axvline(y1[covered[-1]], ymin=1e-12, ymax=1e6, color='k', ls='--', lw=2)
	covered = np.where(tau_y_pt5 > 0.5)
	covered = covered[0]
	axes['4'].axvline(y1[covered[0]], ymin=1e-12, ymax=1e6, color='k', ls='--', lw=1)
	axes['4'].axvline(y1[covered[-1]], ymin=1e-12, ymax=1e6, color='k', ls='--', lw=1)

	# in ipython, view optical depth profile by doing
	# >> plot(y1, tau_y_pt5)
	
	
if inputs.SHOW_IMAGE:	# PLOT density image
	d_min,d_max = 3e-1,1.2e0
	if inputs.Y_OBSERVER:
		im = axes['img'].imshow(d2D_8.transpose(),vmin=d_min, vmax=d_max, cmap='Spectral',extent=axes['img_box'])	
	else:
		im = axes['img'].imshow(d2D_8,vmin=d_min, vmax=d_max, cmap='Spectral',extent=axes['img_box'])
	cbar = plt.colorbar(im, cax=axes['cax'])
	fig_setup.fix_cbar_ticks(axes['cax'])


# save fig --------------------------------------------------------------------
# -----------------------------------------------------------------------------
save_path = 'output/' + inputs.Sim1['output_folder'] + '/pngs/' 
if not os.path.exists(save_path):
	os.makedirs(save_path)
fig_title =  save_path + 'breakdown_' + str(inputs.dumps[n_dump]) + '.png'

fig.savefig(fig_title, dpi=160)
plt.show()
