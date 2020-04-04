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
import imp
import inputs1D as inputs; imp.reload(inputs)
import fig_setup; imp.reload(fig_setup)


# specify runs to compare
n_dumps = [1200,-1]
run_names = ['hep17','Bx1avg'] 
dump_tags = [{'1':'agn1.hep17.merged.out1.'},{'1':'ave-model-b-x1.out1.'}]
colors = ['c','k']

# specify opacity table
line_data = inputs.line_data	# atomic line data

# figure setup ---------------------------------------------------------------
# make the figure template
fig = figure(figsize=(6,7.))
fig,axes = fig_setup.vel_space_plot1D(fig,line_data)
					
# load Athena data -----------------------------------------------------------------------


# load variables
vars = ['d','T','p','v_x','v_y']

# loop over all outputs
for n_dump,run_name,dump_tag,col in zip(n_dumps,run_names,dump_tags,colors):
	
	# unpack previously calculated quantities 
	default_paths = XI.get_default_paths(n_dump,run_name,line_data['ion'])
	flux_dic = pickle.load(open(default_paths['fluxes'],"rb"))
	tau_dic = pickle.load(open(default_paths['taus'],"rb"))
	
	v = flux_dic['vel_kms']
	Ib = flux_dic[line_data['lines'][0]]
	It = flux_dic['total']
	tau_b = tau_dic[line_data['lines'][0]]
	tau_t = tau_dic['total']
	
	if len(line_data['lines']) > 1:
		Ir = flux_dic[line_data['lines'][1]] # LP of red line in doublet
		tau_r = tau_dic[line_data['lines'][1]]
	
	# plots  ----------------------------------------------------------------------
	# -----------------------------------------------------------------------------	
	
	
	# 1st panel ---------------------------
	axes['1'].semilogy(v,tau_b,ls='-',color=col,label=run_name)
	if len(line_data['lines']) > 1:
		axes['1'].semilogy(v,tau_r,'r-')
	axes['1'].set_ylim(5e-3,1e3)
	axes['1'].set_ylabel(r'$\tau_\nu$')

	# 2nd panel ---------------------------
	axes['2'].plot(v,Ib,ls='-',color=col)
	if len(line_data['lines']) > 1:
		axes['2'].plot(v,Ir,'r-')
		axes['2'].plot(v,It,color='purple',ls='-')
		axes['2'].plot(v,Ir-Ib,'k-')
		axes['1'].plot(v,Ir-Ib,'k-')
	axes['2'].set_ylim(0.,1.)
	axes['2'].set_ylabel(r'$I_\nu$')
	
	axes['1'].legend()


# save fig --------------------------------------------------------------------
# -----------------------------------------------------------------------------
save_path = 'output/' + inputs.run_names[0] + '/' + line_data['ion'] + '/pngs/' 
if not os.path.exists(save_path):
	os.makedirs(save_path)
fig_title =  save_path + inputs.run_names[0] + '_comparison.png'

# fig.savefig('figs/LPcalc_profiles.png', dpi=160)
fig.savefig(fig_title, dpi=160)
plt.show()
