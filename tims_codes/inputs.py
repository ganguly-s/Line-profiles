import os
import sys
import pickle
import array as ar
import PW14classes.ti as TImethods


data_path = '/home/tim/Documents/projects/LPs/LPs2020/'
#---------------------------------------------------------------------
# INPUTS

Y_OBSERVER = 0		# observer is in y-dir instead of x-dir
SHOW_IMAGE = 1		# plot image of simulation on top panel
RECALC = 0			# force a recalc of LPs (previous calc loaded by default)

# specify simulation dumps
# dumps = ar.array('i',(i for i in range(0,305,2)))
dumps = [20,27]
n_dump = 0			# dump index; n_dump = 0 will use dumps[0]

# Atomic data
line_data = {}
if 1:	
	line_data['ion'] = 'OVIII'
	line_data['m_i'] = 16.
	line_data['A_Z'] = 0.7820E-03
	line_data['xstar_file'] = 'ion_opacities/xstar_data/o_viii_19.dat'
	line_data['lines'] = ['18.9671','18.9725']
	line_data['fs'] = [2.771e-01,1.385e-01]

# Sim template --------------------------------------------------------
# Generic dictionary for Athena data
Sim = {}

Sim['simtag'] = 'ti'
Sim['filetype'] = 'tab'
Sim['nametag'] = ''
Sim['simfile1'] = 0
Sim['simfileN'] = 9999 #184
Sim['gamma'] = 5./3.

# Actual sim data --------------------------------------------------------

# Resolution
Nx = 1024
Ny = 1024

Sim1 = Sim.copy()
Sim1['simpath'] = data_path + 'sim_data/cRFLDX_pdv5/'
Sim1['simfolder'] = 'noPD_square/'
Sim1['imgfolder'] = 'yt_dat/'
Sim1['dt'] = 1.

# output_folder within dir 'output'
Sim1['output_folder'] = Sim1['simfolder'] + line_data['ion']  

if 0:  # load bulk velocities (cloud tracking runs)
	Sim1['vc_file'] = 'vc_data/pd0.2_amp0.2_veldata.tab'
	vc_data = np.loadtxt(Sim1['vc_file'])
	vc_data[frames[0],1]

# store all sims in an array
Sims = [Sim1]

# Get physical units
xi0 = 190.
Blondin = TImethods.HCcurves(Sim1['gamma'],nT=1e13)
eqvals = Blondin.getEquilibriumValues(xi0)
athena = Blondin.getCodeUnits(eqvals)
