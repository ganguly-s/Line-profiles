import os
import sys
import pickle
import array as ar


data_path = ''  # path to working directory

# flags -------------------------------------------------------------------

flag1D = 1			# for 1D data, use flag1D = 1

# optimization features 
RECALC_LUT = 0		# force a recalc of opacity LUT
RELOAD_GRID = 0		# force a recalc of hydro grid data
DO_PICKLE = 1		# set to 1 to pickle the fluxes and tau-profiles

# Hydro settings --------------------------------------------------------

# specify simulation dumps
# dumps = ar.array('i',(i for i in range(0,305,2)))
n_dump = 0

# specify paths to the runs 
ext = 'tab' # choices are 'athdf','vtk', and 'tab'
run_paths = {}; dump_tags = {}
uov_tags = None

run_paths['1'] = '/Users/brazen/Documents/projects/LPsRandall/1Druns/'

if 1: # Model A
	dumps = [1200]
	run_names = ['hep18'] 
	dump_tags['1'] = 'agn1.hep18.merged.out1.'

if 1: # Model B
	dumps = [1200]
	run_names = ['hep17'] 
	dump_tags['1'] = 'agn1.hep17.merged.out1.'

if 0: # high-res Model A 
	dumps = [1750]
	run_names = ['ddAx8'] 
	dump_tags['1'] = 'agn1.out1.'
		
if 0: # high-res Model B
	dumps = [1950] 
	run_names = ['ddBx8'] 
	dump_tags['1'] = 'agn1.out1.'

if 0: # low-res avg Model A
	dumps = [-1] 
	run_names = ['Ax1avg'] 
	dump_tags['1'] = 'ave-model-a-x1.out1.'
	
if 0: # high-res avg Model A
	dumps = [-1] 
	run_names = ['Ax8avg'] 
	dump_tags['1'] = 'ave-model-a-x8.out1.'
	
if 0: # low-res avg Model B
	dumps = [-1] 
	run_names = ['Bx1avg'] 
	dump_tags['1'] = 'ave-model-b-x1.out1.'
	
if 0: # high-res avg Model B
	dumps = [-1] 
	run_names = ['Bx8avg'] 
	dump_tags['1'] = 'ave-model-b-x8.out1.'

if 0: # low-res avg Model C
	dumps = [-1] 
	run_names = ['Cx1avg'] 
	dump_tags['1'] = 'ave-model-c-x1.out1.'
	
if 0: # high-res avg Model C
	dumps = [-1] 
	run_names = ['Cx8avg'] 
	dump_tags['1'] = 'ave-model-c-x8.out1.'

if 0: # low-res avg Model D
	dumps = [-1] 
	run_names = ['Dx1avg'] 
	dump_tags['1'] = 'ave-model-d-x1.out1.'
	
if 0: # high-res avg Model D
	dumps = [-1] 
	run_names = ['Dx8avg'] 
	dump_tags['1'] = 'ave-model-d-x8.out1.'
			
zoom_region = None
default_grid_name = run_names[0]
grid_path = run_paths['1'] + run_names[0] + '/'


# Opacity table settings --------------------------------------------------------

# Atomic data
line_data = {}

if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'fe_xxv_6.7keV.dat'
	line_data['ion'] = 'fe_xxvi'
	line_data['m_i'] = 52.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1.78344','1.77802'] # corresponding energies are [6.95196E+03,6.97317E+03] keV
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'fe_xxv_6.7keV.dat'
	line_data['ion'] = 'fe_xxv'
	line_data['m_i'] = 52.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1.8505'] #  this is the resonance K-alpha transition at 6.70003E+03 keV
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'fe_x_6375.dat'
	line_data['ion'] = 'fe_x'
	line_data['m_i'] = 52.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1.8505']
if 1:	
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'c_iv_1550.dat'
	line_data['ion'] = 'c_iv'
	line_data['m_i'] = 12.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1548.2','1550.77']
if 0:	
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'c_vi_33.dat'
	line_data['ion'] = 'c_vi'
	line_data['m_i'] = 12.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['33.7342','33.7396']	
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'o_vi_1034.dat'
	line_data['ion'] = 'o_vi'
	line_data['m_i'] = 16.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1031.91'] 
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'o_vii_22.dat'
	line_data['ion'] = 'o_vii'
	line_data['m_i'] = 16.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['21.602'] #,'21.8044','21.807','21.8044']	
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'o_viii_19.dat'
	line_data['ion'] = 'o_viii'
	line_data['m_i'] = 16.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['18.9671','18.9725']
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'he_ii_304.dat'
	line_data['ion'] = 'he_ii'
	line_data['m_i'] = 4.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['303.78','303.786']
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'mg_ii_2798.dat'
	line_data['ion'] = 'mg_ii'
	line_data['m_i'] = 24.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['2798.75']
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'n_v_1240.dat'
	line_data['ion'] = 'n_v'
	line_data['m_i'] = 14.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['1238.82','1242.8']
if 0:
	line_data['path'] = 'ion_opacities/xstar_data/'
	line_data['SED'] = 'AGN1'
	line_data['table'] = 'ne_viii_780.dat'
	line_data['ion'] = 'ne_viii'
	line_data['m_i'] = 20.
	line_data['A_Z'] = 1.
	line_data['lines'] = ['780.324']
			            

# Problem units --------------------------------------------------------
 
units = {}
units['xi'] = 1.37045e19
units['T'] = 7.268808495659315e-09  # this is mbar/kb with mbar = 0.6*m_p
units['n'] = 9.964388669908163e+23  # this is 1./mbar with mbar = 0.6*m_p
units['rho'] = 1.
units['L'] = 1.
units['v'] = 1.

# other parameters
dt = 1. 
vc = 0. # comoving frame velocity

