# Specify desired unit conversion from code units
import numpy as np
# import athenaSetup_blondin as inputs; reload(inputs)
#
""" Intermediate conversions (not used anywhere but here) """

# import athenaSetup_ZUnits as inputs; reload(inputs)

# # Specify any parameters that determine units
# eps = inputs.eps # disk scale height 

# # conversion to cgs units
# xi_i = inputs.xi_i


""" Actual units used by load_vars() in SimDics """
d = 1.  # density units
v = 1.  # r-velocity units
p = 1.  # pressure/energy units
R = 1.  # distance units
T = 1.  # temperature units
xi = 1.
cs = 1. 

kappa = 15772099.50128901 #inputs.kappa_iso