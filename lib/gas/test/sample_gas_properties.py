# sample_gas_properties.py
# PJ 21-Sep-2010

import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from libprep3 import *

gmodel = create_gas_model('ideal-air.lua')
gas = Gas_data(gmodel)
gas.p = 10.935e3 # Pa
gas.T[0] = 4516.25 # degree K
gas.massf[0] = 1.0 # only one species
gmodel.eval_thermo_state_pT(gas)
gmodel.eval_transport_coefficients(gas)
gas.print_values()

