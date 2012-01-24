# sample_thermally_perfect-air.py
# PJ 05-Oct-2011
# To explore problems with PJ's port of Brendan's thesis code.

import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from libprep3 import *

gmodel = create_gas_model('thermally-perfect-air.lua')
gas = Gas_data(gmodel)
gas.p = 10.935e3 # Pa
gas.T[0] = 4516.25 # degree K
set_molef(gas, gmodel, {'O2':0.21, 'N2':0.79})
gmodel.eval_thermo_state_pT(gas)
gmodel.eval_transport_coefficients(gas)
gas.print_values()

print 50*'-'
print "Make a new object, presumably using the C++ copy constructor."
gas2 = Gas_data(gas)
gas2.print_values()

print "Done."


