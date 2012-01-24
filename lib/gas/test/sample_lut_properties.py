# sample_lut_properties.py
# Try out the LUT model of air at rather high temperatures that seem to be
# typical of X-tube simulations.
#
# In the past, I had been using the -extrapolate option when generating
# tables for the gas model.  This has bitten me, yet again, because the
# tabulated data seems to be poorly extrapolated and the iterative solve
# for thermo properties given p and T seems to get stuck and then give
# really crap estimates.  The symptom in L1d3 was then poor Riemann solver
# output which was then passed onto the time update of momentum and energy.
# That crappy update, in turn, made a mess of the decoded properties for 
# the cell where we ended up with crazy values.
#
# This particular test p=180679 Pa, T=30018 K for air is an example of the
# bad data being fed to the EOS-from-pT function.  That function barfing
# quietly on the side and then providing a value of gamma very close to 1.0
# was leading to the generation of inf and nan numbers that propagated 
# through the L1d3 solution.  
#
# So, after all that, the principal cure is to NOT use -extrapolate when
# building look-up tables with the CEA code.  And, since it burned my 
# entire Sunday, I have written this little script to remind me to test...
#
# PJ 16-Jan-2011

import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from libprep3 import *

gmodel = create_gas_model('cea-lut-air.lua.gz')
gas = Gas_data(gmodel)
gas.p = 180679.298101  # Pa
gas.T[0] = 30017.9892315 # degree K
gas.massf[0] = 1.0 # only one species
gmodel.eval_thermo_state_pT(gas)
gmodel.eval_transport_coefficients(gas)
gas.print_values()
print "gamma=", gmodel.gamma(gas)

