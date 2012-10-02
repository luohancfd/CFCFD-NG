# test_cea2_gas_h2ne.py
# Weed out troublesome CEA2 results from an example provided by Chris James.
# PJ, 02-Oct-2012

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name

# g = make_gas_from_name('h2ne', outputUnits='moles') # not same fractions
if 0:
    # This one works OK for cea2 
    g = Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles',
            outputUnits='massf', with_ions=True)
else:
    # but this one gives problems until today, specifically
    # "Cannot make a float from this string:  ******e-3"
    # The value for H+ should be around 0.0099998
    g = Gas(reactants={'H2':0.85, 'Ne':0.15}, inputUnits='moles',
            onlyList=['H2', 'H', 'H+', 'e-', 'Ne'],
            outputUnits='massf', with_ions=True)

rho = 2.754229e-05
e = 8314.0 * 1.548053e+04
g.set_rhoe(rho, e)
g.write_state(sys.stdout)


