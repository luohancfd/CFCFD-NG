#! /usr/bin/env python
# test_setting_moles.py
# PJ & Umar 12-Feb-2011

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from libprep3 import *

species=['O', 'N', 'N2', 'O2', 'NO', 'N_plus', 'O_plus', 'N2_plus',
         'O2_plus', 'NO_plus', 'e_minus', 'Ar', 'Ar_plus']
create_gas_file(model='thermally perfect gas', species=species, fname='gas-model.lua')

# Now, pick up the gas-model file and use it...
gmodel = create_gas_model('gas-model.lua') 
gas = Gas_data(gmodel)

# The gas conditions come from Umar's standing-shock experiment design for X2
# and the mole fractions from RGM's workbook.
gas.p = 4464.0 # Pa
gas.T[0] = 10140.42 # degree K
X = {'O':1.6936e-1, 'N':5.9784e-1, 'N2':6.9757e-5, 'O2':4.7543e-8, 'NO':2.5654e-3, 
     'N_plus':9.6331e-2, 'O_plus':1.7562e-2, 'N2_plus':7.7688e-6, 'O2_plus':5.0837e-8,
     'NO_plus':1.4459e-5, 'e_minus':1.1436e-1, 'Ar':4.0026e-3, 'Ar_plus':4.4835e-4}
gmodel.set_molef(gas, X)

gmodel.eval_thermo_state_pT(gas)
gmodel.eval_transport_coefficients(gas)
gas.print_values()

