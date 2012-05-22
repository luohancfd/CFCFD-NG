#!/usr/bin/env python
"""
test_real_gas_REFPROP.py -- test the REFPROP gas model.

.. Author: Peter Blyton
.. Version: 21/05/2012
"""

import unittest
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

#gmodel = create_gas_file("real gas REFPROP", ["AIR.MIX", "single phase"])
#gmodel = create_gas_file("real gas REFPROP", ["AIR.PPF", "single phase"])
gmodel = create_gas_file("real gas REFPROP", ["R134A.FLD", "two phase"])

# Still need to set up a unittest for REFPROP
gas = Gas_data(gmodel)
gas.p = 101.3e3
gas.T[0] = 300.0
gas.rho = 10.347
gmodel.eval_thermo_state_rhoT(gas)
print gas.rho
print gas.quality
