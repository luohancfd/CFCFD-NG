"""
Unittest to compare the behaviour of the Bender real gas model for CO2
with that of REFPROP.

Author: Peter Blyton
Version: 22/03/2012
"""

import unittest
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

class BenderTest(unittest.TestCase):
    def setUp(self):
        create_gas_file("Bender real gas", ["CO2"])
        self.gmodel = create_gas_model("gas-model.lua")
        self.gas = Gas_data(self.gmodel)
        set_molef(self.gas, self.gmodel, {'CO2':1})

    def test_press_from_rhoT(self):
        self.gas.rho = 268.58
        self.gas.T[0] = 300.0 # Sub critical, on the saturated vapor line
        self.gas.p = 0.0
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p, 6.713e6, delta=1.0e3)

    def test_temp_from_rhop(self):
        self.gas.rho = 40.0
        self.gas.T[0] = 100.0
        self.gas.p = 1.58e6 # Superheated vapour
        self.gmodel.eval_thermo_state_rhop(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 250.0, delta=0.1)

    def test_rho_from_pT(self):
        self.gas.rho = 0.0
        self.gas.T[0] = 350.0
        self.gas.p = 11.071e6 # Supercritical
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(self.gas.rho, 270.0, delta=0.3)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(BenderTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

