#!/usr/bin/env python
"""
test_Bender_real_gas.py -- unittest to compare the behaviour
of the Bender real gas model for CO2 with results from REFPROP.

.. Author: Peter Blyton
.. Version: 22/03/2012
"""

import unittest
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

# Accepted percentage errors / 100
TEMP_TOL = 0.002
PRES_TOL = 0.0001
DENS_TOL = 0.002
ENER_TOL = 0.003
ENTR_TOL = 0.003

class BenderTest(unittest.TestCase):
    def setUp(self):
        create_gas_file("Bender real gas", ["CO2"])
        self.gmodel = create_gas_model("gas-model.lua")
        self.gas = Gas_data(self.gmodel)
        set_molef(self.gas, self.gmodel, {'CO2':1})

    def test_temp_from_rhoe(self):
        self.gas.rho = 46.644
        self.gas.e[0] = 398.77e3  # Sub critical, on the saturated vapor line
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 250.0, delta=TEMP_TOL*self.gas.T[0])
        self.assertAlmostEqual(self.gas.p, 1785.0e3, delta=ENER_TOL*self.gas.p)
        s = self.gmodel.total_entropy(self.gas)
        self.assertAlmostEqual(s, 1.9641e3, delta=ENTR_TOL*s)

    def test_rho_from_pT(self):
        self.gas.T[0] = 350.0
        self.gas.p = 11.071e6 # Supercritical
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(self.gas.rho, 270.0, delta=DENS_TOL*self.gas.rho)
        self.assertAlmostEqual(self.gas.e[0], 410.85e3, delta=ENER_TOL*self.gas.e[0])
        s = self.gmodel.total_entropy(self.gas)
        self.assertAlmostEqual(s, 1.7719e3, delta=ENTR_TOL*s)

    def test_press_from_rhoT(self):
        self.gas.rho = 268.58
        self.gas.T[0] = 300.0 # Sub critical, on the saturated vapor line
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p, 6.713e6, delta=PRES_TOL*self.gas.p)
        self.assertAlmostEqual(self.gas.e[0], 362.09e3, delta=ENER_TOL*self.gas.e[0])
        s = self.gmodel.total_entropy(self.gas)
        self.assertAlmostEqual(s, 1.6215e3, delta=ENTR_TOL*s)

    def test_temp_from_rhop(self):
        self.gas.rho = 40.0
        self.gas.p = 1.58088e6 # Superheated vapour
        self.gmodel.eval_thermo_state_rhop(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 250.0, delta=TEMP_TOL*self.gas.T[0])
        self.assertAlmostEqual(self.gas.e[0], 401.88e3, delta=ENER_TOL*self.gas.e[0])
        s = self.gmodel.total_entropy(self.gas)
        self.assertAlmostEqual(s, 2.0004e3, delta=ENTR_TOL*s)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(BenderTest)
    unittest.TextTestRunner(verbosity=2).run(suite)
