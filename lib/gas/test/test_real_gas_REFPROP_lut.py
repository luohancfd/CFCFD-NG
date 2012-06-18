#!/usr/bin/env python
"""
test_real_gas_REFPROP_lut.py -- test the REFPROP gas model look-up-table by
comparing with results from REFPROP.

.. Author: Peter Blyton
.. Version: 18/06/2012
"""

import unittest
import sys, os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

class REFPROPlutTestR134A_FLD(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_model("REFPROP-lut-R134A.FLD.lua.gz")
        self.gas = Gas_data(self.gmodel)

    def test_props_from_rhoe(self):
        self.gas.rho = 1.0
        self.gas.e[0] = 443.153781283e3
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 350.0, places=2)
        self.assertAlmostEqual(self.gas.p/1e6, 0.0284307497818, places=5)
        self.assertAlmostEqual(self.gmodel.Cv(self.gas)/1e3, 0.842546300634, places=2)
        self.assertAlmostEqual(self.gmodel.Cp(self.gas)/1e3, 0.925398673641, places=2)
        self.gmodel.eval_transport_coefficients(self.gas)
        self.assertAlmostEqual(self.gas.mu/1e6, 13.8232015201, places=3)
        self.assertAlmostEqual(self.gas.k[0]*1e3, 17.5169869898, places=3)

if __name__ == "__main__":
    os.system("build-REFPROP-lut.py --fluid=R134A.FLD --bounds='300.0,450.0,-1.0,2.0' --divisions='250,250'")
    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPlutTestR134A_FLD)
    unittest.TextTestRunner(verbosity=2).run(suite)
