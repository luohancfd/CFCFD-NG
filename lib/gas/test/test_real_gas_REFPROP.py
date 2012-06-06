#!/usr/bin/env python
"""
test_real_gas_REFPROP.py -- test the REFPROP gas model by comparing to
data from the REFPROP GUI interface.

.. Author: Peter Blyton
.. Version: 21/05/2012
"""

import unittest
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

class REFPROPTestR134A_FLD(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_file("real gas REFPROP", ["R134A.FLD", "two phase"])
        self.gas = Gas_data(self.gmodel)

    def test_temp_from_rhoe(self):
        # Within 2 phase region, speed of sound undefined
        self.gas.rho = 66.4905623976
        self.gas.e[0] = 314.657214669e3
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 300.0, places=2)
        self.assertAlmostEqual(self.gas.p/1e6, 0.702820647167, places=3)
        self.assertAlmostEqual(self.gas.quality, 0.5, places=4)

    def test_rho_from_pT(self):
        # Supercritical
        self.gas.T[0] = 400.0
        self.gas.p = 4.2e6
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(self.gas.rho, 201.900034379, places=9)
        self.assertAlmostEqual(self.gas.e[0]/1e3, 452.009183568, places=9)
        self.assertAlmostEqual(self.gas.a, 135.781351458, places=9)
        self.assertEqual(self.gas.quality, 999) # 999 means supercritical
        # Test calculation of transport coefficients
        self.gmodel.eval_transport_coefficients(self.gas)
        self.assertAlmostEqual(self.gas.mu/1e6, 19.4290589130, places=9)
        self.assertAlmostEqual(self.gas.k[0]*1e3, 27.0975847050, places=9)

    def test_press_from_rhoT(self):
        # Saturated vapor
        self.gas.T[0] = 300.0
        self.gas.rho = 34.1928366481
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p/1e6, 0.702820647167, places=9)
        self.assertAlmostEqual(self.gas.e[0]/1e3, 392.711079900, places=9)
        self.assertAlmostEqual(self.gas.quality, 1.0, places=9)

    def test_other_functions(self):
        # Superheated vapor
        self.gas.T[0] = 400.0
        self.gas.rho = 3.08911589927 # corresponds to pressure of 0.1 MPa
        self.assertAlmostEqual(self.gmodel.Cv(self.gas)/1e3, 0.925309754316, places=9)
        self.assertAlmostEqual(self.gmodel.internal_energy(self.gas, 1)/1e3, 486.926162331, places=9)
        self.assertAlmostEqual(self.gmodel.enthalpy(self.gas, 1)/1e3, 519.297883970, places=9)
        self.assertAlmostEqual(self.gmodel.entropy(self.gas, 1)/1e3, 2.17396506192, places=9)

class REFPROPTestR134A_FLD_singlephase(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_file("real gas REFPROP", ["R134A.FLD", "single phase"])
        self.gas = Gas_data(self.gmodel)

    def test_rho_from_pT(self):
        # Supercritical
        self.gas.T[0] = 400.0
        self.gas.p = 4.2e6
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(self.gas.rho, 201.900034379, places=9)
        self.assertAlmostEqual(self.gas.e[0]/1e3, 452.009183568, places=9)
        self.assertAlmostEqual(self.gas.a, 135.781351458, places=9)
        self.assertEqual(self.gas.quality, 999) # 999 means supercritical
        # Test calculation of transport coefficients
        self.gmodel.eval_transport_coefficients(self.gas)
        self.assertAlmostEqual(self.gas.mu/1e6, 19.4290589130, places=9)
        self.assertAlmostEqual(self.gas.k[0]*1e3, 27.0975847050, places=9)
    
    def test_SS_from_rhoe(self):
        # Test out the fast, single phase sound speed calculation
        self.gas.rho = 4.17309524172
        self.gas.e[0] = 402.163750920e3
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.p/1e6, 0.1, places=9)
        self.assertAlmostEqual(self.gas.T[0], 300.0, places=9)
        self.assertAlmostEqual(self.gas.a, 162.07, places=2)

class REFPROPTestAIR_PPF(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_file("real gas REFPROP", ["AIR.PPF", "two phase"])
        self.gas = Gas_data(self.gmodel)

    def test_press_from_rhoT(self):
        self.gas.T[0] = 300.0
        self.gas.rho = 1.16159962683
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p/1e6, 0.1, places=9)
        self.assertAlmostEqual(self.gas.e[0]/1e3, 340.212604297, places=9)
        self.assertAlmostEqual(self.gas.a, 347.318504450, places=9)

class REFPROPTestAIR_MIX(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_file("real gas REFPROP", ["AIR.MIX", "two phase"])
        self.gas = Gas_data(self.gmodel)

    def test_press_from_rhoT(self):
        self.gas.T[0] = 300.0
        self.gas.rho = 1.16129872304
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p/1e6, 0.1, places=9)
        # Different default thermodynamic reference state used by mixture
        self.assertAlmostEqual(self.gas.e[0]/1e3, 214.201724633, places=9)
        self.assertAlmostEqual(self.gas.a, 347.369957553, places=9)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestR134A_FLD)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestR134A_FLD_singlephase)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestAIR_PPF)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestAIR_MIX)
    unittest.TextTestRunner(verbosity=2).run(suite)
