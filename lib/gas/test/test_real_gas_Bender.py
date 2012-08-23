#!/usr/bin/env python
"""
test_real_gas_Bender.py -- test the Bender gas model for the avaliable fluids.

.. Author: Peter Blyton
.. Version: 22/03/2012
"""

import unittest
import sys
import os
import math
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from gaspy import *

class BenderTestCO2(unittest.TestCase):
    """Compare the CO2 implementation against the Bender generated thermodynamic
    table of 'Reynolds, WC (1979). Thermodynamic Properties in SI.'"""
    def setUp(self):
        self.gmodel = create_gas_file("real gas Bender", ["CO2"])
        self.gas = Gas_data(self.gmodel)
        set_molef(self.gas, self.gmodel, {'CO2':1})

    def test_temp_from_rhoe(self):
        # Saturated vapor
        self.gas.rho = 1.0/0.02140
        self.gas.e[0] = 358.59e3 - 1.788e6/self.gas.rho
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 250.0, places=1)
        self.assertAlmostEqual(self.gas.p/1e6, 1.788, places=3)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s/1e3, 1.45, places=3)
        # Comparisons to REFPROP to test other parts of the model
        self.assertAlmostEqual(self.gas.a, 221.22, delta=1.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 745.91, delta=25.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 1236.6, delta=40.0)

    def test_rho_from_pT(self):
        # Supercritical
        self.gas.T[0] = 400.0
        self.gas.p = 20.0e6
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(1.0/self.gas.rho, 0.00262, places=5)
        self.assertAlmostEqual((self.gas.e[0] + self.gas.p/self.gas.rho)/1e3, 403.03, places=1)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s/1e3, 1.2634, places=4)
        # Comparisons to REFPROP to test other parts of the model
        self.assertAlmostEqual(self.gas.a, 310.75, delta=2.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 888.27, delta=8.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 1886.8, delta=25.0)

    def test_press_from_rhoT(self):
        # Saturated vapor, near critical point
        self.gas.T[0] = 300.0
        self.gas.rho = 1.0/0.003749
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p/1e6, 6.706, places=3)
        self.assertAlmostEqual((self.gas.e[0] + self.gas.p/self.gas.rho)/1e3, 309.95, places=1)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s/1e3, 1.1120, places=4)
        # Comparisons to REFPROP to test other parts of the model
        self.assertAlmostEqual(self.gas.a, 185.33, delta=10.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 1247.6, delta=240.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 11921.0, delta=2200.0)

    def test_temp_from_rhop(self):
        # Superheated vapor
        self.gas.rho = 1.0/0.05379
        self.gas.p = 1.0e6
        self.gmodel.eval_thermo_state_rhop(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 300.0, places=1)
        self.assertAlmostEqual((self.gas.e[0] + self.gas.p/self.gas.rho)/1e3, 419.95, places=1)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s/1e3, 1.7737, places=4)
        # Comparisons to REFPROP to test other parts of the model
        self.assertAlmostEqual(self.gas.a, 262.42, delta=0.05)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 682.17, delta=2.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 920.89, delta=2.0)

class BenderTestH2O(unittest.TestCase):
    """Compare the H2O implementation against data from REFPROP.
    Polt and Maurer (1992) do not provide validation data for their Bender fit.
    Note that REFPROP uses the IAPWS Helmholtz formulation which is much more
    accurate, especially around the critical point."""
    def setUp(self):
        self.gmodel = create_gas_file("real gas Bender", ["H2O"])
        self.gas = Gas_data(self.gmodel)
        set_molef(self.gas, self.gmodel, {'H2O':1})

    def test_temp_from_rhoe(self):
        # Superheated vapor
        self.gas.rho = 3.1305
        self.gas.e[0] = 3002.3e3
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 700.0, delta=6.0)
        self.assertAlmostEqual(self.gas.p, 1.0e6, delta=4e4)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s, 7.5504e3, delta=10.0)
        self.assertAlmostEqual(self.gas.a, 640.55, delta=3.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 1644.7, delta=40.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 2136.8, delta=80.0)

    def test_rho_from_pT(self):
        # Superheated vapor
        self.gas.T[0] = 700.0
        self.gas.p = 1.0e6
        self.gmodel.eval_thermo_state_pT(self.gas)
        self.assertAlmostEqual(self.gas.rho, 3.1305, delta=0.1)
        self.assertAlmostEqual(self.gas.e[0], 3002.3e3, delta=1e4)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s, 7.5504e3, delta=10.0)
        self.assertAlmostEqual(self.gas.a, 640.55, delta=5.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 1644.7, delta=30.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 2136.8, delta=60.0)

    def test_press_from_rhoT(self):
        # Saturated vapor
        self.gas.T[0] = 300.0
        self.gas.rho = 0.025590
        self.gmodel.eval_thermo_state_rhoT(self.gas)
        self.assertAlmostEqual(self.gas.p, 3536.8, delta=10.0)
        self.assertAlmostEqual(self.gas.e[0], 2411.6e3, delta=60.0)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s, 8.5174e3, delta=1.0)
        self.assertAlmostEqual(self.gas.a, 427.89, delta=2.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 1442.2, delta=40.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 1914.1, delta=60.0)

    def test_temp_from_rhop(self):
        # Saturated vapor nearing critical point, note Tc = 647K
        self.gas.rho = 1.3694
        self.gas.p = 0.24577e6
        self.gmodel.eval_thermo_state_rhop(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 400.0, delta=15.0)
        self.assertAlmostEqual(self.gas.e[0], 2536.2e3, delta=1500.0)
        s = self.gmodel.mixture_entropy(self.gas)
        self.assertAlmostEqual(s, 7.0581e3, delta=11.0)
        self.assertAlmostEqual(self.gas.a, 484.67, delta=5.0)
        Cv = self.gmodel.Cv(self.gas)
        self.assertAlmostEqual(Cv, 1643.6, delta=220.0)
        Cp = self.gmodel.Cp(self.gas)
        self.assertAlmostEqual(Cp, 2218.3, delta=330.0)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(BenderTestCO2)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(BenderTestH2O)
    unittest.TextTestRunner(verbosity=2).run(suite)
