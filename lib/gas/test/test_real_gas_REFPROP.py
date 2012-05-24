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

class REFPROPTestR134A_FLD(unittest.TestCase):
    def setUp(self):
        self.gmodel = create_gas_file("real gas REFPROP", ["R134A.FLD", "two phase"])
        self.gas = Gas_data(self.gmodel)

    def test_temp_from_rhoe(self):
        # Within 2 phase region
        self.gas.rho = 66.491
        self.gas.e[0] = 314.66e3
        self.gmodel.eval_thermo_state_rhoe(self.gas)
        self.assertAlmostEqual(self.gas.T[0], 300.0, places=2)
        self.assertAlmostEqual(self.gas.p, 300.0, places=2)

#    def test_rho_from_pT(self):
#        # Supercritical
#        self.gas.T[0] = 400.0
#        self.gas.p = 20.0e6
#        self.gmodel.eval_thermo_state_pT(self.gas)
#        self.assertAlmostEqual(self.gas.rho, 380.58, delta=0.02)
#        self.assertAlmostEqual((self.gas.e[0] + e_offset)/1e3, 432.24, delta=0.2)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s, 1.7875e3 - s_offset, delta=1.0)
#        #self.assertAlmostEqual(self.gas.a, 310.33, delta=2.0)
#        Cv = self.gmodel.Cv(self.gas)
#        self.assertAlmostEqual(Cv, 881.90, delta=2.0)
#        #Cp = self.gmodel.Cp(self.gas) # Something wrong with calculation of Cp
#        #self.assertAlmostEqual(Cp, 1878.4, delta=5.0)

#    def test_press_from_rhoT(self):
#        # Saturated vapor, near critical point
#        self.gas.T[0] = 300.0
#        self.gas.rho = 264.89
#        self.gmodel.eval_thermo_state_rhoT(self.gas)
#        self.assertAlmostEqual(self.gas.p/1e3, 6703.7, delta=0.5)
#        self.assertAlmostEqual((self.gas.e[0] + e_offset)/1e3, 364.39, delta=0.1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s, 1.6303e3 - s_offset, delta=1.0)
#        #self.assertAlmostEqual(self.gas.a, 191.02, places=2)
#        Cv = self.gmodel.Cv(self.gas)
#        self.assertAlmostEqual(Cv, 1094.1, delta=1.0)
#        #Cp = self.gmodel.Cp(self.gas)
#        #self.assertAlmostEqual(Cp, 1006.3, delta=5.0)

#    def test_temp_from_rhop(self):
#        # Superheated vapor
#        self.gas.rho = 18.576
#        self.gas.p = 1.0e6
#        self.gmodel.eval_thermo_state_rhop(self.gas)
#        self.assertAlmostEqual(self.gas.T[0], 300.0, delta=0.01)
#        self.assertAlmostEqual((self.gas.e[0] + e_offset)/1e3, 445.85, delta=0.1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s, 2.2923e3 - s_offset, delta=1.0)
#        self.assertAlmostEqual(self.gas.a, 262.37, delta=0.05)
#        Cv = self.gmodel.Cv(self.gas)
#        self.assertAlmostEqual(Cv, 684.66, delta=1.0)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp, 923.55, delta=1.0)

#class REFPROPTestAIR_PPF(unittest.TestCase):
#    def setUp(self):
#        self.gmodel = create_gas_file("real gas REFPROP", ["AIR.PPF", "single phase"])
#        self.gas = Gas_data(self.gmodel)

#    def test_temp_from_rhoe(self):
#        # Triple point properties of vapor
#        self.gas.rho = 0.028197 # Data of McLinden sensitive to more than 3dp.
#        self.gas.e[0] = 321.25e3
#        self.gmodel.eval_thermo_state_rhoe(self.gas)
#        self.assertAlmostEqual(self.gas.T[0], 273.15 - 103.30, places=1)
#        self.assertAlmostEqual(self.gas.p, 390.0, places=1)
#        self.assertAlmostEqual(self.gas.a, 127.0, places=0)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.9638, places=3)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp/1e3, 0.585, places=3)

#    def test_rho_from_pT(self):
#        # Saturated vapor
#        self.gas.T[0] = 273.15
#        self.gas.p = 292.69e3
#        self.gmodel.eval_thermo_state_pT(self.gas)
#        self.assertAlmostEqual(self.gas.rho, 14.420, places=3)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 378.38, places=2)
#        self.assertAlmostEqual(self.gas.a, 147.0, places=1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.7274, places=3)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp/1e3, 0.883, places=3)

#    def test_press_from_rhoT(self):
#        # Saturated vapor. At this point, Cv is right but Cp wrong, also causing
#        # SS to be wrong. See Span (2000) Multiparameter equations of state, pg 43.
#        # Cp of a saturation state is that in the limit when T -> Ts, as Cp
#        # is meaningless at the saturation boundary as phase change is a constant
#        # temperature process under isobaric conditions.
#        self.gas.rho = 267.60
#        self.gas.T[0] = 273.15 + 95.0
#        self.gmodel.eval_thermo_state_rhoT(self.gas)
#        self.assertAlmostEqual(self.gas.p/1e6, 3.5916, places=3)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 407.17, places=2)
#        #self.assertAlmostEqual(self.gas.a, 102.0, places=1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.6490, places=3)
#        #Cp = self.gmodel.Cp(self.gas)
#        #self.assertAlmostEqual(Cp, 4942.0, places=1)

#    def test_temp_from_rhop(self):
#        # Superheated vapor, comparison to REFPROP as McLinden only provide
#        # saturation tables.
#        self.gas.rho = 10.0
#        self.gas.p = 318.76e3
#        self.gmodel.eval_thermo_state_rhop(self.gas)
#        self.assertAlmostEqual(self.gas.T[0], 400.0, delta=0.02)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 485.59, delta=0.1)
#        self.assertAlmostEqual(self.gas.a, 185.03, delta=0.03)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s, 2.0762e3, delta=1.0)
#        Cv = self.gmodel.Cv(self.gas)
#        self.assertAlmostEqual(Cv, 929.94, delta=6.0)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp, 1021.6, delta=6.0)

#class REFPROPTestAIR_MIX(unittest.TestCase):
#    def setUp(self):
#        self.gmodel = create_gas_file("real gas REFPROP", ["AIR.MIX", "single phase"])
#        self.gas = Gas_data(self.gmodel)

#    def test_temp_from_rhoe(self):
#        # Triple point properties of vapor
#        self.gas.rho = 0.028197 # Data of McLinden sensitive to more than 3dp.
#        self.gas.e[0] = 321.25e3
#        self.gmodel.eval_thermo_state_rhoe(self.gas)
#        self.assertAlmostEqual(self.gas.T[0], 273.15 - 103.30, places=1)
#        self.assertAlmostEqual(self.gas.p, 390.0, places=1)
#        self.assertAlmostEqual(self.gas.a, 127.0, places=0)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.9638, places=3)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp/1e3, 0.585, places=3)

#    def test_rho_from_pT(self):
#        # Saturated vapor
#        self.gas.T[0] = 273.15
#        self.gas.p = 292.69e3
#        self.gmodel.eval_thermo_state_pT(self.gas)
#        self.assertAlmostEqual(self.gas.rho, 14.420, places=3)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 378.38, places=2)
#        self.assertAlmostEqual(self.gas.a, 147.0, places=1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.7274, places=3)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp/1e3, 0.883, places=3)

#    def test_press_from_rhoT(self):
#        # Saturated vapor. At this point, Cv is right but Cp wrong, also causing
#        # SS to be wrong. See Span (2000) Multiparameter equations of state, pg 43.
#        # Cp of a saturation state is that in the limit when T -> Ts, as Cp
#        # is meaningless at the saturation boundary as phase change is a constant
#        # temperature process under isobaric conditions.
#        self.gas.rho = 267.60
#        self.gas.T[0] = 273.15 + 95.0
#        self.gmodel.eval_thermo_state_rhoT(self.gas)
#        self.assertAlmostEqual(self.gas.p/1e6, 3.5916, places=3)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 407.17, places=2)
#        #self.assertAlmostEqual(self.gas.a, 102.0, places=1)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s/1e3, 1.6490, places=3)
#        #Cp = self.gmodel.Cp(self.gas)
#        #self.assertAlmostEqual(Cp, 4942.0, places=1)

#    def test_temp_from_rhop(self):
#        # Superheated vapor, comparison to REFPROP as McLinden only provide
#        # saturation tables.
#        self.gas.rho = 10.0
#        self.gas.p = 318.76e3
#        self.gmodel.eval_thermo_state_rhop(self.gas)
#        self.assertAlmostEqual(self.gas.T[0], 400.0, delta=0.02)
#        self.assertAlmostEqual(self.gas.e[0]/1e3, 485.59, delta=0.1)
#        self.assertAlmostEqual(self.gas.a, 185.03, delta=0.03)
#        s = self.gmodel.total_entropy(self.gas)
#        self.assertAlmostEqual(s, 2.0762e3, delta=1.0)
#        Cv = self.gmodel.Cv(self.gas)
#        self.assertAlmostEqual(Cv, 929.94, delta=6.0)
#        Cp = self.gmodel.Cp(self.gas)
#        self.assertAlmostEqual(Cp, 1021.6, delta=6.0)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestR134A_FLD)
    unittest.TextTestRunner(verbosity=2).run(suite)

#    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestAIR_PPF)
#    unittest.TextTestRunner(verbosity=2).run(suite)

#    suite = unittest.TestLoader().loadTestsFromTestCase(REFPROPTestAIR_MIX)
#    unittest.TextTestRunner(verbosity=2).run(suite)
