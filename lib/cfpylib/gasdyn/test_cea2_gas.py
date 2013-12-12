#!/usr/bin/env python
"""
test_cea2_gas.py -- test script
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

import unittest
from cfpylib.gasdyn.cea2_gas import Gas, make_gas_from_name


class TestLibGasGas(unittest.TestCase):

    def test_low_temperature_air(self):
        testGas = make_gas_from_name('air')
        p = 100.0e3 # Pa
        T = 300.0 # degrees K
        testGas.set_pT(p, T)
        gam = testGas.gam
        R = testGas.R
        Cp = testGas.C_p
        self.assertAlmostEqual(gam, 1.4, delta=0.02)
        self.assertAlmostEqual(R, 287.0, delta=2)
        rho = testGas.rho
        self.assertAlmostEqual(rho, p/(R*T), delta=0.01)
        import math
        self.assertAlmostEqual(testGas.a, math.sqrt(gam*R*T), delta=2.0)
        s = testGas.s
        p2 = 150.0e3
        T2 = 350.0
        testGas.set_pT(p2, T2)
        s2 = testGas.s
        delta_s = Cp*math.log(T2/T) - R*math.log(p2/p)
        self.assertAlmostEqual(s2-s, delta_s, delta=0.3)
        testGas.set_ps(p2, s)
        rho2 = testGas.rho
        self.assertAlmostEqual(p2/pow(rho2,gam), p/pow(rho,gam), delta=100.0)
       
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestLibGasGas)
    unittest.TextTestRunner(verbosity=2).run(suite)
