#!/usr/bin/env python
"""
test_gas_flow_with_ideal_gas.py -- test script

Demonstrate that the gas_flow functions can work for ideal_gas 
as well as cea2_gas.

PJ, 02-Apr-2012
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

import unittest
from cfpylib.gasdyn.ideal_gas import Gas
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import p2_p1, T2_T1, u2_u1, theta_obl


class TestGasFlowIdeal(unittest.TestCase):
    def test_shock_given_V(self):
        s1 = Gas()
        s1.set_pT(1.0e5, 300.0)
        self.assertAlmostEqual(s1.rho, 1.1612, delta=0.001)
        s2 = s1.clone()
        V1 = 3000.0
        M1 = V1/s1.a
        V2,Vg = normal_shock(s1, V1, s2)
        p2_ideal = s1.p * p2_p1(M1)
        T2_ideal = s1.T * T2_T1(M1)
        V2_ideal = V1 * u2_u1(M1)
        Vg_ideal = V1 - V2_ideal
        self.assertAlmostEqual(s2.p, p2_ideal, delta=1.0)
        self.assertAlmostEqual(s2.T, T2_ideal, delta=1.0)
        self.assertAlmostEqual(V2, V2_ideal, delta=1.0)
        self.assertAlmostEqual(Vg, Vg_ideal, delta=1.0)
        return

    def test_shock_given_p_ratio(self):
        s1 = Gas()
        s1.set_pT(1.0e5, 300.0)
        self.assertAlmostEqual(s1.rho, 1.1612, delta=0.001)
        V1, V2, Vg, s2 = normal_shock_p2p1(s1, p2_p1(3000/s1.a))
        self.assertAlmostEqual(V1, 3000.0, delta=1.0)
        return

    def test_finite_wave(self):
        V1 = 0.0
        s9 = Gas()
        s9.set_pT(1.0e5, 320.0) # keep gas close to ideal
        Jplus = V1 + 2*s9.a/(1.4-1)
        V2, s10 = finite_wave_dp('cplus', V1, s9, 60.0e3)
        ideal_V2 = Jplus - 2*s10.a/(1.4-1)
        self.assertAlmostEqual(V2, ideal_V2, delta=2.0)
        V1 = 0.0
        s9.set_pT(1.0e5, 320.0)
        Jplus = V1 + 2*s9.a/(1.4-1)
        V2, s10 = finite_wave_dv('cplus', V1, s9, 125.0)
        self.assertAlmostEqual(Jplus, V2 + 2*s10.a/(1.4-1), delta=2.0)
        return

    def test_oblique_shock(self):
        # essentially ideal gas at these conditions
        s1 = Gas()
        s1.set_pT(100.0e3, 300.0)
        M1 = 1.5
        beta = 45.0 * math.pi/180
        V1 = M1 * s1.a
        theta, V2, s2 = theta_oblique(s1, V1, beta)
        self.assertAlmostEqual(s2.p, 114580, delta=50)
        self.assertAlmostEqual(theta, theta_obl(M1, beta), delta=0.001)
        beta2 = beta_oblique(s1, V1, theta)
        self.assertAlmostEqual(beta2*180/math.pi, 45.0, delta=0.01)
        return

    def test_taylor_maccoll_given_beta(self):
        # At these conditions, expect essentially-ideal gas flow.
        # M1=1.5, beta=49deg, expect theta=20deg from NACA1135.
        s1 = Gas()
        s1.set_pT(100.0e3, 300.0)
        M1 = 1.5
        V1 = M1 * s1.a
        beta = 49.0 * math.pi/180
        theta_c, V_c, s_c = theta_cone(s1, V1, beta)
        self.assertAlmostEqual(theta_c*180.0/math.pi, 20.0, delta=0.05)
        surface_pressure_coeff = (s_c.p - s1.p)/(0.5*s1.rho*V1*V1) 
        self.assertAlmostEqual(surface_pressure_coeff, 0.385, delta=0.01)
        # M1=1.8, beta=45deg, theta=24deg from NACA1135.
        M1 = 1.8
        V1 = M1 * s1.a
        beta = 45.0 * math.pi/180
        theta_c, V_c, s_c = theta_cone(s1, V1, beta)
        self.assertAlmostEqual(theta_c*180.0/math.pi, 24.0, delta=0.1)
        surface_pressure_coeff = (s_c.p - s1.p)/(0.5*s1.rho*V1*V1) 
        self.assertAlmostEqual(surface_pressure_coeff, 0.466, delta=0.01)
        return

    def test_taylor_maccoll_given_theta(self):
        # Conical shock from cone with half-angle 20deg.
        s1 = Gas()
        s1.set_pT(100.0e3, 300.0)
        M1 = 1.5
        V1 = M1 * s1.a
        beta = beta_cone(s1, V1, 20.0*math.pi/180)
        self.assertAlmostEqual(beta*180/math.pi, 49.0, delta=0.05)
        return

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGasFlowIdeal)
    unittest.TextTestRunner(verbosity=2).run(suite)
