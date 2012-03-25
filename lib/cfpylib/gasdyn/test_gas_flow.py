#!/usr/bin/env python
"""
test_gas_flow.py -- test script
"""

import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))

import unittest
from cfpylib.gasdyn.cea2_gas import Gas
from cfpylib.gasdyn.gas_flow import *
from cfpylib.gasdyn.ideal_gas_flow import theta_obl


class TestGasFlow(unittest.TestCase):
    def test_shock_given_V(self):
        s1 = Gas({'Air':1.0})
        s1.set_pT(1.0e5, 300.0)
        self.assertAlmostEqual(s1.rho, 1.1612, delta=0.001)
        s2 = s1.clone()
        V2,Vg = normal_shock(s1, 3000.0, s2)
        self.assertAlmostEqual(s2.p, 9.1779e+06, delta=1.0)
        self.assertAlmostEqual(s2.T, 3572.97, delta=1.0)
        self.assertAlmostEqual(V2, 394.09, delta=1.0)
        self.assertAlmostEqual(Vg, 2605.9, delta=1.0)
        return

    def test_shock_given_p_ratio(self):
        s1 = Gas({'Air':1.0})
        s1.set_pT(1.0e5, 300.0)
        self.assertAlmostEqual(s1.rho, 1.1612, delta=0.001)
        s2 = s1.clone()
        V1, V2, Vg, s2 = normal_shock_p2p1(s1, 9.1779e+06/s1.p)
        V2,Vg = normal_shock(s1, 3000.0, s2)
        self.assertAlmostEqual(V1, 3000.0, delta=1.0)
        self.assertAlmostEqual(s2.T, 3572.97, delta=1.0)
        self.assertAlmostEqual(V2, 394.09, delta=1.0)
        self.assertAlmostEqual(Vg, 2605.9, delta=1.0)
        return

    def test_reflected_shock(self):
        s1 = Gas({'Air':1.0})
        s1.set_pT(1.0e5, 300.0)
        s2 = s1.clone()
        V2,Vg = normal_shock(s1, 3000.0, s2)
        s5 = s1.clone()
        Vr_b = reflected_shock(s2, Vg, s5)
        self.assertAlmostEqual(s5.p, 8.4329e+07, delta=10.0)
        self.assertAlmostEqual(s5.T, 6121.63, delta=1.0)
        self.assertAlmostEqual(s5.rho, 43.909, delta=0.01)
        self.assertAlmostEqual(Vr_b, 656.69, delta=0.1)
        return

    def test_expand_from_stagnation(self):
        s1 = Gas({'Air':1.0})
        s1.set_pT(1.0e5, 300.0)
        s2 = s1.clone()
        V2,Vg = normal_shock(s1, 3000.0, s2)
        s5 = s1.clone()
        Vr_b = reflected_shock(s2, Vg, s5)
        s6, V = expand_from_stagnation(0.0025, s5)
        self.assertAlmostEqual(s6.p, 210820.0, delta=10.0)
        self.assertAlmostEqual(s6.T, 2331.29, delta=1.0)
        self.assertAlmostEqual(s6.rho, 0.31475, delta=0.0001)
        self.assertAlmostEqual(V, 3766.2, delta=0.1)
        s7 = total_condition(s6, V)
        self.assertAlmostEqual(s7.p, 8.4013e+07, delta=10.0)
        self.assertAlmostEqual(s7.T, 6120.72, delta=1.0)
        self.assertAlmostEqual(s7.rho, 43.748, delta=0.01)
        s8 = pitot_condition(s6, V)
        self.assertAlmostEqual(s8.p, 4.3691e+06, delta=10.0)
        self.assertAlmostEqual(s8.T, 5450.51, delta=1.0)
        self.assertAlmostEqual(s8.rho, 2.4222, delta=0.01)
        return

    def test_finite_wave(self):
        V1 = 0.0
        s9 = Gas({'Air':1.0})
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
        s1 = Gas({'Air':1.0})
        s1.set_pT(100.0e3, 300.0)
        M1 = 1.5
        beta = 45.0 * math.pi/180
        V1 = M1 * s1.a
        theta, V2, s2 = theta_oblique(s1, V1, beta)
        self.assertAlmostEqual(s2.p, 114620, delta=10)
        self.assertAlmostEqual(theta, theta_obl(M1, beta), delta=0.001)
        beta2 = beta_oblique(s1, V1, theta)
        self.assertAlmostEqual(beta2*180/math.pi, 45.0, delta=0.01)
        return

    def test_taylor_maccoll_given_beta(self):
        # At these conditions, expect essentially-ideal gas flow.
        # M1=1.5, beta=49deg, expect theta=20deg from NACA1135.
        s1 = Gas({'Air':1.0})
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
        s1 = Gas({'Air':1.0})
        s1.set_pT(100.0e3, 300.0)
        M1 = 1.5
        V1 = M1 * s1.a
        beta = beta_cone(s1, V1, 20.0*math.pi/180)
        self.assertAlmostEqual(beta*180/math.pi, 49.0, delta=0.05)
        return

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestGasFlow)
    unittest.TextTestRunner(verbosity=2).run(suite)
