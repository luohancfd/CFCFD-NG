#! /usr/bin/env python
"""
libgas_gas.py: access the gas models from the libgas library using the
cfpylib/gasdyn interface.

.. Author: Peter J Blyton
.. Version: 21/06/2012
"""

from ..nm.zero_solvers import secant
import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
try:
    from gaspy import *
    libgas_ok = True
except:
    libgas_ok = False

class Gas(object):
    """
    Provides the place to hold the libgas gas data object and gas model object.
    """
    def __init__(self, name='CO2.FLD', gasModelType='real gas REFPROP'):
        """
        Set up the libgas model from the specified gas model type.

        :param species: species name for the particular gas model type
        :param gasModelType: the particular libgas gas model type
        """
        if not libgas_ok:
            raise ImportError("Cannot use this gas model as libgas cannot be found!")
        self.name = name
        self.gasModelType = gasModelType
        self.gasModel = create_gas_file(gasModelType, [name])
        self.gasData = Gas_data(self.gasModel)
        set_molef(self.gasData, self.gasModel, {name:1})
        self.set_pT(100.0e3, 300.0)
        return

    def clone(self):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        other = Gas(self.name, self.gasModelType)
        other.set_pT(self.p, self.T)
        return other

    def set_pT(self, p, T, transProps=True):
        """
        Compute the thermodynamic state from given pressure and temperature.

        :param p: pressure, Pa
        :param T: temperature, K
        :param transProps: if True, compute transport properties as well.
        """
        self.p = p
        self.gasData.p = p
        self.T = T
        self.gasData.T[0] = T
        # Calculate density, sound speed, internal energy and quality if available
        self.gasModel.eval_thermo_state_pT(self.gasData)
        self.rho = self.gasData.rho
        self.a = self.gasData.a
        self.son = self.a
        self.u = self.gasData.e[0]
        self.e = self.u
        self.quality = self.gasData.quality
        # Manually call methods to calculate other thermodynamic properties
        self.h = self.gasModel.enthalpy(self.gasData, 0)
        self.s = self.gasModel.entropy(self.gasData, 0)
        self.R = self.gasModel.R(self.gasData)
        self.C_p = self.gasModel.Cp(self.gasData)
        self.C_v = self.gasModel.Cv(self.gasData)
        self.gam = self.gasModel.gamma(self.gasData)
        if transProps:
            self.gasModel.eval_transport_coefficients(self.gasData)
            self.mu = self.gasData.mu
            self.k = self.gasData.k[0]
        else:
            self.mu = 0.0
            self.k = 0.0
        return

    def set_rhoT(self, rho, T, transProps=True):
        """
        Compute the thermodynamic state from given density and temperature.

        :param rho: density, kg/m**3
        :param T: temperature, K
        :param transProps: if True, compute transport properties as well.
        """
        self.gasData.rho = rho
        self.gasData.T[0] = T
        self.gasModel.eval_thermo_state_rhoT(self.gasData)
        return self.set_pT(self.gasData.p, T, transProps)

    def set_ps(self, p, s, transProps=True):
        """
        Compute the thermodynamic state from given pressure and entropy

        :param p: pressure, Pa
        :param s: entropy, J/(kg.K)
        :param transProps: if True, compute transport properties as well.
        """
        # The libgas library does not have a pressure-entropy thermodynamic
        # state solver, so we need to do the iterative calculation ourselves.
        def entropy_solve(temp):
            self.set_pT(p, temp) # calculate density
            entropy = self.gasModel.entropy(self.gasData, 0) # entropy from temp and density
            return s - entropy
        T = secant(entropy_solve, 250.0, 260.0, tol=1.0e-4)
        if T == "FAIL": raise Exception("Secant solver failed, bailing out!")
        return self.set_pT(p, T, transProps)

    def write_state(self, strm):
        """
        Writes the gas state data to the specified stream.
        """
        strm.write('    p: %g Pa, T: %g K, rho: %g kg/m**3, e: %g J/kg, h: %g J/kg, a: %g m/s\n'
                   % (self.p, self.T, self.rho, self.u, self.h, self.son) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.C_p, self.mu, self.k) )
        strm.write('    name: %s\n' % self.name)
        return
