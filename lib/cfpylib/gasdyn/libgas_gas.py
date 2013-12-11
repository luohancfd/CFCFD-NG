#! /usr/bin/env python
"""
libgas_gas.py: access the gas models from the libgas library using the
cfpylib/gasdyn interface.

.. Author: Peter J Blyton
.. Version: 21/06/2012
.. Version: 11-Dec-2013 generalised a little by PeterJ
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
    def __init__(self, fname='gas-model.lua', massf=None, molef=None):
        """
        Set up the libgas model from the generic input file.

        :param fname: gas-model config file
        :param massf: optional dictionary of mass fractions
        :param molef: optional dictionary of mole fractions

        Rowan's thermochemistry module uses the Lua file to define
        the gas model, in detail.  There are so many options for 
        this input file that we whimp out and delegate the construction
        of a suitable file to other tools.  One such tool is gasmodel.py
        which, in turn, delegates all of it's work to Rowan's Lua
        program gasfile.lua.
        """
        if not libgas_ok:
            raise ImportError("Cannot use libgas_gas model because gaspy cannot be found.")
        self.fname = fname
        self.gasModel = create_gas_model(fname)
        self.gasData = Gas_data(self.gasModel)
        if massf is None and molef is None:
            name0 = self.gasModel.species_name(0)
            if name0 == "LUT": # [todo] we really need to fix the look-up-table code.
                set_massf(self.gasData, self.gasModel, [1.0,])
            else:
                set_molef(self.gasData, self.gasModel, {name0:1.0})
        elif (type(massf) is dict) or (type(massf) is list):
            set_massf(self.gasData, self.gasModel, massf)
        elif (type(molef) is dict) or (type(molef) is list):
            set_molef(self.gasData, self.gasModel, molef)
        self.set_pT(100.0e3, 300.0)
        return

    def clone(self):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        other = Gas(self.fname)
        nsp = self.gasModel.get_number_of_species()
        other.gasData.massf = self.gasData.massf
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
        self.gasData.T[0] = T # [todo] consider all modes
        # Calculate density, sound speed, internal energy and quality if available
        self.gasModel.eval_thermo_state_pT(self.gasData)
        self.rho = self.gasData.rho
        self.a = self.gasData.a
        self.e = self.gasData.e[0] # [todo] consider all modes
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
            self.k = self.gasData.k[0] # [todo] sum over all modes
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
                   % (self.p, self.T, self.rho, self.e, self.h, self.a) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.C_p, self.mu, self.k) )
        strm.write('    name: %s\n' % self.name)
        return

def make_gas_from_name(gasName):
    """
    Manufacture a Gas object from a small library of options.

    :param gasName: one of the names for the special cases set out below.
    """
    if gasName.lower() in ['co2-refprop']:
        os.system('gasmodel.py --model="real gas REFPROP"'+
                  ' --species="CO2.FLD" --output="co2-refprop.lua"')
        return Gas('co2-refprop.lua')
    elif gasName.lower() in ['co2-bender']:
        os.system('gasmodel.py --model="real gas Bender"'+
                  ' --species="CO2" --output="co2-bender.lua"')
        return Gas('co2-bender.lua')
    elif gasName.lower() in ['air-thermally-perfect']:
        os.system('gasmodel.py --model="thermally perfect gas" --species="N2 O2"')
        return Gas('gas-model.lua', molef={'O2':0.21, 'N2':0.79})
    elif gasName.lower() in ['r134a-refprop']:
        os.system('gasmodel.py --model="real gas REFPROP"'+
                  ' --species="R134A.FLD" --output="r134a-refprop.lua"')
        return Gas('r134a-refprop.lua')
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName
