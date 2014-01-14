#! /usr/bin/env python
"""
ideal_gas.py: Thermodynamic properties of an ideal gas.

This module provides a look-alike Gas class for use in 
the gas flow functions.  Whereever cea2_gas works, so should this.

.. Author: 
   PA Jacobs
   School of Mechanical Engineering
   The University of Queensland

.. Versions:
   02-Apr-12: first cut from cea2_gas.py
"""

import sys, math

R_universal = 8314.0;  # J/kgmole.K

class Gas(object):
    """
    Provides the place to keep property data for the ideal gas.
    """
    def __init__(self, Mmass=28.96, gamma=1.4, name='air', 
                 s1=0.0, T1=298.15, p1=101.325e3,
                 mu_ref=1.716e-5, T_ref=273.0, S_mu=111.0,
                 Prandtl=0.71):
        """
        Set up a new object, from either a name of species list.

        :param Mmass: molecular mass, g/mole
        :param gamma: ratio of specific heats
        :param name: string name of gas (something like a species name in cea2_gas)
        :param s1: reference entropy, J/kg/K
        :param T1: temperature for reference entropy, K
        :param p1: pressure for reference entropy, Pa
        :param mu_ref: reference viscosity for Sutherland expression, Pa.s
        :param T_ref: reference temperature for Sutherland expression, degree K
        :param S_mu: constant (degree K) in Sutherlans expression
        :param Prandtl: mu.C_p/k
        """
        assert gamma > 1.0 and gamma <= 2.0, ('odd value: gamma=%g' % gamma)
        assert Mmass > 1.0 and Mmass < 1000.0, ('odd value: Mmass=%g' % Mmass)
        self.Mmass = Mmass
        self.R = R_universal / Mmass
        self.gam = gamma
        self.C_v = self.R / (gamma - 1)
        self.C_p = self.R + self.C_v
        self.name = name    
        # reference entropy
        self.s1 = s1
        self.T1 = T1
        self.p1 = p1
        # Data for transport properties, based on Sutherland variation.
        self.mu_ref = mu_ref
        self.T_ref = T_ref
        self.S_mu = S_mu
        self.Prandtl = Prandtl
        # set default thermo conditions
        self.set_pT(100.0e3, 300.0)
        return

    def clone(self):
        """
        Clone the current Gas object to make another, just the same.

        :returns: the new Gas object.
        """
        other = Gas(self.Mmass, self.gam, self.name, 
                    s1=self.s1, T1=self.T1, p1=self.p1,
                    mu_ref=self.mu_ref, T_ref=self.T_ref, S_mu=self.S_mu,
                    Prandtl=self.Prandtl)
        other.set_pT(self.p, self.T)
        return other

    def set_pT(self, p, T, transProps=True):
        """
        Fills out gas state from given pressure and temperature.

        :param p: pressure, Pa
        :param T: temperature, K
        :param transProps: if True, compute transport properties as well.
        """
        self.p = p
        self.T = T
        self.rho = p / (self.R * T)
        self.a = math.sqrt(self.gam * self.R * T)
        self.e = self.C_v * T
        self.h = self.C_p * T
        self.s = self.s1 + self.C_p * math.log(T/self.T1) - self.R * math.log(p/self.p1)
        if transProps:
            self.mu = self.mu_ref * (T/self.T_ref)**1.5 * (self.T_ref+self.S_mu)/(T+self.S_mu)
            self.k = self.mu * self.C_p / self.Prandtl
        else:
            self.mu = 0.0
            self.k = 0.0
        return

    def set_rhoT(self, rho, T, transProps=True):
        """
        Fills out gas state from given density and temperature.

        :param rho: density, kg/m**3
        :param T: temperature, K
        """
        p = rho * self.R * T
        return self.set_pT(p, T, transProps)

    def set_ps(self, p, s, transProps=True):
        """
        Fills out gas state from given pressure and specific entropy.

        :param p: pressure, Pa
        :param s: entropy, J/(kg.K)
        """
        cp_ln_TT1 = s - self.s1 + self.R * math.log(p/self.p1)
        T = self.T1 * math.exp(cp_ln_TT1 / self.C_p)
        return self.set_pT(p, T, transProps)

    def set_ph(self, p, h, transProps=True):
        """
        Fills out gas state from given pressure and enthalpy.

        :param p: pressure, Pa
        :param h: enthalpy, J/kg
        """
        T = h / self.C_p
        return self.set_pT(p, T, transProps)

    def write_state(self, strm):
        """
        Writes the gas state data to the specified stream.
        """
        strm.write('    p: %g Pa, T: %g K, rho: %g kg/m**3, e: %g J/kg, h: %g J/kg, a: %g m/s, s: %g J/(kg.K)\n'
                   % (self.p, self.T, self.rho, self.e, self.h, self.a, self.s) )
        strm.write('    R: %g J/(kg.K), gam: %g, Cp: %g J/(kg.K), mu: %g Pa.s, k: %g W/(m.K)\n'
                   % (self.R, self.gam, self.C_p, self.mu, self.k) )
        strm.write('    name: %s\n' % self.name)
        return

def make_gas_from_name(gasName):
    """
    Manufacture a Gas object from a small library of options.

    :param gasName: one of the names for the special cases set out below
    """
    if gasName in ['air', 'Air', 'air5species']:
        return Gas()
    elif gasName in ['n2', 'N2', 'nitrogen']:
        return Gas(Mmass=28.0, gamma=1.4, name='N2', 
                   s_1=0.0, T_1=298.15, p_1=101.325e3,
                   mu_ref=1.663e-5, T_ref=273.0, S_mu=107.0,
                   Prandtl=0.71)
    elif gasName in ['co2', 'CO2', 'carbon dioxide', 'carbon-dioxide']:
        return Gas(Mmass=44.0, gamma=1.301, name='CO2', 
                   s_1=0.0, T_1=298.15, p_1=101.325e3,
                   mu_ref=1.370e-5, T_ref=273.0, S_mu=222.0,
                   Prandtl=0.72)
    else:
        raise Exception, 'make_gas_from_name(): unknown gasName: %s' % gasName

def list_gas_names():
    """
    :returns: the list of gases available in make_gas_from_name()
    """
    return ['air', 'n2', 'co2']

# --------------------------------------------------------------

if __name__ == '__main__':
    print 'Test/demonstrate the Gas class...'
    print 'gases available in make_gas_from_name():'
    for name in list_gas_names():
        print "   ", name
    #
    print '\nDefault constructor with Air as the test gas.'
    a = Gas()
    a.set_pT(100.0e3, 300.0)
    a.write_state(sys.stdout)
    print 'and the same Air at a higher temperature'
    a.set_pT(100.0e3, 4000.0)
    a.write_state(sys.stdout)
    #
    print '\nCheck enthalpy specification'
    b = make_gas_from_name('air')
    b.set_ph(a.p, a.h)
    b.write_state(sys.stdout)
    #
    print '\nCheck entropy specification'
    b = make_gas_from_name('air')
    b.set_ps(a.p, a.s)
    b.write_state(sys.stdout)
    #
    print 'End of test.'

