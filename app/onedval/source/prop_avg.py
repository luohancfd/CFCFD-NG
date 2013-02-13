# Author: Rowan J. Gollan
# Place: UQ, Brisbane, Queensland, Australia
# Date: 26-Jun-2012
#
# This module provides methods to give
# one-dimensionalised (averaged) flow properties
# based on a variety of techniques.
#

import sys

from math import sqrt
from libprep3 import Vector3, dot
from libprep3 import Gas_data, set_massf
from cfpylib.nm.zero_solvers import secant
from cfpylib.nm.nelmin import minimize

def compute_fluxes(cells, var_map):
    f_mass = 0.0
    f_mom = Vector3(0.0, 0.0, 0.0)
    f_energy = 0.0
    N = Vector3(0.0, 0.0, 0.0)
    A = 0.0
    rholabel = var_map['rho']
    plabel = var_map['p']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    h0label = var_map['h0']
    for c in cells:
        dA = c.area()*R0*R0
        n = c.normal()
        rho = c.get(rholabel)
        p = c.get(plabel)
        vel = Vector3(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        h0 = c.get(h0label)
        # Add increments
        f_mass = f_mass + rho*dot(vel, n)*dA
        f_mom = f_mom + (rho*dot(vel, n)*vel + p*n)*dA
        f_energy = f_energy + (rho*dot(vel, n)*h0)*dA
    return {'mass': f_mass, 'mom': f_mom, 'energy': f_energy}

def area_weighted_avg(cells, props, var_map):
    phis = dict.fromkeys(props, 0.0)
    area = 0.0
    for c in cells:
        dA = c.area()
        area = area + dA
        for p in props:
            label = p
            if p in var_map:
                label = var_map[p]
            phis[p] = phis[p] + c.get(label)*dA
    
    for p in props:
        phis[p] = phis[p]/area
    
    return phis

def mass_flux_weighted_avg(cells, props, var_map):
    phis = dict.fromkeys(props, 0.0)
    f_mass = 0.0
    rholabel = var_map['rho']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    for c in cells:
        dA = c.area()
        rho = c.get(rholabel)
        vel = Vector3(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        n = c.normal()
        w = rho*dot(vel, n)
        f_mass = f_mass + w*dA
        for p in props:
            label = p
            if p in var_map:
                label = var_map[p]
            phis[p] = phis[p] + c.get(label)*w*dA
    
    for p in props:
        phis[p] = phis[p]/f_mass

    return phis


def stream_thrust_avg(cells, props, var_map, species, gmodel):
    f_mass = 0.0
    f_mom = Vector3(0.0, 0.0, 0.0)
    f_energy = 0.0
    f_sp = [0.0,]*len(species)
    nsp = gmodel.get_number_of_species()
    N = Vector3(0.0, 0.0, 0.0)
    A = 0.0
    rholabel = var_map['rho']
    plabel = var_map['p']
    ulabel = var_map['u']
    vlabel = var_map['v']
    wlabel = var_map['w']
    h0label = var_map['h0']
    for c in cells:
        dA = c.area() 
        n = c.normal()
        rho = c.get(rholabel)
        p = c.get(plabel)
        vel = Vector3(c.get(ulabel),
                      c.get(vlabel),
                      c.get(wlabel))
        h0 = c.get(h0label)
        # Add increments
        u_n = dot(vel, n)
        f_mass = f_mass + rho*u_n*dA
        f_mom = f_mom + (rho*u_n*vel + p*n)*dA
        f_energy = f_energy + (rho*u_n*h0)*dA
        if nsp > 1:
            for isp, sp in enumerate(species):
                massf = c.get(sp)
                f_sp[isp] = f_sp[isp] + rho*massf*u_n*dA
        A = A + dA
        N = N + n*dA

    N = N / A
    f_mom_s = dot(f_mom, N)
    Q = Gas_data(gmodel)
    if nsp > 1:
        massf = f_sp/f_mass
        set_massf(Q, gmodel, massf)
    else:
        Q.massf[0] = 1.0
    

    def f_to_minimize(x):
        rho, T, u = x
        # Use equation of state to compute other thermo quantities
        Q.rho = rho
        Q.T[0] = T
        gmodel.eval_thermo_state_rhoT(Q)
        h = gmodel.mixture_enthalpy(Q)
        p = Q.p
        # Compute errors
        fmass_err = abs(f_mass - rho*u*A)/(abs(f_mass) + 1.0)
        fmom_err = abs(f_mom_s - (rho*u*u + p)*A)/(abs(f_mom_s) + 1.0)
        fe_err = abs(f_energy - (rho*u*A*(h + 0.5*u*u)))/abs(f_energy + 1.0)
        # Total error is the sum
        return fmass_err + fmom_err + fe_err

    # Compute an initial guess based on mass-flux weighted averages
    mfw_props = mass_flux_weighted_avg(cells, ['rho', 'T', 'u', 'v', 'w'], var_map)
    print "mfw_props= ", mfw_props
    u = sqrt(mfw_props['u']**2 + mfw_props['v']**2 + mfw_props['w']**2)
    guess = [mfw_props['rho'], mfw_props['T'], u]
    x, fx, conv_flag, nfe, nres = minimize(f_to_minimize, guess, [0.01, 10.0, 10.0])
    rho, T, u = x
    Q.rho = rho; Q.T[0] = T
    gmodel.eval_thermo_state_rhoT(Q)
    Q.print_values()
    sys.exit(1)
    return
    
    
    

    
    
        

