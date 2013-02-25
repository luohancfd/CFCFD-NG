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
from scipy.optimize import minimize

def area(cells):
    A = 0.0
    for c in cells:
        A = A + c.area()
    return A

def avg_pos(cells, var_map):
    x = 0.0; y = 0.0; z = 0.0
    xlabel = var_map['x']
    ylabel = var_map['y']
    zlabel = var_map['z']
    A = area(cells)
    for c in cells:
        dA = c.area()
        x += c.get(xlabel)*dA
        y += c.get(ylabel)*dA
        z += c.get(zlabel)*dA
    x /= A
    y /= A
    z /= A
    return Vector3(x, y, z)

def oriented_normal(n):
    """Orients the normal in a consistent direction.

    Tecplot does not list the corners of cells in a consistent
    manner, so the orientation (outward or inward facing)
    of the computed normals can differ between adjacent cells.
    Presently, we'll enforce the assumption
    of a positive-x sense to all normals.
    """
    return Vector3(abs(n.x), n.y, n.z)
    

def compute_fluxes(cells, var_map, species, gmodel):
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
        n = oriented_normal(c.normal())
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
        for isp, sp in enumerate(species):
            # Try to grab mass fraction or set to 1.0 if not found
            massf = c.get(sp, 1.0)
            f_sp[isp] = f_sp[isp] + rho*massf*u_n*dA
    return {'mass': f_mass, 'mom': f_mom, 'energy': f_energy, 'species':f_sp}

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
        n = oriented_normal(c.normal())
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
        n = oriented_normal(c.normal())
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
        massf = [ f_isp/f_mass for f_isp in f_sp ]
        set_massf(Q, gmodel, massf)
    else:
        Q.massf[0] = 1.0
    
    LARGE_PENALTY = 1.0e6

    def f_to_minimize(x):
        rho, T, u = x
        if rho < 0.0 or T < 0.0:
            # Give a big penalty
            return LARGE_PENALTY
        # Use equation of state to compute other thermo quantities
        Q.rho = rho
        Q.T[0] = T
        flag = gmodel.eval_thermo_state_rhoT(Q)
        if flag != 0:
            # If there are problems, then these are NOT good values
            # so return a large error
            return LARGE_PENALTY
        h = gmodel.mixture_enthalpy(Q)
        p = Q.p
        # Compute errors
        fmass_err = abs(f_mass - rho*u*A)/(abs(f_mass) + 1.0)
        fmom_err = abs(f_mom_s - (rho*u*u + p)*A)/(abs(f_mom_s) + 1.0)
        fe_err = abs(f_energy - (rho*u*A*(h + 0.5*u*u)))/(abs(f_energy) + 1.0)
        # Total error is the sum
        return fmass_err + fmom_err + fe_err

    # Compute an initial guess based on mass-flux weighted averages
    mfw_props = mass_flux_weighted_avg(cells, ['rho', 'p', 'T', 'u', 'v', 'w', 'M'], var_map)
    u = sqrt(mfw_props['u']**2 + mfw_props['v']**2 + mfw_props['w']**2)
    guess = [mfw_props['rho'], mfw_props['T'], u]
    result = minimize(f_to_minimize, guess, method='Nelder-Mead')
    rho, T, u = result.x
    # Check the results make sense
    if rho < 0.0 or T < 0.0:
        result.success = False

    if not result.success:
        # Try again but use area-weigthed average as starting point
        aw_props =  area_weighted_avg(cells, ['rho', 'p', 'T', 'u', 'v', 'w', 'M'], var_map)
        u = sqrt(aw_props['u']**2 + aw_props['v']**2 + aw_props['w']**2)
        guess = [aw_props['rho'], aw_props['T'], u]
        result = minimize(f_to_minimize, guess, method='Nelder-Mead')
        rho, T, u = result.x
        if not result.success:
            print "The minimizer had difficulty finding the best set of one-d values: rho, T, u"
            print "result [rho, T, u]= ", result.x

        if rho < 0.0 or T < 0.0:
            print "The minimizer claims a successful finish, but the gas state is non-physical."
            print "result [rho, T, u]= ", result.x

    Q.rho = rho; Q.T[0] = T
    gmodel.eval_thermo_state_rhoT(Q)

    phis = dict.fromkeys(props, 0.0)

    for k in phis:
        if k == 'rho':
            phis[k] = rho
        elif k == 'p':
            phis[k] = Q.p
        elif k == 'T':
            phis[k] = T
        elif k == 'M':
            phis[k] = u/Q.a
        elif k == 'u':
            phis[k] = u
        elif k in species:
            isp = gmodel.get_isp_from_species_name(k)
            phis[k] = massf[isp]
        else:
            print "Do not know what to do for flux-conserved average of:", k
    
    return phis
    
    
    

    
    
        

