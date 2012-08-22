"""
nenzfr_stats.py -- Flow statistics functions needed by the main program.

VERSION: 
26-06-2012
    Updated to calculate the nozzle exit statistics using conservation of
    mass, momentum and energy. 
    Area and mass weighted statistics are also calculated and written out.
22-Aug-2012
    Brute-force approach to estimating the CMME 1D flow conditions.
"""

import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

from numpy import pi, sqrt, zeros, array
from cfpylib.nm.zero_solvers import secant
from cfpylib.nm.nelmin import minimize
from libprep3 import Vector, vabs, dot
from libprep3 import create_gas_model, Gas_data, set_massf

#---------------------------------------------------------------

def print_stats(sliceFileName,jobName,coreRfraction,gmodelFile):

    # Area weighted statistics...
    # Output filename: jobName-exit.statsArea
    print_stats_MoA(sliceFileName,jobName,coreRfraction,'area')

    # Mass weighted statistics...
    # Output filename: jobName-exit.statsMass
    print_stats_MoA(sliceFileName,jobName,coreRfraction,'mass')

    # Conserve mass, momemtum, energy...
    # Output filename: jobName-exit.stats
    print_stats_CMME(sliceFileName,jobName,coreRfraction,gmodelFile)
    return


def get_slice_data(sliceFileName):
    """
    Returns the list of variables found in the slice and
    a dictionary containing the flow data across a single slice.
    For each variable, there is a list of data values at points across the slice.
    """
    fp = open(sliceFileName, 'r')
    # Keep a list of variables in order of appearance.
    varLine = fp.readline().strip()
    items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # print "variable_list=", variable_list
    # Store the data in lists against these names.
    data = {}
    for var in variable_list:
        data[var] = []
    for line in fp.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        assert len(items) == len(variable_list)
        for i in range(len(items)):
            data[variable_list[i]].append(float(items[i]))
    fp.close()
    return variable_list, data


def print_stats_MoA(sliceFileName,jobName,coreRfraction,weight):
    """
    (MoA = Mass- or Area-weighted)
    Display either the mass-weighted or area-weighted statistics
    of flow properties at the nozzle exit.
    """
    print weight.capitalize()+"-Weighted Nozzle-exit statistics:"
    variable_list, data = get_slice_data(sliceFileName)
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    xs = data['pos.x']
    y_edge = ys[-1] * coreRfraction
    #
    # Compute and print area-weighted-average core flow values.
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.stats'+weight.capitalize(),'w')
    fout.write('CoreRadiusFraction: %10.6f\n' % coreRfraction)
    fout.write('%15s  %12s   %10s  %10s %10s\n' % 
               ("variable","mean-value","minus","plus","std-dev"))
    fout.write(65*'-'+'\n')
    #
    print "%15s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 65*'-'
    u = data['vel.x']; v = data['vel.y']; rho = data['rho']
    for var in variable_list:
        if var in exclude_list: continue
        MassFlux = 0.0; F = 0.0;
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            y1 = 0.5*(ys[j]+ys[j+1])
            x1 = 0.5*(xs[j]+xs[j+1])
            if j == 0:
                y0 = 0.0
                dx = xs[j]-x1
                x0 = xs[j]+dx
            else:
                y0 = 0.5*(ys[j-1]+ys[j])
                x0 = 0.5*(xs[j-1]+xs[j])
            # Area element...
            #d_Area = y1**2 - y0**2
            d_Area = pi*(y0+y1)*sqrt((y1-y0)**2+(x0-x1)**2)
            # Unit normal vector...
            edge = Vector(x0,y0,0.0)-Vector(x1,y1,0.0)
            nhat = Vector(-edge.y, edge.x, 0.0)/vabs(edge)
            # Velocity vector...
            vel = Vector(u[j],v[j],0.0)
            # Weighting factor...
            if weight in ['area']:
                weighting = 1.0
            elif weight in ['mass']:
                weighting = rho[j]*dot(vel,nhat)
            else:
                print "Unknown weighting option given"
                break
            # Accumulate total...
            F += data[var][j] * weighting * d_Area
            MassFlux += weighting * d_Area
        mean = F/MassFlux
        # Identify low and high values.
        diff_minus = 0.0
        diff_plus = 0.0
        count = 0.0
        stddev = 0.0
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            diff = data[var][j] - mean
            diff_minus = min(diff, diff_minus)
            diff_plus = max(diff, diff_plus)
            count += 1
            stddev += diff**2
        # Calculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%15s  %12.5g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%15s  %12.5g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #
    print 65*'-'
    fout.write(65*'-'+'\n')
    fout.close()
    return

#-------------------------------------------------------------------

def print_stats_CMME(sliceFileName,jobName,coreRfraction,gmodelFile):
    """
    Display statistics of flow properties at the nozzle exit using
    conservation of mass, momentum and energy (CMME) method.
 
    The implementation is loosely based on the paper:
    Baurle, R.A. and Gaffney, R.L. (2008)
    "Extraction of One-Dimensional Flow Properties from Multidimensional Data Sets",
    Journal of Propulsion and Power, vol. 24, no. 4, pg. 704

    Equations numbers given in the comments of the code refer to the paper.

    PJ, 22-Aug-2012: Use brute-force to get the effective 1D flow properties.
    """
    print "Nozzle-exit statistics (CMME):"
    variable_list, data = get_slice_data(sliceFileName)
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    xs = data['pos.x']
    y_edge = ys[-1] * coreRfraction
    #
    # Over the core region of interest we now want to calculate: 
    # (1) total area; 
    # (2) the unit normal;
    # (3) total mass flux for each species; 
    # (4) total momentum flux; and 
    # (5) total energy flux.
    #
    # First, we need to do a bit of fiddling with variable names in order
    # to get the species names.
    speciesKeys = [k for k in variable_list if k.startswith("mass")]
    # Get a list of just the species names (without the "massf[i]-" prefix)...
    speciesNames = [name.split('-')[1] for name in speciesKeys]
    nsp = len(speciesKeys)
    #
    # Initialise the totals across the slice:
    Area = 0.0 #...area
    nhatTotal = Vector(0.0) #...unit normal
    speciesMassFluxes = zeros((nsp,))
    massFlux = 0.0
    momentumFlux = Vector(0.0)
    energyFlux = 0.0
    # We also need to calculate the mass-weighted Mach number for later use. 
    mach = 0.0
    # Mass-weighted turbulent parameters will also be reported in the stats file.
    tke = 0.0
    omega = 0.0
    dt_chem = 0.0
    #
    # Define the following for ease of use in the integration process below.
    rho = data['rho']
    u = data['vel.x']
    v = data['vel.y']
    p = data['p']
    h0 = data['total_h'] # this includes tke
    a = data['a']
    M = data['M_local']
    #
    # Integrate the conserved properties across the slice.
    for j in range(len(ys)):
        if ys[j] > y_edge: break
        y1 = 0.5*(ys[j]+ys[j+1])
        x1 = 0.5*(xs[j]+xs[j+1])
        if j == 0:
            y0 = 0.0
            dx = xs[j]-x1
            x0 = xs[j]+dx
        else:
            y0 = 0.5*(ys[j-1]+ys[j])
            x0 = 0.5*(xs[j-1]+xs[j])
        # We assume the data is axi-symmetric but no
        # assumption is made about the straightness
        # of the plane over which we are integrating.
        # Hence we use a truncated cone for the area.
        d_Area = pi*(y0+y1)*sqrt((y1-y0)**2+(x0-x1)**2)
        Area += d_Area
        # Calculate the outward pointing unit
        # normal vector to the differential area
        # element...
        edge = Vector(x0,y0,0)-Vector(x1,y1,0)
        nhat = Vector(-edge.y, edge.x, 0)/vabs(edge)
        nhatTotal += nhat*d_Area
        # Velocity vector and weighting function...
        vel = Vector(u[j],v[j],0.0)
        weighting = rho[j]*dot(vel,nhat)
        # Total mass, momentum and energy fluxes...
        # Eq'ns (2a), (2b), and (2c) in paper.
        momentumFlux += (weighting*vel + p[j]*nhat)*d_Area
        #h0 = e0[j] + p[j]/rho[j] + 0.5*dot(vel,vel)
        #energyFlux += weighting*h0*d_Area
        energyFlux += weighting*h0[j]*d_Area
        speciesFractions = array([data[species][j] for species in speciesKeys])
        speciesMassFluxes += weighting*speciesFractions*d_Area
        massFlux += weighting*d_Area
        #
        # Mass-weighted parameters...
        mach += M[j]*weighting*d_Area
        if 'tke' in data.keys():
            tke += data['tke'][j]*weighting*d_Area
        if 'omega' in data.keys():
            omega += data['omega'][j]*weighting*d_Area
        if 'dt_chem' in data.keys():
            dt_chem += data['dt_chem'][j]*weighting*d_Area
    # Check the total mass flux...Eq'n (A2)
    assert abs(massFlux - sum(speciesMassFluxes)) <= 1e-7
    #
    # Now, compute the effective 1D flow properties
    # to have the same integral properties.
    #
    # Calculate the overall unit normal vector using Eq'n (10)
    g_nhat = nhatTotal/Area
    # and a scalar momentum flux of the effective 1D flow, Eq'n (A8)
    momentumFluxScalar = dot(momentumFlux,g_nhat)
    # 1D species mass fractions...Eq'n (A4)
    massFrac = speciesMassFluxes/massFlux
    # Mass-weighted Mach number
    aveMach = mach/massFlux
    # and mass-weighted turbulent parameters.
    ave_tke = tke/massFlux #...Eq'n (20) and (A2)
    ave_omega = omega/massFlux
    ave_dt_chem = dt_chem/massFlux
    # and reconstruct the gas state
    gmodel = create_gas_model(gmodelFile)
    gdata = Gas_data(gmodel)
    massfDict = dict([(k,v) for k,v in zip(speciesNames,massFrac)])
    set_massf(gdata,gmodel,massfDict)
    gdata.rho = data['rho'][0]
    gdata.T[0] = data['T[0]'][0]
    gmodel.eval_thermo_state_rhoT(gdata)
    vx = massFlux/(gdata.rho*Area)
    # print "Before optimizer, vx=", vx, "p=", gdata.p, "T=", gdata.T[0], "rho=", gdata.rho
    #
    def error_estimate(params, 
                       gasData=gdata, gasModel=gmodel,
                       fm=massFlux, fp=momentumFluxScalar, fe=energyFlux,
                       tke=ave_tke, Area=Area, mach=aveMach):
        """
        Estimate the badness the current guess for 1D flow properties.
        params : [rho, T, v] 1D flow properties to be evaluated
        """
        rho, T, vx = params
        gasData.T[0] = T
        gasData.rho = rho
        gmodel.eval_thermo_state_rhoT(gasData)
        p = gasData.p
        h = gasModel.total_enthalpy(gasData) #...1D static enthalpy
        M = vx/gasData.a
        # Realtive errors in each conserved quantity.
        # The "+1.0" items are to avoid (unlikely) problems with zero values
        # for the given quantities.
        fm_err = abs(fm - rho*vx*Area)/(abs(fm)+1.0)
        fp_err = abs(fp - (rho*vx*vx*Area + p*Area))/(abs(fp)+1.0)
        fe_err = abs(fe - rho*vx*Area*(h+0.5*vx*vx))/(abs(fe)+1.0)
        mach_err = abs(M - mach)/(abs(mach)+1.0)
        # The overall error estimate is a weighted sum.
        return fm_err + fp_err + fe_err + 0.1*mach_err
    #
    print "Optimize estimate of 1D flow properties"
    flow_params, fx, conv_flag, nfe, nres = minimize(error_estimate, 
                                                     [gdata.rho, gdata.T[0], vx],
                                                     [0.01, 10.0, 10.0])
    rho, T, vx = flow_params
    # print "fx=", fx
    # print "convergence-flag=", conv_flag
    # print "number-of-fn-evaluations=", nfe
    # print "number-of-restarts=", nres
    gdata.T[0] = T
    gdata.rho = rho
    gmodel.eval_thermo_state_rhoT(gdata)
    p = gdata.p
    # print "After optimizer, vx=", vx, "p=", p, "T=", T, "rho=", rho
    g = gmodel.gamma(gdata)
    R = gmodel.R(gdata)
    M = vx/gdata.a
    total_h = energyFlux/massFlux
    #
    gmodel.eval_transport_coefficients(gdata)
    properties = {}
    properties['rho'] = rho
    properties['vel.x'] = vx
    properties['vel.y'] = 0.0
    properties['p'] = p
    properties['a'] = gdata.a
    properties['mu'] = gdata.mu
    properties['k[0]'] = gdata.k[0]
    properties['e[0]'] = gdata.e[0]
    properties['T[0]'] = T
    properties['M_local'] = M
    properties['total_h'] = total_h
    properties['gamma'] = g
    properties['R'] = R
    for k in range(nsp):
        properties[speciesKeys[k]] = massFrac[k]
    # Calculate Pitot pressure using Rayleigh formula...
    if M > 1.0:
        t1 = (g+1)*M*M/2
        t2 = (g+1)/(2*g*M*M - (g-1));
        pitot_p = p * pow(t1,(g/(g-1))) * pow(t2,(1/(g-1)));
    else:
        t1 = 1 + 0.5*(g-1)*M*M
        pitot_p = p * pow(t1,(g/(g-1)))
    properties['pitot_p'] = pitot_p
    # Calculate total pressure (as isentropic process)...
    t1 = 1 + 0.5*(g-1)*M*M
    total_p = p * pow(t1,(g/(g-1)))
    properties['total_p'] = total_p
    #
    # Calculate the turbulent viscosity and thermal
    # conductivity using definitions given in:
    #    Chan, W.Y.K., Jacobs, P.A., Nap, J.P., Mee, D.J.,
    #    Kirchhartz, R.M. (2011)
    #    "The k-w turbulence model in Eilmer3: User guide
    #    test cases", Research Report Number 2010/11
    properties['mu_t'] = properties['rho']*ave_tke/ave_omega
    Pr_t = 8.0/9.0 #...turbulent Prandtl number.
    Cp = gmodel.Cp(gdata)
    properties['k_t'] =  Cp*properties['mu_t']/Pr_t
    # Finally, add the remaining mass-weighted turbulent
    # and chemistry parameters...
    properties['tke'] =  ave_tke
    properties['omega'] = ave_omega
    properties['dt_chem'] = ave_dt_chem
    #
    # Now we have all the one-dimensionalised flow properties,
    # calculate the statistics and write a summary.
    #
    fout = open(jobName+'-exit.stats','w')
    fout.write('CoreRadiusFraction: %10.6f\n' % coreRfraction)
    fout.write('%15s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(65*'-')
    fout.write('\n')
    print "%15s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 65*'-'
    #
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    for var in variable_list:
        if var in exclude_list: continue
        # Identify the low and high values
        diff_minus = 0.0
        diff_plus = 0.0
        count = 0.0
        stddev = 0.0
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            diff = data[var][j] - properties[var]
            diff_minus = min(diff, diff_minus)
            diff_plus = max(diff, diff_plus)
            count += 1
            stddev += diff**2
        # Calculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%15s  %12.5g    %10.3g %10.3g %10.3g" % \
              (var, properties[var], diff_minus, diff_plus, stddev)
        fout.write('%15s  %12.5g    %10.3g %10.3g %10.3g\n' % \
              (var, properties[var], diff_minus, diff_plus, stddev))
    #
    print 65*'-'
    fout.write(65*'-'+'\n')
    fout.close()
    return
