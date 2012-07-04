"""
nenzfr_stats.py -- Flow statistics functions needed by the main program.

VERSION: 26-06-2012
    Updated to calculate the nozzle exit statistics using conservation
    of mass, momentum and energy. 

    Area and mass weighted statistics are also calculated and written
    out.
"""

import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

from numpy import pi, sqrt, zeros, array
from libprep3 import Vector, vabs, dot, get_gas_model_ptr
from cfpylib.nm.zero_solvers import secant
from e3prep import select_gas_model
from gaspy import *

#---------------------------------------------------------------
def print_stats(sliceFileName,jobName,coreRfraction):

    # Area weighted statistics...
    # Output filename: jobName-exit.statsArea
    print_stats_MoA(sliceFileName,jobName,coreRfraction,'area')

    # Mass weighted statistics...
    # Output filename: jobName-exit.statsMass
    print_stats_MoA(sliceFileName,jobName,coreRfraction,'mass')

    # Conserve mass, momemtum, energy...
    # Output filename: jobName-exit.stats
    print_stats_CMME(sliceFileName,jobName,coreRfraction)
    return

def print_stats_MoA(sliceFileName,jobName,coreRfraction,weight):
    """
    (MoA = Mass- or Area-weighted)
    Display either the mass-weighted or area-weighted statistics
    of flow properties at the nozzle exit.
    """
    print weight.capitalize()+"-Weighted Nozzle-exit statistics:"
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
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
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
            #dA = y1**2 - y0**2
            dA = pi*(y0+y1)*sqrt((y1-y0)**2+(x0-x1)**2)
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
            F += data[var][j] * weighting * dA
            MassFlux += weighting * dA
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
        print "%10s  %12.5g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.5g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return

#-------------------------------------------------------------------
def print_stats_CMME(sliceFileName,jobName,coreRfraction):
    """
    Display statistics of flow properties at the nozzle exit
    using conservation of mass, momentum and energy method. 
    The implementation is based on the paper:
          Baurle, R.A. and Gaffney, R.L. (2008)
          "Extraction of One-Dimensional Flow Properties
          from Multidimensional Data Sets", Journal of
          Propulsion and Power, vol. 24, no. 4, pg. 704
    Equations numbers given in the comments of the code
    refer to the paper.
    """
    print "Nozzle-exit statistics:"
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
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    xs = data['pos.x']
    y_edge = ys[-1] * coreRfraction
    #
    #
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.stats','w')
    fout.write('CoreRadiusFraction: %10.6f\n' % coreRfraction)
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
    #
    # Over the core region of interest we now want to
    # calculate: (1) total area; (2) the unit normal;
    # (3) total mass flux for each species; (4) total
    # momentum flux; and (5) total energy flux.
    #
    # First initialise required values:
    A = 0.0 #...area
    nhatTotal = Vector(0.0) #...unit normal
    speciesKeys = [k for k in variable_list if k.startswith("mass")]
    massFluxes = zeros([len(speciesKeys)])
    momentumFlux = Vector(0.0)
    energyFlux = 0.0
    # We also need to calculate the mass-weighted Mach
    # number for later use. Mass-weighted turbulent
    # parameters will also be reported in the stats file.
    mach = 0.0; totalMassFlow = 0.0
    tke = 0.0; omega = 0.0; dt_chem = 0.0
    # Define the following for ease of use...
    rho = data['rho']; u = data['vel.x']
    v = data['vel.y']; p = data['p']
    h0 = data['total_h'] # this includes tke
    a = data['a']; M = data['M_local']
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
        dA = pi*(y0+y1)*sqrt((y1-y0)**2+(x0-x1)**2)
        A += dA
        # Calculate the outward pointing unit
        # normal vector to the differential area
        # element...
        edge = Vector(x0,y0,0)-Vector(x1,y1,0)
        nhat = Vector(-edge.y, edge.x, 0)/vabs(edge)
        nhatTotal += nhat*dA
        # Velocity vector and weighting function...
        vel = Vector(u[j],v[j],0.0)
        weighting = rho[j]*dot(vel,nhat)
        # Total mass, momentum and energy fluxes...
        # Eq'ns (2a), (2b), and (2c) in paper.
        momentumFlux += (weighting*vel + p[j]*nhat)*dA
        #h0 = e0[j] + p[j]/rho[j] + 0.5*dot(vel,vel)
        #energyFlux += weighting*h0*dA
        energyFlux += weighting*h0[j]*dA
        speciesFractions = [data[species][j] for species in speciesKeys]
        massFluxes += weighting*array(speciesFractions)*dA

        # Mass-weighted parameters...
        mach += M[j]*weighting*dA
        totalMassFlow += weighting*dA
        tke += data['tke'][j]*weighting*dA
        omega += data['omega'][j]*weighting*dA
        dt_chem += data['dt_chem'][j]*weighting*dA

    # Calculate the overall unit normal vector using
    # Eq'n (10) from the paper...
    g_nhat = nhatTotal/A
    # Total Mass flux...Eq'n (A2)
    totalMassFlux = sum(massFluxes)
    # 1D species mass fractions...Eq'n (A4)
    massFrac = massFluxes/totalMassFlux
    # Mass-weighted Mach number...
    aveMach = mach/totalMassFlow
    # Mass-weighted turbulent parameters...
    ave_tke = tke/totalMassFlow #...Eq'n (20) and (A2)
    ave_omega = omega/totalMassFlow
    ave_dt_chem = dt_chem/totalMassFlow

    # Get a list of just the species names (without the "massf[i]-"
    # prefix)...
    speciesList = [name.split('-')[1] for name in speciesKeys]
    speciesDict = dict([(k,v) for k,v in zip(speciesList,massFrac)])
    # Now determine what type of gas model has been used...
    if len(speciesList)==1 and speciesList[0] in ['LUT']:
        print "I don't know how to deal with an equlibrium"+\
              " LUT yet."
        return
    # Create a gas model to use...this won't work to equilibrium LUT(?)
    select_gas_model(fname='gas-model.lua')
    gmodel = get_gas_model_ptr()
    gdata = Gas_data(gmodel)
    # Set mass-fractions...
    set_massf(gdata,gmodel,speciesDict)
    # Calculate gas constant for the mixture...
    R = gmodel.R(gdata)
    # Calculate 'maximum' temperature...
    momentumFluxScalar = dot(momentumFlux,g_nhat) #...Eq'n (A8)
    Tmax = (momentumFluxScalar/totalMassFlux)**2/(4*R) #...Eq'n (A13)

    # We now need to define functions based on Eq'n (A12)
    def energy_error_pos(x,gasData=gdata,gasModel=gmodel,\
                     fpvec=momentumFlux,fp=momentumFluxScalar,\
                     fm=totalMassFlux,fe=energyFlux,tke=ave_tke):
        """
        x    : temperature
        gasData  : gas model data object
        gasModel : gas model object
        fpvec: momentum Flux vector
        fp   : momentum Flux scalar (Eq'n A8)
        fm   : mass flux
        fe   : energy flux
        tke  : mass-averaged tke
        """
        T = x #...1D temperature
        R = gasModel.R(gasData) #...gas constant
        gasData.T[0] = T #...update gas data temperature
        h = gasModel.total_enthalpy(gasData) #...1D static enthalpy
        # NB. 'total_enthalpy' command above is to get sum the static
        # enthalpy of ALL species.
        #
        pFluxPerMass = fp/fm
        #...Eq'n (A11), positive
        vel_dot_nhat = 0.5*( pFluxPerMass + sqrt(pFluxPerMass**2 - 4*R*T) )
        # dot product of velocity with itself...
        vel_dot_vel = 1/(fm*fm)*dot(fpvec,fpvec) - 2*R*T - (R*T/vel_dot_nhat)**2
        #...Eq'n (A12). We explicitly account for tke as per Eq'n (19)...
        error = h + 0.5*vel_dot_vel - fe/fm + tke
        return error
    def energy_error_neg(x,gasData=gdata,gasModel=gmodel,\
                     fpvec=momentumFlux,fp=momentumFluxScalar,\
                     fm=totalMassFlux,fe=energyFlux,tke=ave_tke):
        """
        x    : temperature
        gasData  : gas model data object
        gasModel : gas model object
        fpvec: momentum Flux vector
        fp   : momentum Flux scalar (Eq'n A8)
        fm   : mass flux
        fe   : energy flux
        tke  : mass-averaged tke
        """
        T = x #...1D temperature
        R = gasModel.R(gasData) #...gas constant
        gasData.T[0] = T #...update gas temperature
        h = gasModel.total_enthalpy(gasData) #...1D static enthalpy
        #
        pFluxPerMass = fp/fm
        #...Eq'n (A11), negative
        vel_dot_nhat = 0.5*( pFluxPerMass - sqrt(pFluxPerMass**2 - 4*R*T) )
        # dot product of velocity with itself...
        vel_dot_vel = (1/fm)**2*dot(fpvec,fpvec) - 2*R*T - (R*T/vel_dot_nhat)**2
        #...Eq'n (A12). We explicitly account for tke as per Eq'n (19)...
        error = h + 0.5*vel_dot_vel - fe/fm + tke
        return error
    #
    # Now calculate the two possible temperatures...
    #
    # Positive root...
    Tpos = secant(energy_error_pos,\
                  300.0,1000.0,limits=[0.0,Tmax],tol=1e-6)
    # Negative root...
    Tneg = secant(energy_error_neg,\
                  300.0,1000.0,limits=[0.0,Tmax],tol=1e-6)
    # We may not always find a solution so only construct a dictionary
    # of properties if we have...
    if Tpos not in ['FAIL','Fail']:
        pos, gdataPos = \
           assemble_data(gmodel,gdata,Tpos,'pos',momentumFluxScalar,\
                  momentumFlux,totalMassFlux,energyFlux,A,g_nhat)
    else: # Set a non-sensical value for M_local for later use
        pos = {}; pos['M_local'] = 1000
    #
    if Tneg not in ['FAIL','Fail']:
        neg, gdataNeg = \
           assemble_data(gmodel,gdata,Tneg,'neg',momentumFluxScalar,\
                  momentumFlux,totalMassFlux,energyFlux,A,g_nhat)
    else:
        neg = {}; neg['M_local'] = 1000
    # In case something really goes wrong...
    if Tpos in ['FAIL','Fail'] and Tneg in ['FAIL','Fail']:
        print "Failed to find 1D Temperature that satisfies CMME.\n"+\
              "Something went wrong."
        return

    # Following the recommendation of Baurle and Gaffney (2008) we
    # select the set of properties for which the Mach number is
    # closest to the mass-weighted average Mach number previously
    # calculated.
    if abs(pos['M_local']-aveMach) < abs(neg['M_local']-aveMach):
        properties = pos
        Gas_data.copy_values_from(gdata,gdataPos)
    else:
        properties = neg
        Gas_data.copy_values_from(gdata,gdataNeg)

    # Add the species mass fraction information to the
    # "properties" dictionary...
    for k in range(len(speciesList)):
        properties[speciesKeys[k]] = massFrac[k]

    # Calculate Pitot pressure using Rayleigh formula...
    g = properties['gamma']; M = properties['M_local']
    p = properties['p']
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

    #----------------------------------------------------
    # We now have all the one-dimensionalised flow
    # properties and calculate the statistics and
    # write a summary.
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
        # Caculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%10s  %12.5g    %10.3g %10.3g %10.3g" % \
              (var, properties[var], diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.5g    %10.3g %10.3g %10.3g\n' % \
              (var, properties[var], diff_minus, diff_plus, stddev))
    #
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return

def assemble_data(gasModel,gasData,T,flag,fp,fpvec,fm,fe,A,g_nhat):
    """
    Small function to calculate as many properties as possible
    for a given gas-model based on an input temperature and
    fluxes. Calculation method follows that outlined in Appendix
    A of the paper:
        Baurle, R.A. and Gaffney, R.L. (2008)
        "Extraction of One-Dimensional Flow Properties from
        Mulitdimensional Data Sets", Journal of Propulsion
        and Power, vol. 24, no. 4, pg. 704

    gasModel : libgas object of gas
    gasData  : ligas object of gas properties
    T        : temperature
    flag     : sets whether to use the postive or negative branch
             : of Eq'n (A11)
    fp       : momentum flux scalar value (Eq'n A8)
    fpvec    : momentum flux vector
    fm       : mass flux
    fe       : energy flux
    A        : total area
    g_nhat   : unit normal (Eq'n 10)
    """
    R = gasModel.R(gasData) #...gas constant
    gasData.T[0] = T
    h = gasModel.enthalpy(gasData,0) #...enthalpy

    pFluxPerMass = fp/fm
    if flag in ['neg']:
        vel_dot_nhat = 0.5*( pFluxPerMass - sqrt(pFluxPerMass**2 - 4*R*T) )
    elif flag in ['pos']:
        vel_dot_nhat = 0.5*( pFluxPerMass + sqrt(pFluxPerMass**2 - 4*R*T) )

    rho = fm/vel_dot_nhat/A #...density Eq'n (A9)
    gasData.rho = rho
    p = rho*R*T             #...pressure
    # now check the pressure...
    gasModel.eval_thermo_state_rhoT(gasData)
    pcheck = gasData.p
    if abs(p - pcheck) > 1.0e-4:
        print "something went wrong. pressure from gasModel doesn't"+\
              "match that from gas law"
        return

    vel = (fpvec - p*A*g_nhat)/fm #...velocity vector, Eq'n (A6)
    M = sqrt(dot(vel,vel))/gasData.a    #...Mach number

    gasModel.Cp(gasData)
    gasModel.Cv(gasData)
    gasModel.eval_transport_coefficients(gasData)

    # Create a dictionary of properties to return...
    out = {}
    out['rho'] = rho; out['vel.x'] = vel.x; out['vel.y'] = vel.y
    out['p'] = p; out['a'] = gasData.a;
    out['mu'] = gasData.mu; out['k[0]'] = gasData.k[0]
    out['e[0]'] = gasData.e[0]; out['T[0]'] = T; out['M_local'] = M
    out['total_h'] = fe/fm #...Eq'n (A5)
    out['gamma'] = gasModel.gamma(gasData)
    out['R'] = R
    # Output a copy of gasData for later use...
    gasDataOut = Gas_data(gasModel)
    Gas_data.copy_values_from(gasDataOut,gasData)
    return out, gasDataOut



#------------------------------------------------------------------
def print_stats_orignal(sliceFileName,jobName,coreRfraction):
    """
    (Original implementation by PJ.)
    Display the area-weighted statistics of flow properties 
    at the nozzle exit.
    """
    print "Nozzle-exit statistics:"
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
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    y_edge = ys[-1] * coreRfraction
    #
    # Compute and print area-weighted-average core flow values.
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.statsArea','w')
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
    for var in variable_list:
        if var in exclude_list: continue
        A = 0.0; F = 0.0;
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            if j == 0:
                y0 = 0.0
            else:
                y0 = 0.5*(ys[j-1]+ys[j])
            y1 = 0.5*(ys[j]+ys[j+1])
            dA = y1**2 - y0**2
            F += data[var][j] * dA
            A += dA
        mean = F/A
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
        print "%10s  %12.4g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.4g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return

