#!/usr/bin/env python
# nenzfr_RSA.py
#
# This script either:
#    1) Builds an approximation to the Response Surface
#    for the freestream properties (w.r.t Vs and pe) using 
#    a number of nenzfr cases generated by 
#    "nenzfr_perturbed.py", 
# or
#    2) Using a nominated RSA file, calculates the 
#    freestream properties for given values of (Vs,pe). 
# 
# Luke Doherty
# School of Mechancial and Mining Engineering
# The University of Queensland

VERSION_STRING = "04-Feb-2013"

import string
import sys, os, gzip
import optparse
from numpy import *
from math import copysign
from nenzfr_utils import run_command, quote, read_case_summary, \
     read_nenzfr_outfile, read_estcj_outfile
from nenzfr_sensitivity import add_extra_variables
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------
def calculate_RS_coefficients(exitVar, DictOfCases, nozzleData, RSAtype):
    """
    Function to calculate the coefficients for the chosen response
    surface by solving a set of linear equations using the results of
    a set of perturbation cases.
    """
    beta = {}
    # Loop through each freestream property
    for var in exitVar:
        Y = zeros((len(DictOfCases),1))
        if RSAtype in ['radial']:
            X = zeros((len(DictOfCases),len(DictOfCases)))
        elif RSAtype in ['2nd-order']:
            X = zeros((len(DictOfCases),6))
        
        # Loop through each perturbation case
        k = 0
        X0 = array(DictOfCases['case00'])
        for case in DictOfCases.keys():
            # Store the value of the exit flow property
            # for the current case
            Y[k] = nozzleData[case][var]
            if RSAtype in ['radial']:
                # Using normalised coordinates for (Vs,pe), 
                # calculate the set of Euclidean norms for 
                # the current case with respect to all other
                # cases.
                X_temp = array(DictOfCases[case])/X0 - \
                         array(DictOfCases.values())/X0
                X[k,:] = [sqrt(dot(X_temp[i],X_temp[i])) for i in range(len(X_temp))]
            elif RSAtype in ['2nd-order']:
                X_temp = array(DictOfCases[case])/X0
                X[k,:] = [1.0, X_temp[0], X_temp[1], X_temp[0]**2, X_temp[1]**2,\
                          X_temp[0]*X_temp[1]]
            k += 1
        
        # Now solve the set of linear equations:
        #     Y = X*B
        B = linalg.lstsq(X,Y)
        beta[var] = B[0].transpose() # Make it a row
        #if RSAtype in ['radial']:
        #    B = linalg.solve(X,Y)
        #    beta[var] = B.transpose() # Make it a row
        #elif RSAtype in ['2nd-order']:
        #    #B = dot(dot(linalg.inv(dot(X.transpose(),X)),X.transpose()),Y)
        #    B = linalg.lstsq(X,Y)
        #    beta[var] = B[0].transpose() # Make it a row
    return beta

def write_RSA_file(beta, exitVar, DictOfCases, RSAtype, FileToWrite):
    """
    Write out a file summarising the response surface coefficients
    for each exit flow property.
    """
    fout = open(FileToWrite,'w')
    # Write out the type of RS as the first line. 
    fout.write('{0:>9}'.format(RSAtype))
    fout.write('\n')
    # Write out title line
    fout.write('{0:>12}\t'.format('variable'))
    for case in DictOfCases.keys():
        fout.write('{0:>15}\t'.format(case))
    fout.write('\n')
    # Write out the Vs values for each case
    fout.write('{0:>12}\t'.format('Vs'))
    for values in DictOfCases.values():
        fout.write('{0:>15.5g}\t'.format(values[0]))
    fout.write('\n')
    # Write out the pe values for each case
    fout.write('{0:>12}\t'.format('pe'))
    for values in DictOfCases.values():
        fout.write('{0:>15.6g}\t'.format(values[1]))
    fout.write('\n')
    # Write out a horizontal line
    for k in range(len(DictOfCases)):
        fout.write('{0:->15}'.format('-'))
    fout.write('{0:->24}'.format('-'))
    fout.write('\n')
    # For 2nd-order RS we write out what each beta is for
    if RSAtype in ['2nd-order']:
        fout.write('{0:>12}\t{1:>15}\t{2:>15}\t{3:>15}\t{4:>15}\t{5:>15}\t{6:>15}'.\
                   format('','1','Vs','pe','Vs**2','pe**2','Vs*pe'))
        fout.write('\n')
        for k in range(6):
            fout.write('{0:->15}'.format('-'))
        fout.write('{0:->21}'.format('-'))
        fout.write('\n')    
    # Now write out the radial basis function coefficients
    # for each freestream property
    for var in exitVar:
        fout.write('{0:>12}\t'.format(var))
        for b in beta[var][0]:
            fout.write('{0:>15.7g}\t'.format(b))
        fout.write('\n')
    fout.close()
    return 0

def read_RSA_file(FileToRead):
    """
    Read in a Response Surface file and return:
    
    exitVar: a list of the exit flow properties that have
        been fitted with surfaces;
    DictOfCases: a dictionary detailing the perturbation 
        cases that were used to create the response
        surface. The keys are the case names while the 
        values are a list of the form [Vs, pe];
    RSAtype: the type of response surface
    beta: a dictionary of coefficients. The keys are
        the exit flow property names while the values
        are a list of floats of coefficients. The type
        of response surface dictates how the coefficients
        should be used.
    """
    fp = open(FileToRead,'r')
    
    # Get the type of response surface that we are loading
    RSAtype = fp.readline().strip()
    # Get case names
    titles = fp.readline().strip().split("\t")
    titleList = [k.strip() for k in titles if k!="" and k!="variable"]
    # Get  Vs value for each case
    values = fp.readline().strip().split("\t")
    caseVs = [float(k) for k in values if k!="" and k!="Vs"]
    # Get pe value for each case
    values = fp.readline().strip().split("\t")
    casePe = [float(k) for k in values if k!="" and k!="pe"]
    # Assemble a Dictionary of Cases
    DictOfCases = {}
    for j in range(len(titleList)):
        DictOfCases[titleList[j]] = [caseVs[j],casePe[j]]
    fp.readline() # This is just a line of "-"
    if RSAtype in ['2nd-order']:
        fp.readline() # 
        fp.readline() # This is just a line of "-"
    # Now read the rest of the data and assemble a 
    # dictionary of beta values
    beta = {}
    fileLines = fp.readlines()
    exitVar = []
    for line in fileLines:
        data = line.strip().split("\t")
        values = [k for k in data if k!=""]
        var = values[0]
        exitVar.append(var)
        del values[0]
        beta[var] = [float(values[k]) for k in range(len(values))]
    fp.close()
    return exitVar, DictOfCases, RSAtype, beta

def extrapolation_check(Vs,pe,DictOfCases):
    """
    Compare the current Vs, pe values to the size of the domain
    for which the response surface was designed. If we are 
    outside the domain we print some warning messages so that 
    the user is aware that the response surface will be 
    extrapolating.
    """
    minVs = min(DictOfCases.values(), key=lambda x:x[0])[0]
    maxVs = max(DictOfCases.values(), key=lambda x:x[0])[0]
    
    minPe = min(DictOfCases.values(), key=lambda x:x[1])[1]
    maxPe = max(DictOfCases.values(), key=lambda x:x[1])[1]
     
    if Vs < minVs or Vs > maxVs:
        dVs = min([abs(x) for x in [Vs-maxVs, Vs-minVs]])
        dVs = copysign(dVs,Vs-maxVs)
        print 
        print "WARNING: Vs is outside the design domain of the current response surface\n"+\
              "         by {0:5.1f} m/s or {1:5.1%} of the nominal.".format(dVs, dVs/DictOfCases['case00'][0])
    if pe < minPe or pe > maxPe:
        dpe = min([abs(x) for x in [pe-maxPe, pe-minPe]])
        dpe = copysign(dpe,pe-maxPe)
        print
        print "WARNING: pe is outside the design domain of the current response surface\n"+\
              "         by {0:>7.1f} Pa or {1:>5.1%} of the nominal.".format(dpe, dpe/DictOfCases['case00'][1])
    print
    return

def calculate_freestream(Vs,pe,exitVar,DictOfCases,RSAtype,beta):
    """
    Calculate and return a dictionary of the freestream values 
    corresponding to the given (Vs,pe) using the response surface 
    coefficients (beta) provided. The equation used is dictated
    by the type of surface being used (RSAtype). The freestream
    properties returned correspond to those in the exitVar list. 
    """
    freeStreamValues = {}
    # Loop through each nozzle property
    for var in exitVar:
        # Get the nominal values for (Vs,pe)
        X0 = array(DictOfCases['case00'])
        # Calculate normalised coordinates
        X = array([Vs,pe])/X0
        if RSAtype in ['radial']:
            Xi = array(DictOfCases.values())/X0
            # Calculate Euclidean distance from the 
            # current (Vs,pe) values to each of the cases
            Xdiff = X - Xi
            R = [sqrt(dot(Xdiff[i],Xdiff[i])) for i in range(len(Xdiff))]
        elif RSAtype in ['2nd-order']:
            R = [1.0, X[0], X[1], X[0]**2, X[1]**2,X[0]*X[1]]
        # Calculate the nozzle property
        freeStreamValues[var] = sum(R*array(beta[var]))
    return freeStreamValues    

def calculate_residuals(freeStreamValues, nozzleData, exitVar, \
                        RSAtype=None, beta=None, DictOfCases=None):
    """
    Function to calculate the residuals (differences between predicted
    freestream values and those returned by nenzfr.

    This function calculates the residuals for two situations. 
    (1) freeStreamValues are given. In this case we are comparing
        the input freeStreamValues to the input nozzleData for 
        each property in the exitVar list.
    (2) no freeStreamValues are given. In this case we assume that
        we are calculating residuals for each perturbation case that
        was used to create the response surface. We must therefore
        specify the type of response surface (RSAtype), the
        coefficients (beta) and the perturbation cases (DictOfCases).

    The function returns the dictionary 'residuals'. The keys of
    this dictionary correspond to the freestream property names in 
    exitVar. The dictionary values are either a single float or a
    list of floats (one value for every perturbation case).
    """
    residuals = {}
    if freeStreamValues is not None:
        # Loop through each exit flow property
        for var in exitVar:
            residuals[var] = (freeStreamValues[var]-nozzleData[var])\
                                                   /nozzleData[var]*100.0
    else:
        # Calculate the residuals for every case
        if beta is None:
            print "Bad input (beta) for 'calculate_residuals'"
            return -2
        if DictOfCases is None:
            print "Bad input (DictOfCases) for 'calculate_residuals'"
            return -2
        if RSAtype is None:
            print "Bad input (RSAtype) for 'calculate_residuals'"
            return -2
        
        # Calculate freestream properties using the response surface
        # for each case
        rsaFreestream = {}
        for case in DictOfCases.keys():
            Vs = DictOfCases[case][0]
            pe = DictOfCases[case][1]
            
            rsaFreestream[case] = calculate_freestream(Vs,pe,exitVar,\
                                                DictOfCases,RSAtype,beta)
            
                
        # Now calculate the residuals (as percentages) and put them into 
        # a dictionary that has a similar structure as to the "beta" 
        # dictionary.
        for var in exitVar:
            residuals[var] = []
            for case in DictOfCases.keys():
                # Return the residuals as PERCENTAGE values
                res = (rsaFreestream[case][var]-nozzleData[case][var])\
                                               /nozzleData[case][var]*100.0
                residuals[var].append(res)
            # The following line is a work around that allows me to later use 
            # the "write_RSA_file" function to write the residuals to a file
            residuals[var] = [residuals[var]]
    return residuals


def write_flow_summary(Vs,pe,valuesDict,outFileName,exitVar):
    """
    Write out a file summarising the predicted freestream
    values.
    """
    fout = open(outFileName,'w')
    
    # Write title line and values of Vs, pe used. Underline
    # this header.
    fout.write('{0:>12}{1:>12}'.format('variable','value'))
    fout.write('\n')
    fout.write('{0:>12}{1:>12.5g}'.format('Vs',Vs))
    fout.write('\n')
    fout.write('{0:>12}{1:>12.6g}'.format('pe',pe))
    fout.write('\n')
    fout.write('{0:->24}'.format('-'))
    fout.write('\n')
    # Loop through each nozzle property 
    for var in exitVar: 
        fout.write('{0:>12}'.format(var))
        fout.write('{0:>12.5g}'.format(valuesDict[var]))
        fout.write('\n')
    fout.close()

def main():
    """
    Examine the command-line options to decide the what to do
    and then either build a Response Surface Approximation or 
    use a nominated RSA and Vs, and pe values to calculate 
    new freestream property values.
    """
    op = optparse.OptionParser(version=VERSION_STRING)

    op.add_option('--create-RSA', dest='createRSA', action='store_true',
                  default=False, help="create the Response Surface "
                  "approximation using perturbation results generated by "
                  "'nenzfr_perturbed.py'. [default: %default]")

    op.add_option('--RSA-type', dest='RSAtype', choices=['radial','2nd-order'],
                  default='2nd-order',help="specify whether the response "
                  "surface should be formed from a summation of radial basis "
                  "functions or from a 2nd-order equation in Vs and pe, "
                  "choices: radial, 2nd-order. [default: %default]")

    op.add_option('--calculate-residuals', dest='calcResiduals', action='store_true',
                  default=False, help="calculate difference between a RSA freestream "
                  "estimate and actual nenzfr results. [default: %default] ")
    
    op.add_option('--RSA-file', dest='RSAfile', default='response_surface.dat',
                  help="specify the name of the reponse surface file that "
                  "is to be either created or used/read. [default: %default]")

    op.add_option('--exitStatsFile', dest='exitStatsFileName',
                  default='nozzle-exit.stats',
                  help="file that holds the averaged nozzle-exit "
                       "data and is to be read in for each perturbation "
                       "case [default: %default]")

    op.add_option('--estcjFile', dest='estcjFile', default='nozzle-estcj.dat',
                  help="file that holds the estcj result and is to be read in "
                       "for each perturbation case. [default: %default]")
    
    op.add_option('--Vs', dest='Vs', default=None, type='float',
                  help=("incident shock speed, in m/s. [default: %default]"))
    op.add_option('--pe', dest='pe', default=None, type='float',
                  help=("equilibrium pressure (after shock reflection), in Pa. "
                        "[default: %default]"))
    op.add_option('--exitFile', dest='exitFileName', default='nozzle-exit.RSAdat',
                  help="file for holding the RSA calculated nozzle exit data "
                       "[default: %default]")
    op.add_option('--add-extra-variables', dest='addExtraVariables', action='store_true',
                  default=False, help=("specify whether q, rho*u_x, Re_u and p/q should "
                       "also be calculated. Not used when creating an RSA. "
                       "[default: %default]"))

    opt, args = op.parse_args()
        
    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    if not opt.createRSA:
        if opt.Vs is None:
            print "Need to supply a value for Vs."
            bad_input = True    
        if opt.pe is None:
            print "Need to supply a value for pe."
            bad_input = True
        if opt.calcResiduals:
            if not os.path.exists(opt.estcjFile):
                print "'"+opt.estcjFile+"' does not exist in current directory"
                bad_input = True
            if not os.path.exists(opt.exitStatsFileName):
                print "'"+opt.exitStatsFileName+"' does not exist in current directory"
                bad_input = True
    if opt.createRSA:
        if not os.path.exists('perturbation_cases.dat'):
            print "'perturbation_cases.dat' does not exist in current directory"
            bad_input = True    
        if opt.addExtraVariables is True:
            print "Ignoring option --add-extra-variables for calculation of RSA"
            opt.addExtraVariables = False
    if bad_input:
        return -2
    
    if opt.createRSA: 
        """
        Create a Response Surface Approximation
        """
        # Read the "perturbation_cases.dat" file
        perturbedVariables, DictOfCases = read_case_summary()
        
        # Load in all the freestream data and the supply temperature
        # and enthalpy
        nozzleData = {}
        for case in DictOfCases.keys():
            nozzleData[case], exitVar = \
                 read_nenzfr_outfile('./'+case+'/'+opt.exitStatsFileName)
            supply = read_estcj_outfile('./'+case+'/'+opt.estcjFile)
            nozzleData[case]['supply_T'] = supply['T']
            nozzleData[case]['supply_h'] = supply['h']
        exitVar.insert(0,'supply_T')
        exitVar.insert(1,'supply_h')
        
        # Calculate the basis function coefficients
        beta = calculate_RS_coefficients(exitVar, DictOfCases, \
                                                    nozzleData, opt.RSAtype)
        # Write out a file summarising the coefficients
        write_RSA_file(beta, exitVar, DictOfCases, opt.RSAtype, opt.RSAfile)
        
        if opt.calcResiduals:
            """ 
            Check the quality of the response surface by calculating 
            the residuals for the cases that were used to form the 
            response surface. Unless something went wrong, these 
            should all be very small.
            """
            ld_exitVar, ld_DictOfCases, opt.RSAtype, ld_beta = \
                                                  read_RSA_file(opt.RSAfile)
            
            residuals = calculate_residuals(None, nozzleData, exitVar, \
                                          opt.RSAtype, ld_beta, DictOfCases)
            
            write_RSA_file(residuals, exitVar, DictOfCases, \
                         opt.RSAtype+' RS residuals as PERCENTAGES', opt.RSAfile+'_residuals')
        
    else: 
        """
        Use the RSA to predict new freestream properties
        """
        # Load in the nominated file
        exitVar, DictOfCases, opt.RSAtype, beta = read_RSA_file(opt.RSAfile)
        
        # Check if the given Vs,pe are within the design domain of the
        # current RSA.
        extrapolation_check(opt.Vs,opt.pe,DictOfCases)

        # Calculate nozzle property values
        freeStreamValues = calculate_freestream(opt.Vs,opt.pe,exitVar,\
                                               DictOfCases,opt.RSAtype,beta)
        
        # If desired we add dynamic pressure, mass flux, unit Reynolds number
        # and p/q to the calculated freestream variables
        if opt.addExtraVariables:
            freeStreamValues, exitVar = \
                               add_extra_variables(freeStreamValues, exitVar)
        
        # Write an output file
        write_flow_summary(opt.Vs,opt.pe,freeStreamValues,\
                                                 opt.exitFileName,exitVar)
        
        if opt.calcResiduals:
            """
            Compare the freestream values calculated using the RS to
            actual simulation results.
            """
            nenzfrData, exitVar = \
                        read_nenzfr_outfile(opt.exitStatsFileName)
            nenzfrSupply = read_estcj_outfile(opt.estcjFile)
            nenzfrData['supply_T'] = nenzfrSupply['T']
            nenzfrData['supply_h'] = nenzfrSupply['h']
            exitVar.insert(0,'supply_T')
            exitVar.insert(1,'supply_h')
            
            if opt.addExtraVariables:
                nenzfrData, exitVar = add_extra_variables(nenzfrData, exitVar)
            
            residuals = calculate_residuals(freeStreamValues,nenzfrData,\
                                                                 exitVar)
            write_flow_summary(opt.Vs,opt.pe,residuals,\
                                   opt.exitFileName+'_residuals',exitVar) 


    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Sensitivity:\n Calculate Shock Tunnel Test Flow Conditions\nusing a Response Surface Approximation"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
