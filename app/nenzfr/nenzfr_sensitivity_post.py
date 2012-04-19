#!/usr/bin/env python
# nenzfr_sensitivity_post.py
#
# This script...
# 
# Luke Doherty
# School of Mechancial and Mining Engineering
# The University of Queensland

VERSION_STRING = "17-April-2012"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
#import numpy
from numpy import array, mean, logical_and, zeros
import copy
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------
    
def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    print "About to run cmd:", cmdText
    args = shlex.split(cmdText)
    p = subprocess.Popen(args)
    # wait until the subprocess is finished
    stdoutData, stderrData = p.communicate() 
    return

def quote(str):
    """
    Put quotes around a string.
    """
    return '"' + str + '"'

def read_case_summary():
    """
    Reads the file "sensitivity_cases.dat" to determine which variables
    have been perturbed and values for the perturbed variables for each
    case consituting the sensitivity calculation.
    
    :returns: perturbedVariables - a list of variable names
            : DictOfCases - a dictionary with the case names as keys 
                            and the values of each perturbed variable.
    """
    fp = open('sensitivity_cases.dat','r')
    
    varList = fp.readline().strip().split(" ")
    perturbedVariables = [k for k in varList if k!="#" and k!=""]
    fp.readline()
    DictOfCases = {}
    for line in fp.readlines():
        caseData = line.strip().split(" ")
        caseName = caseData[0]
        DictOfCases[caseName] = [float(k) for k in caseData if k!=caseName and k!=""]
    fp.close()
    return perturbedVariables, DictOfCases

def read_nenzfr_outfile(FileToRead):
    """
    Reads the nenzfr out-file containing the statistics of the 
    exit flow properties.
    
    :FileToRead: the name of the file to be read. Default should be 
                 something like "nozzle-exit.stats"
    :returns: a list of all the exit flow properties and a dictionary
              of the mean-values.
    """
    fp = open(FileToRead,'r')
    # Collumn titles
    titles = fp.readline().strip().split(" ")
    titleList = [k for k in titles if k!="" and k!="variable"]
    
    fp.readline() # This is a row of "-"
    exitDataDict = {}
    fileLines = fp.readlines()
    del fileLines[-1] # Get rid of the last line which is just a row of "-"
    exitProperty = []
    # Now read in the rest of the data
    for line in fileLines:
        data = line.strip().split(" ")
        values = [k for k in data if k!=""]
        variable = values[0]
        exitProperty.append(variable)
        exitDataDict[variable] = float(values[1])
        #variable = data[0]
        #exitProperty.append(variable)
        #exitDataDict[variable] = [float(k) for k in data if k!=variable and k!=""]
    fp.close()
    return exitDataDict, exitProperty

def read_estcj_outfile(FileToRead):
    """
    Read just the line in the estcj output file that contains the nozzle
    supply temperature and enthalpy.
    """
    supplyDict = {}
    lineToRead = 0
    
    fp = open(FileToRead,'r')
    for line in fp.readlines():
        # We only want to read one line in the file
        if lineToRead == 1: 
            for data in line.split(','):
                # Data is a list of strings with one string for each property.
                # The form of the strings is 'property: value unit'. We split
                # list and create "values" which is a list of lists of the form
                #     values = [[property],[value, unit]]
                # We then put this into a dictionary.
                values = [x.split() for x in data.split(':')]
                supplyDict[values[0][0]] = float(values[1][0])
            lineToRead = 0
            
        if line.strip() in ['State 5s: equilibrium condition (relaxation to pe)',]:
            # We want to read the next line
            lineToRead = 1
    fp.close()
    return supplyDict

def get_values(dict, propertyList):
    """
    Adapted from Peter Jacob's function "get_fractions" which may be
    found within "cea2_gas.py".
    """
    valueList = []
    #print propertyList
    for s in propertyList:
        #print s
        if s in dict.keys():
            #print s
            #print dict[s]
            valueList.append(dict[s])
        else:
            valueList.append(0.0)
            print "WARNING: "+s+"was not found in the current case dictionary."
    return valueList

def write_sensitivity_summary(sensitivity, perturbedVariables, exitVar, nominalData, abs=0):
    """
    Write out a file summarising the sensitivity of each exit flow parameter to each of
    the input parameters. Depending on the value of the "abs" flag, the values are written
    as either absolute or relative values (to the nominal exit flow property value).
    """
    titleFormatDict = {'p1':'{0:{fill}>13}', 'T1':'{0:{fill}>10}', 'Vs':'{0:{fill}>11}',
                       'pe':'{0:{fill}>15}', 'Tw':'{0:{fill}>10}', 'BLTrans':'{0:{fill}>10}',
                       'TurbVisRatio':'{0:{fill}>14}', 'TurbInten':'{0:{fill}>11}',
                       'CoreRadiusFraction':'{0:{fill}>20}'}
    
    if abs == 1:
        fout = open('sensitivities_abs.dat','w')
        
        formatDict = {'p1':'{0:13.2f}', 'T1':'{0:>10.2f}', 'Vs':'{0:>11.2f}',
                      'pe':'{0:>15.2f}', 'Tw':'{0:>10.2f}', 'BLTrans':'{0:>10.2f}',
                      'TurbVisRatio':'{0:>14.2f}', 'TurbInten':'{0:>11.2f}',
                      'CoreRadiusFraction':'{0:>20.2f}'}
    else:
        fout = open('sensitivities_rel.dat','w')
        
        formatDict = {'p1':'{0:13.2%}', 'T1':'{0:>10.2%}', 'Vs':'{0:>11.2%}',
                      'pe':'{0:>15.2%}', 'Tw':'{0:>10.2%}', 'BLTrans':'{0:>10.2%}',
                      'TurbVisRatio':'{0:>14.2%}', 'TurbInten':'{0:>11.2%}',
                      'CoreRadiusFraction':'{0:>20.2%}'}
    
    # Write header information
    fout.write('{0:>13}'.format('variable'))
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format(k,fill=''))
    fout.write('{0:>10}'.format('total'))
    fout.write('\n')
    for k in perturbedVariables:
        fout.write(titleFormatDict[k].format('-',fill='-'))
    fout.write('{0:->23}'.format('-'))
    fout.write('\n')
    
    for k in exitVar:
        fout.write('{0:>13}'.format(k))
        X_total = 0.0
        for kk in perturbedVariables:
            X_kk = sensitivity[kk][exitVar.index(k)]
            if abs != 1:
                X_kk = X_kk/nominalData[k]
            
            #X_kk = sensitivity[kk][exitVar.index(k)]/nominalData[k] *100.0
            fout.write(formatDict[kk].format(X_kk))
            X_total += X_kk**2
        if abs == 1:
            fout.write('{0:10.2f}'.format(X_total**0.5))
        else:
            fout.write('{0:10.2%}'.format(X_total**0.5))
        fout.write('\n')
    
    fout.close()


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by Nenzfr.
    """
    op = optparse.OptionParser(version=VERSION_STRING)

    op.add_option('--runCMD', dest='runCMD', default='./',
                  help=("command used to execute the run-script file "
                        "[default: %default]"))
    op.add_option('--Cluster', dest='Cluster', default='Mango',
                  choices =['Mango', 'Barrine'],
                  help=("specify on which cluster the computations are to be ran. "
                        "This is used to define which run template script will "
                        "be used. Options: "
                        "Mango; Barrine [default: %default]"))
    
    op.add_option('--gradient', dest='gradient', default='linear',
                  choices=['linear','quadratic'],
                  help=("specify whether the gradient is to be calculated "
                        "using a linear or quadratic equation. Quadratic is "
                        "only  available with the --create-LUT option. "
                        "Options: linear, quadratic [default: %default]"))
    
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")

    op.add_option('--exitStatsfile', dest='exitStatsFileName',
                  default='nozzle-exit.stats',
                  help="file that holds the averaged nozzle-exit data and is too be "
                       "read in for each perturbation case [default: %default]")
     
    op.add_option('--create-LUT', dest='createLUT', action='store_true',
                  default=False, 
                  help="create a LUT by perturbing only p1, T1, Vs and pe")

    opt, args = op.parse_args()

    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    if opt.createLUT is False:
        if opt.gradient == "quadratic":
            opt.gradient = 'linear'
            print "Sensitivity will be calculated using a linear slope."
    if bad_input:
        return -2
    
    # Read the sensitivity_case_summary file to get the perturbed variables and 
    # their various values
    perturbedVariables, DictOfCases = read_case_summary()
    
    # Define the name of the nominal case and load the exit plane data
    nominal = 'case000'
    nominalData, exitVar = read_nenzfr_outfile('./'+nominal+'/'+opt.exitStatsFileName)
    # Load the nozzle supply data
    nominalSupply = read_estcj_outfile('./'+nominal+'/'+opt.jobName+'-estcj.dat')
    # Now add the relevant supply data (T, h) to the nominalData dictionary
    nominalData['supply_T'] = nominalSupply['T']
    nominalData['supply_h'] = nominalSupply['h']
    # Add the supply variables to the exitVar list
    exitVar.insert(0,'supply_T')
    exitVar.insert(1,'supply_h')
    
    nominalValues = get_values(nominalData, exitVar)
    
    # Loop through each of the perturbed variables
    sensitivity = {}
    for k in range(len(perturbedVariables)):
        var = perturbedVariables[k]
        
        # Define the name of the relevant perturbed cases and load the 
        # associated data
        high = 'case'+"{0:02}".format(k)+'1'
        highData, dontNeed  = read_nenzfr_outfile('./'+high+'/'+opt.exitStatsFileName)
        highSupply = read_estcj_outfile('./'+high+'/'+opt.jobName+'-estcj.dat')
        highData['supply_T'] = highSupply['T']
        highData['supply_h'] = highSupply['h']
        
        low = 'case'+"{0:02}".format(k)+'2'
        lowData, dontNeed = read_nenzfr_outfile('./'+low+'/'+opt.exitStatsFileName)
        lowSupply = read_estcj_outfile('./'+low+'/'+opt.jobName+'-estcj.dat')
        lowData['supply_T'] = lowSupply['T']
        lowData['supply_h'] = lowSupply['h']
        
        #print low, nominal, high
    
        if opt.gradient == "linear":
            highValues = get_values(highData,exitVar)
            lowValues = get_values(lowData,exitVar)
            
            highX = DictOfCases[high][perturbedVariables.index(var)]
            lowX = DictOfCases[low][perturbedVariables.index(var)]
            nominalX = DictOfCases[nominal][perturbedVariables.index(var)]
            
            highWeighting = (nominalX - lowX)/(highX - lowX)
            lowWeighting = (highX - nominalX)/(highX - lowX)
            
            sensitivity[var] = ( highWeighting*(array(highValues)-array(nominalValues))/\
                                  (highX - nominalX) + \
                               lowWeighting*(array(nominalValues) - array(lowValues))/\
                                  (nominalX - lowX) ) * nominalX/array(nominalValues)
            #sensitivity[var] = (array(highValues)-array(lowValues))/(highX-lowX)*\
            #                    nominalX/array(nominalValues)
            
            #weightA = (nominalX - lowX)/(highX - lowX)
            #print "weightA=",weightA
            #slopeA = (array(highValues)-array(nominalValues))/(highX - nominalX)
            #print "slopeA=",slopeA
            #
            #weightB = (highX - nominalX)/(highX - lowX)
            #print "weightB=",weightB
            #slopeB = (array(nominalValues)-array(lowValues))/(nominalX - lowX)
            #print "slopeB=",slopeB
            #print "" 
            #print "sensitivity_1=",weightA*slopeA + weightB*slopeB
            #print "sensitivity_2=",(array(highValues)-array(lowValues))/(highX-lowX)
            #print "sensitivity_3=", sensitivity[var]
            #print "difference=",weightA*slopeA + weightB*slopeB-(array(highValues)-array(lowValues))/(highX-lowX)
        else:
            # For quadratic curve fit we have additional cases 
            # that need to be loaded
            tooHigh = 'case'+"{0:02}".format(k)+'3'
            tooHighData,dontNeed = \
                  read_nenzfr_outfile('./'+tooHigh+'/'+opt.exitStatsFileName)
            tooLow = 'case'+"{0:02}".format(k)+'4'
            tooLowData,dontNeed = \
                  read_nenzfr_outfile('./'+tooLow+'/'+opt.exitStatsFileName)

    #print sensitivity
   
    # Now write out a file of the sensitivities as both absolute and relative numbers
    write_sensitivity_summary(sensitivity, perturbedVariables, exitVar, nominalData, 1)
    write_sensitivity_summary(sensitivity, perturbedVariables, exitVar, nominalData, 0) 
    
    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Sensitivity:\n Calculate Sensitivity of Shock Tunnel Test Flow Conditions for a varying inputs"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
