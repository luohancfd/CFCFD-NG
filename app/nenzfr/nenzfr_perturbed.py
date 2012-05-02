#!/usr/bin/env python
# nenzfr_perturbed.py
#
# This script coordinates the running of a series of nenzfr 
# calculations that are perturbations around a nominal condition.
# The results may be used to create either a LUT or calculate
# the sensitivities of the freestream properties to the inputs.
# 
# Luke Doherty
# School of Mechancial and Mining Engineering
# The University of Queensland

VERSION_STRING = "02-May-2012"

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

def prepare_run_script(substituteDict, jobName, Cluster):
    """
    Prepare the actual run file for Nenzfr from a template.
    """
    templateFileName = "run_template_"+Cluster +".sh"
    # Check that the templateFileName is in the current directory,
    # if it isn't copy it from "E3BIN/nenzfr_data_files/"
    if not os.path.exists(templateFileName):
        command_text = 'cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./'
        #print command_text
        run_command('cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./')
    scriptFileName = "run_" + jobName + ".sh"
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.safe_substitute(substituteDict)
    fp = open(scriptFileName, 'w')
    fp.write(text)
    fp.close()
    return scriptFileName
    
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


def configure_input_dictionary(perturbedVariables, perturbedDict,\
                               defaultPerturbations, levels, perturb): 
    """
    Take the input dictionary, for which the values corresponding to each key are 
    strings and convert the strings into a list of floats. Return the adjusted input
    dictionary.
    """
    # If necessary, remove BLTrans variable from the perturbedDict so that 
    # we don't try to make a float of it later.
    if 'BLTrans' in perturbedDict.keys():
        if perturbedDict['BLTrans'] == "x_c[-1]*1.1":
            del perturbedDict['BLTrans']

    for i in range(len(perturbedDict)):
        var = perturbedDict.keys()[i]
        temp =  perturbedDict.values()[i]
        temp = temp.strip(']').strip('[').split(',')
        temp = [float(temp[j]) for j in range(len(temp))]
        if var in perturbedVariables:
            perturbedDict[var] = set_perturbed_values(var, temp, \
                                defaultPerturbations, perturb, levels)
        else:
            perturbedDict[var] = [temp[0]]
    return perturbedDict

def set_perturbed_values(var, temp, defaultPerturbations, perturb, levels):
    """
    Caculate a list of perturbed values for variable of interest based on 
    input data.
    """
    if levels == 3:
        if len(temp) == 1:
            # have just the value
            print "Using default relative perturbation of "+\
                   str(defaultPerturbations[var])+"% for "+var
            temp = [temp[0], temp[0]*(1+defaultPerturbations[var]/100.0),\
                             temp[0]*(1-defaultPerturbations[var]/100.0)]
        elif len(temp) == 2:
            # have either [value, delta]         or
            #             [value, percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1]]
        elif len(temp) == 3:
            # have either [value, +delta, -delta]                  or 
            #             [value, +percent_delta, -percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0+temp[2]/100.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]+temp[2]]
        else:
            print "Too many values given for "+var
            return -2
    elif levels == 5:
        if len(temp) == 1:
            print "Using default relative perturbation of "+\
                           str(defaultPerturbations[var])+"% for "+var
            temp = [temp[0], temp[0]*(1+defaultPerturbations[var]/100.0),\
                             temp[0]*(1-defaultPerturbations[var]/100.0),\
                             temp[0]*(1+defaultPerturbations[var]/100.0*2.0),\
                             temp[0]*(1-defaultPerturbations[var]/100.0*2.0)]
        elif len(temp) == 2:
            # have either [value, delta]            or
            #              [value, +percent_delta]
            # we assume that "2"*delta == 2*delta
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0),\
                                 temp[0]*(1.0+temp[1]/100.0*2.0),\
                                 temp[0]*(1.0-temp[1]/100.0*2.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1],\
                                 temp[0]+2.0*temp[1],\
                                 temp[0]-2.0*temp[1]]
        elif len(temp) == 3:
            # have either [value, delta, "2"*delta]                  or
            #             [value, percent_delta, "2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0),\
                                 temp[0]*(1.0+temp[2]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1],\
                                 temp[0]+temp[2],\
                                 temp[0]-temp[2]]
        elif len(temp) == 4:
            # have either [value, +delta, -delta, "2"*delta]                          or
            #             [value, +percent_delta, -percent_delta, "2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0),\
                                 temp[0]*(1.0+temp[3]/100.0),\
                                 temp[0]*(1.0-temp[3]/100.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[2],\
                                 temp[0]+temp[3],\
                                 temp[0]-temp[3]]
        elif len(temp) == 5:
            # have either [value, +delta, -delta, +"2"*delta, -"2"*delta] or
            #             [value +percent_delta, -percent_delta, +"2"*percent_delta, -"2"*percent_delta]
            if perturb == "rel":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0),\
                                 temp[0]*(1.0+temp[3]/100.0),\
                                 temp[0]*(1.0-temp[4]/100.0)]
            else:
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[2],\
                                 temp[0]+temp[3],\
                                 temp[0]-temp[4]]
        else:
            print "Too many values given for "+var
            return -2
    return temp

def set_case_running(caseString, caseDict, textString):
    """
    A short function to set a given case running in an appropriately
    named sub-directory.
    """
    print 60*"-"
    print caseString
    print textString
    
    # Create sub-directory for the current case
    run_command('mkdir ./'+caseString)
     
    # Set up the run script for Nenzfr
    scriptFileName = prepare_run_script(caseDict, \
        caseDict['jobName'].strip('"')+'_'+caseString, caseDict['Cluster'])
    
    # Move the run script to its sub-directory
    command_text = 'mv '+scriptFileName+' ./'+caseString+'/'+scriptFileName
    run_command(command_text)

    # If required, copy the nozzle.timing file to the sub-directory
    if caseDict['blockMarching'] in ["--block-marching",]:
        if os.path.exists('nozzle.timing'):
            command_text = 'cp nozzle.timing ./'+caseString+'/'
            run_command(command_text)
    
    # If require, copy the equilibrium gas LUT to the sub-directory
    if caseDict['chemModel'] in ['"eq"',]:
        command_text = 'cp ./'+caseDict['gmodelFile']+' ./'+caseString+'/'
        run_command(command_text)
    
    # Change into the sub-directory, ensure the run script is exectuable and
    # then run it
    os.chdir(caseString)
    run_command('chmod u+x '+scriptFileName)
    print ""
    print caseDict['runCMD']+scriptFileName
    print ""
    # I am not sure how to replace the following line with the run_command function
    os.system(caseDict['runCMD']+scriptFileName)
    os.chdir('../')
    return

def write_case_summary(varList,caseDict,caseString,newfile):
    """
    A short function to write the values of the perturbed variables 
    for the current case to a summary file.
    """
    # Define some dictionaries with the formating of each of the possible 
    # perturbed variables
    formatDict = {'p1':'{0:>13.5g}', 'T1':'{0:>10.5g}', 'Vs':'{0:>11.5g}',
                  'pe':'{0:>15.6g}', 'Tw':'{0:>10.5g}', 'BLTrans':'{0:>10.5g}',
                  'TurbVisRatio':'{0:>14.5g}', 'TurbInten':'{0:>11.5g}',
                  'CoreRadiusFraction':'{0:>20.5g}'}
    titleFormatDict = {'p1':'{0:{fill}>13}', 'T1':'{0:{fill}>10}', 
                       'Vs':'{0:{fill}>11}', 'pe':'{0:{fill}>15}', 
                       'Tw':'{0:{fill}>10}', 'BLTrans':'{0:{fill}>10}',
                       'TurbVisRatio':'{0:{fill}>14}', 'TurbInten':'{0:{fill}>11}',
                       'CoreRadiusFraction':'{0:{fill}>20}'}
    # For the first time we create a new file and write the header information.
    # Each following case is appended to the existing file.
    if newfile == 1:
        fout = open('perturbation_cases.dat','w')
        # Write title line
        fout.write('{0:>7}'.format('#'))
        for k in tuple(varList):
            fout.write(titleFormatDict[k].format(k,fill=''))
        fout.write('\n')
        # Underline the title
        for k in varList:
            fout.write(titleFormatDict[k].format('-',fill='-'))
        fout.write('{0:->7}'.format('-'))
        fout.write('\n')
    else:
        fout = open('perturbation_cases.dat','a')
    # Now write out the data for the current case
    fout.write('{0:>7}'.format(caseString))
    for k in varList:
        fout.write(formatDict[k].format(caseDict[k]))
    fout.write('\n')
    fout.close()

def write_case_config(caseDict):
    """
    """
    fout = open('nominal_case.config','w')
    for k in range(len(caseDict)):
        fout.write('{0}'.format(caseDict.keys()[k]+':   '))
        fout.write('{0}'.format(caseDict.values()[k]))
        fout.write('\n')
    fout.close()    

def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate a series of Nenzfr calculations with 
    various inputs perturbed around the input nominal condition.
    """
    op = optparse.OptionParser(version=VERSION_STRING)

    op.add_option('--runCMD', dest='runCMD', default='./',
                  help=("command used to execute the run-script file "
                        "[default: %default]"))
    op.add_option('--cluster', dest='Cluster', default='Mango',
                  choices =['Mango', 'Barrine'],
                  help=("specify on which cluster the computations are to be ran. "
                        "This is used to define which run template script will "
                        "be used. Options: "
                        "Mango; Barrine [default: %default]"))
    op.add_option('--perturb', dest='perturb', default='rel', choices=['rel','abs'],
                  help=("specify whether the given perturbations are absolute "
                        "or relative values. Options: rel, abs [default: %default]"))
    op.add_option('--levels', dest='levels', default=3,choices=['3','3-reduced','5'],
                  help=("specify how many values are to be used for each variable, "
                        "including the nominal. Options: 3, 3-reduced, 5 "
                        "[default: %default]"))
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")

    # Multitude of options required by nenzfr.
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"))

    group = optparse.OptionGroup(op, "Inputs that may be perturbed:", 
                  "The following inputs my be perturbed by providing a list in one of "
                  "the following forms: "
                  " (1)  [nominal, +delta];"
                  " (2)  [nominal, +delta_1, -delta_2];"
                  " (3)  [nominal, +delta_1, -delta_2, +delta_3];"
                  " (4)  [nominal, +delta_1, -delta_2, +delta_3, -delta_4]. "
                  "Note that the deltas may be given as absolute values or "
                  "as relative values in the form of percentages.")
    group.add_option('--p1', dest='p1', default=None, 
                  help=("shock tube fill pressure (in Pa) and its perturbation/s "
                        "as a list. [default delta: 5%]" ))
    group.add_option('--T1', dest='T1', default=None,
                  help=("shock tube fill temperature, in degrees K and its perturbation/s "
                        "as a list. [default delta: 5%]"))
    group.add_option('--Vs', dest='Vs', default=None,
                  help=("incident shock speed, in m/s and its perturbation/s as a list. "
                        "[default delta: 5%]"))
    group.add_option('--pe', dest='pe', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa and its "
                        "perturbation/s as a list. [default delta: 5%]"))

    group.add_option('--Twall', dest='Tw', default='300.0',
                  help=("Nozzle wall temperature, in K and its perturbation/s as a list. "
                        "[default: %default, delta: 5%]"))
    group.add_option('--BLTrans', dest='BLTrans', default="x_c[-1]*1.1",
                  help=("Transition location for the Boundary layer and its perturbation/s "
                        "as a list. Used to define the turbulent portion of the nozzle. "
                        "[default: >nozzle length i.e. laminar nozzle, delta: 5%]"))
    group.add_option('--TurbVisRatio', dest='TurbVisRatio', default='100.0',
                  help=("Turbulent to Laminar Viscosity Ratio and its perturbation/s "
                        "as a list. [default: %default, delta: 5%]"))
    group.add_option('--TurbIntensity', dest='TurbInten', default='0.05',
                  help=("Turbulence intensity at the throat and its perturbation/s "
                  "[default: %default, delta: 5%]"))
    group.add_option('--CoreRadiusFraction', dest="coreRfraction", default='0.6666666667',
                  help=("Radius of core flow as a fraction of "
                  "the nozzle exit radius and its perturbation/s "
                  "[default: %default, delta: 5%]"))
    op.add_option_group(group)
    
    op.add_option('--chem', dest='chemModel', default='eq',
                  choices=['eq', 'neq', 'frz', 'frz2'],
                  help=("chemistry model: " "eq=equilibrium; "
                        "neq=non-equilibrium; " "frz=frozen "
                        "[default: %default]"))
    op.add_option('--area', dest='areaRatio', default=1581.165,
                  help=("nozzle area ratio. only used for estcj calc. "
                        "use when --cfile(--gfile) are "
                        "specified. [default: %default]"))
    op.add_option('--cfile', dest='contourFileName',
                  default='contour-t4-m10.data',
                  help="file containing nozzle contour [default: %default]")
    op.add_option('--gfile', dest='gridFileName', default='None',
                  help="file containing nozzle grid. "
                  "overrides --cfile if both are given "
                  "[default: %default]")
    op.add_option('--exitfile', dest='exitSliceFileName',
                  default='nozzle-exit.data',
                  help="file for holding the nozzle-exit data [default: %default]")
    op.add_option('--block-marching', dest='blockMarching', action='store_true',
                  default=False, help="run nenzfr in block-marching mode")    
    op.add_option('--create-RSA', dest='createRSA', action='store_true',
                  default=False, 
                  help="perturb only Vs and pe in preparation for calculating "
                       "a Response Surface Approximation")

    # The following defaults suit a Mach 10 Nozzle calculation.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells"))
    op.add_option('--nnj', dest='nnj', type='int', default=100,
                  help=("number of radial cells"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk)"))
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help=("number of radial blocks"))
    op.add_option('--bx', dest='bx', type='float', default=1.05,
                  help=("clustering in the axial direction"))
    op.add_option('--by', dest='by', type='float', default=1.002,
                  help=("clustering in the radial direction"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow"))
    op.add_option('--max-step', dest='max_step', type='int', default=80000,
                  help=("maximum simulation steps allowed"))

    opt, args = op.parse_args()
        
    # Set the default relative perturbation values (as percentages)
    defaultPerturbations = {'p1':2.5, 'T1':2.5, 'Vs':2.5, 'pe':2.5, 
                           'Tw':2.5, 'BLTrans':2.5, 'TurbVisRatio':2.5,
                           'TurbInten':2.5, 'CoreRadiusFraction':2.5}

    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    if opt.p1 is None:
        print "Need to supply a value for p1."
        bad_input = True
    if opt.T1 is None:
        print "Need to supply a value for T1."
        bad_input = True
    if opt.Vs is None:
        print "Need to supply a value for Vs."
        bad_input = True    
    if opt.pe is None:
        print "Need to supply a value for pe."
        bad_input = True
    if opt.blockMarching is True:
        opt.blockMarching = "--block-marching"
    else:
        opt.blockMarching = ""
    if bad_input:
        return -2
    # Convert to integer. NB. I don't specify the type as int in the option
    # as this negates the ability to have choices.
    if opt.levels in ['3-reduced']:
        if opt.createRSA is True:
            opt.levels = 2.5
        else:
            print "the --levels=3-reduced is only valid when creating a response surface."+
                  " changing to levels=3" 
            opt.levels = 3
    else:
       opt.levels = int(opt.levels)
    
    # Set up a list of the perturbed variables.
    if opt.createRSA:
        print "Create Response Surface (RS) option selected. Perturbing only Vs and pe"
        perturbedVariables = ['Vs','pe']
    else:
        perturbedVariables = ['p1','T1','Vs','pe','Tw','BLTrans','TurbVisRatio',
                              'TurbInten','CoreRadiusFraction',]
        if opt.BLTrans == "x_c[-1]*1.1": # Laminar nozzle
            # TODO: Remove TurbInten once it has been shown not to an effect.
            # For some reason it currently (as of 25/04/2012) does.
            print "Nozzle is assumed laminar (based on input BLTrans), therefore\n"+\
                  "BLTrans and TurbVisRatio are not perturbed"
            del perturbedVariables[perturbedVariables.index('BLTrans')]
            del perturbedVariables[perturbedVariables.index('TurbVisRatio')]
            #del perturbedVariables[perturbedVariables.index('TurbInten')]
            
    # Now set up a dictionary of all variables that could possibly be perturbed. 
    # Then configure this dictionary. For each dictionary key, if the key is not
    # in the perturbedVariables list, its value is a list with a single float 
    # element, otherwise it is a list of floats corresponding to the 
    # perturbed values.
    perturbedDict =  {'p1':opt.p1, 'T1':opt.T1, 'Vs':opt.Vs, 'pe':opt.pe,
                      'Tw':opt.Tw, 'BLTrans':opt.BLTrans,
                      'TurbVisRatio':opt.TurbVisRatio,
                      'TurbInten':opt.TurbInten,
                      'CoreRadiusFraction':opt.coreRfraction}
    perturbedDict = configure_input_dictionary(perturbedVariables, perturbedDict,\
                             defaultPerturbations, opt.levels, opt.perturb)
     
    # Set up a parameter dictionary 
    # (with all variables at their nominal value)
    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':opt.pe,
                 'chemModel':quote(opt.chemModel),
                 'contourFileName':quote(opt.contourFileName),
                 'gridFileName':quote(opt.gridFileName),
                 'exitSliceFileName':quote(opt.exitSliceFileName),
                 'areaRatio':opt.areaRatio, 'blockMarching':opt.blockMarching,
                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'nbj':opt.nbj,
                 'bx':opt.bx, 'by':opt.by,
                 'max_time':opt.max_time, 'max_step':opt.max_step,
                 'Tw':opt.Tw, 'TurbVisRatio':opt.TurbVisRatio,
                 'TurbInten':opt.TurbInten, 'BLTrans':opt.BLTrans,
                 'CoreRadiusFraction':opt.coreRfraction,
                 'Cluster':opt.Cluster,'runCMD':opt.runCMD}
    for var in perturbedDict.keys():
        paramDict[var] = perturbedDict[var][0]
        
    # As building an equilibrium gas LUT is so time consuming, we do it here
    # and then copy the resulting LUT into each case sub-directory. The following
    # lines are copied almost verbatim from "nenzfr.py"
    if opt.chemModel in ['eq']:
        if opt.gasName in ['n2']:
            eqGasModelFile = 'cea-lut-'+upper(opt.gasName)+'.lua.gz'
        else:
            eqGasModelFile = 'cea-lut-'+opt.gasName+'.lua.gz'
        if not os.path.exists(eqGasModelFile):
            run_command('build-cea-lut.py --gas='+opt.gasName)
        paramDict['gmodelFile'] = eqGasModelFile

    if opt.createRSA is not True: # Perturbing for a sensitivity calculation
        # Calculate Nominal condition
        caseString = 'case'+"{0:02}".format(0)+"{0:01}".format(0)
        textString = "Nominal Condition"
        caseDict = copy.copy(paramDict)
        # Run the nominal case and write the values of the perturbed variables
        # to a summary file
        set_case_running(caseString, caseDict, textString)
        write_case_summary(perturbedVariables,caseDict,caseString,1)

        # Now run all the perturbed conditions
        for k in range(len(perturbedVariables)):
            var = perturbedVariables[k]
            perturbCount = opt.levels

            for kk in range(perturbCount):
                if kk != 0:
                    caseString = 'case'+"{0:02}".format(k)+"{0:01}".format(kk)
                    textString = var+" perturbed to "+str(perturbedDict[var][kk])
                    # Perturbation of the "CoreRadiusFraction" input may be done
                    # by (re)post processing the nominal case solution using the
                    # --just-stats option in nenzfr. We therefore don't need to
                    # run any separate cases for this variable. The perturbation
                    # is handled by "nenzfr_sensitivity.py".
                    caseDict = copy.copy(paramDict)
                    caseDict[var] = perturbedDict[var][kk]
                    if var != 'CoreRadiusFraction':
                        # Run the current case
                        set_case_running(caseString, caseDict, textString)
                    # Write current case to the summary file
                    write_case_summary(perturbedVariables,caseDict,\
                                       caseString,0)

    else: # Perturbing to create a LUT
        var1 = perturbedVariables[0] # 'Vs'
        var2 = perturbedVariables[1] # 'pe'
        
        if opt.levels == 2.5
            casesToRun = [(2,1),      (1,1),
                                (0,0),
                          (2,2),      (1,2)]
        elif opt.levels == 3:
            casesToRun = [(2,1),(0,1),(1,1),
                          (2,0),(0,0),(1,0),
                          (2,2),(0,2),(1,2)]
        elif opt.levels == 5:
            casesToRun = [(4,3),      (0,3),      (3,3),
                                (2,1),      (1,1),      
                          (4,0),      (0,0),      (3,0),
                                (2,2),      (1,2),      
                          (4,4),      (0,4),      (3,4)]
            #casesToRun = [      (2,3),      (1,3),     
            #              (4,1),      (0,1),      (3,1),
            #                    (2,0),      (1,0),     
            #              (4,2),      (0,2),      (3,2),
            #                    (2,4),      (1,4)      ]
        
        # Run the nominal case first
        caseString = 'case'+"{0:01}{0:01}".format(0,0)
        caseDict = copy.copy(paramDict)
        textString = "Nominal Case: "+var1+"="+str(perturbedDict[var1][0])+\
                                 "; "+var2+"="+str(perturbedDict[var2][1])
        write_case_config(caseDict)
        set_case_running(caseString, caseDict, textString)
        write_case_summary(perturbedVariables, caseDict, caseString, 1)
         
        # Now run all other cases
        for case in casesToRun:
            if case != (0,0):
                caseString = 'case'+"{0:01}{1:01}".format(case[0],case[1])
                textString = var1+" perturbed to "+str(perturbedDict[var1][case[0]])+                                "\n"+var2+" perturbed to "+str(perturbedDict[var2][case[1]])
                caseDict = copy.copy(paramDict)
                caseDict[var1] = perturbedDict[var1][case[0]]
                caseDict[var2] = perturbedDict[var2][case[1]]

                set_case_running(caseString, caseDict, textString)
                write_case_summary(perturbedVariables, caseDict, caseString, 0)

    
    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Perturbed:\n Calculate Shock Tunnel Test Flow Conditions for inputs perturbed about the nominal"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
