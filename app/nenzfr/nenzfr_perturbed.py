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

VERSION_STRING = "24-May-2012"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
from numpy import array, mean, logical_and, zeros
from nenzfr_utils import run_command, quote, prepare_run_script
from nenzfr_input_utils import input_checker, nenzfr_perturbed_input_checker
import copy
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------
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
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0)]
            elif perturb == 'absolute':
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1]]
        elif len(temp) == 3:
            # have either [value, +delta, -delta]                  or 
            #             [value, +percent_delta, -percent_delta]
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0+temp[2]/100.0)]
            elif perturb == 'absolute':
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
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0),\
                                 temp[0]*(1.0+temp[1]/100.0*2.0),\
                                 temp[0]*(1.0-temp[1]/100.0*2.0)]
            elif perturb == 'absolute':
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1],\
                                 temp[0]+2.0*temp[1],\
                                 temp[0]-2.0*temp[1]]
        elif len(temp) == 3:
            # have either [value, delta, "2"*delta]                  or
            #             [value, percent_delta, "2"*percent_delta]
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[1]/100.0),\
                                 temp[0]*(1.0+temp[2]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0)]
            elif perturb == 'absolute':
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[1],\
                                 temp[0]+temp[2],\
                                 temp[0]-temp[2]]
        elif len(temp) == 4:
            # have either [value, +delta, -delta, "2"*delta]                          or
            #             [value, +percent_delta, -percent_delta, "2"*percent_delta]
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0),\
                                 temp[0]*(1.0+temp[3]/100.0),\
                                 temp[0]*(1.0-temp[3]/100.0)]
            elif perturb == 'absolute':
                temp = [temp[0], temp[0]+temp[1],\
                                 temp[0]-temp[2],\
                                 temp[0]+temp[3],\
                                 temp[0]-temp[3]]
        elif len(temp) == 5:
            # have either [value, +delta, -delta, +"2"*delta, -"2"*delta] or
            #             [value +percent_delta, -percent_delta, +"2"*percent_delta, -"2"*percent_delta]
            if perturb == "relative":
                temp = [temp[0], temp[0]*(1.0+temp[1]/100.0),\
                                 temp[0]*(1.0-temp[2]/100.0),\
                                 temp[0]*(1.0+temp[3]/100.0),\
                                 temp[0]*(1.0-temp[4]/100.0)]
            elif perturb == 'absolute':
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
    scriptFileName, cfgFileName = prepare_run_script(caseDict, \
        caseDict['jobName'].strip('"')+'_'+caseString, caseDict['Cluster'])
    
    # Move the run script to its sub-directory
    command_text = 'mv '+scriptFileName+' ./'+caseString+'/'+scriptFileName
    run_command(command_text)
    
    # Move the cfg file to its sub-directory too
    command_text = 'mv '+cfgFileName+' ./'+caseString+'/'+cfgFileName
    run_command(command_text)
    
    # If require, copy the equilibrium gas LUT to the sub-directory
    if caseDict['chemModel'] == 'eq':
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

def main(cfg={}):
    """
    Examine the command-line options to decide the what to do
    and then coordinate a series of Nenzfr calculations with 
    various inputs perturbed around the input nominal condition.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename for the config file"))
    opt, args = op.parse_args()
    config_file = opt.config_file
       
    if not cfg: #if the configuration dictionary has not been filled up already, load it from a file
        try: #from Rowan's onedval program
            execfile(config_file, globals(), cfg)
        except IOError as e:
            print "Error {0}".format(str(e))
            print "There was a problem reading the config file: '{0}'".format(config_file)
            print "Check that it conforms to Python syntax."
            print "Bailing out!"
            sys.exit(1)
            
    #check inputs using original nenzfr input checker first
 
    cfg['bad_input'] = False
   
    cfg = input_checker(cfg)
    
    # add default pertubations and check new nenzfr perturbed inputs
    
    # Set the default relative perturbation values (as percentages)
    cfg['defaultPerturbations'] = {'p1':2.5, 'T1':2.5, 'Vs':2.5, 'pe':2.5, 
                           'Tw':2.5, 'BLTrans':2.5, 'TurbVisRatio':2.5,
                           'TurbInten':2.5, 'CoreRadiusFraction':2.5}
    
    cfg = nenzfr_perturbed_input_checker(cfg)
    
    #bail out here if there is an issue
    if cfg['bad_input']:
        return -2
        
    perturbedDict = cfg['perturbedDict']
        
    for var in perturbedDict.keys():
        cfg[var] = perturbedDict[var][0]
        
    # As building an equilibrium gas LUT is so time consuming, we do it here
    # and then copy the resulting LUT into each case sub-directory. The following
    # lines are copied almost verbatim from "nenzfr.py"
    if cfg['chemModel'] in ['eq']:
        if cfg['gasName'] in ['n2']:
            eqGasModelFile = 'cea-lut-'+upper(cfg['gasName'])+'.lua.gz'
        else:
            eqGasModelFile = 'cea-lut-'+cfg['gasName']+'.lua.gz'
        if not os.path.exists(eqGasModelFile):
            run_command('build-cea-lut.py --gas='+cfg['gasName'])
        cfg['gmodelFile'] = eqGasModelFile

    if not cfg['createRSA']: # Perturbing for a sensitivity calculation
        # Calculate Nominal condition
        caseString = 'case'+"{0:02}".format(0)+"{0:01}".format(0)
        textString = "Nominal Condition"
        caseDict = copy.copy(cfg)
        caseDict['caseName'] = caseString
        write_case_config(caseDict)
        # Run the nominal case and write the values of the perturbed variables
        # to a summary file
        set_case_running(caseString, caseDict, textString)
        write_case_summary(cfg['perturbedVariables'],caseDict,caseString,1)

        # Now run all the perturbed conditions
        for k in range(len(cfg['perturbedVariables'])):
            var = cfg['perturbedVariables'][k]
            perturbCount = cfg['levels']

            for kk in range(perturbCount):
                if kk != 0:
                    caseString = 'case'+"{0:02}".format(k)+"{0:01}".format(kk)
                    textString = var+" perturbed to "+str(perturbedDict[var][kk])
                    # Perturbation of the "CoreRadiusFraction" input may be done
                    # by (re)post processing the nominal case solution using the
                    # --just-stats option in nenzfr. We therefore don't need to
                    # run any separate cases for this variable. The perturbation
                    # is handled by "nenzfr_sensitivity.py".
                    caseDict = copy.copy(cfg)
                    caseDict[var] = perturbedDict[var][kk]
                    caseDict['caseName'] = caseString
                    if var != 'CoreRadiusFraction':
                        # Run the current case
                        set_case_running(caseString, caseDict, textString)
                    # Write current case to the summary file
                    write_case_summary(cfg['perturbedVariables'],caseDict,\
                                       caseString,0)

    else: # Perturbing to create a LUT
        var1 = cfg['perturbedVariables'][0] # 'Vs'
        var2 = cfg['perturbedVariables'][1] # 'pe'
        
        if cfg['levels'] == 2.5:
            casesToRun = [(2,1),      (1,1),
                                (0,0),
                          (2,2),      (1,2)]
        elif cfg['levels'] == 3:
            casesToRun = [(2,1),(0,1),(1,1),
                          (2,0),(0,0),(1,0),
                          (2,2),(0,2),(1,2)]
        elif cfg['levels'] == 5:
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
        caseDict['caseName'] = caseString
        textString = "Nominal Case: "+var1+"="+str(perturbedDict[var1][0])+\
                                 "; "+var2+"="+str(perturbedDict[var2][0])
        write_case_config(caseDict)
        set_case_running(caseString, caseDict, textString)
        write_case_summary(cfg['perturbedVariables'], caseDict, caseString, 1)
         
        # Now run all other cases
        for case in casesToRun:
            if case != (0,0):
                caseString = 'case'+"{0:01}{1:01}".format(case[0],case[1])
                textString = var1+" perturbed to "+str(perturbedDict[var1][case[0]])+\
                        "\n"+var2+" perturbed to "+str(perturbedDict[var2][case[1]])
                caseDict = copy.copy(paramDict)
                caseDict['caseName'] = caseString
                caseDict[var1] = perturbedDict[var1][case[0]]
                caseDict[var2] = perturbedDict[var2][case[1]]

                set_case_running(caseString, caseDict, textString)
                write_case_summary(cfg['perturbedVariables'], caseDict, caseString, 0)

    
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
