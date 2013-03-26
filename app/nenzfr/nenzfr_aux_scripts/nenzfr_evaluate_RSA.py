#!/usr/bin/env python
# nenzfr_evaluate_RSA.py
# 
# Used when evaluating the Mach 4 nozzle simulations.
# Script needs a clean up/commenting corrected and added.
# 
# Luke Doherty
# School of Mechancial and Mining Engineering
# The University of Queensland
#
# Version: 08-May-2012 - Initial completed and working code
#          04-Feb-2013 - Added --RSAfile and --add-extra-variables
#                        options. 
#          05-Feb-2013 - Added Root-Sum-Square to summary file 
VERSION_STRING = "04-Feb-2013"

import shlex, subprocess, string
from subprocess import PIPE
import sys, os, gzip
import optparse
#from numpy import array, mean, logical_and, zeros, dot, sqrt, linalg
from numpy import *
import copy
from nenzfr_sensitivity import read_case_summary
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

def read_outfile(FileToRead):
    """
    Reads the out-file containing the residuals of the
    exit flow properties.

    :FileToRead: the name of the file to be read. Default should be
                 something like "case00-exit.RSAdat_residuals"
    :returns: a list of all the exit flow property names and a 
              dictionary of the mean-values.
    """
    fp = open(FileToRead,'r')
    fp.readline() # Names. Don't need.
    fp.readline() # Vs row. Already have this data
    fp.readline() # pe row. Already have this data
    fp.readline() # This is just a row of "-"
    
    residualDict = {}
    fileLines = fp.readlines()
    exitProperty = []
    # Now read in the rest of the data
    for line in fileLines:
        data = line.strip().split(" ")
        values = [k for k in data if k!=""]
        variable = values[0]
        exitProperty.append(variable)
        residualDict[variable] = float(values[1])
    fp.close()
    return residualDict, exitProperty

def write_summaryFile(residuals, exitVar, DictOfCases, casesList=None):
    """
    """
    fout = open('residual_summary.dat','w')
    
    if casesList is None:
        casesList = DictOfCases.keys()
    if 'case00' in casesList: 
        del casesList[casesList.index('case00')]
     
    # Write out title line
    fout.write('{0:>12}'.format('variable'))
    for case in casesList:
        fout.write('{0:>15}'.format(case))
    fout.write('{0:>15}'.format('Root-Sum-'))
    fout.write('\n')
    # Write out the Vs values for each case
    fout.write('{0:>12}'.format('Vs'))
    for case in casesList:
        fout.write('{0:>15.5g}'.format(DictOfCases[case][0]))
    fout.write('{0:>15}'.format('Square'))
    fout.write('\n')
    # Write out the pe values for each case
    fout.write('{0:>12}'.format('pe'))
    for case in casesList:
        fout.write('{0:>15.6g}'.format(DictOfCases[case][1]))
    fout.write('{0:>15}'.format(''))
    fout.write('\n')
    # Write out a horizontal line
    for k in range(len(DictOfCases)):
        fout.write('{0:->15}'.format('-'))
    fout.write('{0:->27}'.format('-'))
    fout.write('\n')
    
    # Write out data for each case and exit variable. We
    # also calculate and write out the root-sum-square
    # residual
    for var in exitVar:
        fout.write('{0:>12}'.format(var))
        sumSquare = 0.0
        for case in casesList:
            sumSquare += float(square(residuals[case][var]))
            fout.write('{0:>15.7g}'.format(residuals[case][var]))
        fout.write('{0:>15.7g}'.format(sqrt(sumSquare)))
        fout.write('\n')
    fout.close()
    return 0
        
def main():
    """
    Examine the command-line options to decide what to do.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    
    op.add_option('--run-defaults', dest='runDefaults', action='store_true',
                  default=True, help="[default: %default]")
    op.add_option('--RSA-file', dest='RSAfile', default='response-surface.dat',
                  help="specify the Response-Surface file name "
                  "[default: %default]")
    op.add_option('--perturb-file', dest='perturbFile', default='perturbation_cases.dat',
                  help="specify the file name holding the perturbation case information "
                  "[default: %default]")
    op.add_option('--save-dir', dest='saveDir', default=None,
                  help="specify the directory to which results should be saved "
                  "[default: %default]")
    op.add_option('--write-summary', dest='writeSummary', action='store_true',
                  default=False, help="[default: %default]")
    op.add_option('--add-extra-variables', dest='addExtraVariables', action='store_true',
                  default=False, help=("specify whether q, rho*u_x, Re_u and p/q should "
                       "also be calculated. [default: %default]"))
    opt, args = op.parse_args()
    
    # Read the perturbation summary  file to get the perturbed
    # variables and their various values
    perturbedVariables, DictOfCases = read_case_summary(opt.perturbFile)
    
    if os.path.exists('./case00') or os.path.exists('./caseA'):
        runCMD = './'
    elif os.path.exists('../case00') or os.path.exists('../caseA'):
        runCMD = '../'
    else:
        print "Cannot find case folders"
        return -2
     
    # Loop through each case
    for case in DictOfCases.keys():
        if case not in ['case00',]:
            #command_text = "nenzfr_RSA.py --Vs="+str(DictOfCases[case][0])+" --pe="+\
            #           str(DictOfCases[case][1])+" --exitFile="+case+"-exit.RSAdat"+\
            #           " --calculate-residuals --estcjFile="+runCMD+case+"/nozzle-estcj.dat"+\
            #           " --exitStatsFile="+runCMD+case+"/nozzle-exit.stats"
            if opt.addExtraVariables is True:
                addExtraVariables = " --add-extra-variables"
            else:
                addExtraVariables = ""
            
            command_text = "nenzfr_RSA.py --Vs="+str(DictOfCases[case][0])+" --pe="+\
                       str(DictOfCases[case][1])+" --exitFile="+case+"-exit.RSAdat"+\
                       " --calculate-residuals --estcjFile="+runCMD+case+"/nozzle-estcj.dat"+\
                       " --exitStatsFile="+runCMD+case+"/nozzle-exit.stats"+\
                       " --RSA-file="+opt.RSAfile+addExtraVariables
            
            #print command_text
            #run_command(command_text)
            print ""
            print "Running case: "+case
            print ""
            run_command(command_text)
            
    for case in DictOfCases.keys():
        if case not in ['case00']:
            if opt.saveDir is not None:
                command_text = "mv "+case+"-exit.RSAdat* "+opt.saveDir
                #print command_text
                os.system(command_text)
                #run_command(command_text)
    
    if opt.writeSummary:
        residuals = {}
        #caseList = [         'case23',         'case13',\
        #            'case41',         'case01',         'case31',\
        #                     'case20',         'case10',\
        #            'case42',         'case02',         'case32',\
        #                     'case24',         'case14',\
        #            'caseA','caseB','caseC','caseD']
        caseList = ['caseA','caseB','caseC','caseD','caseE','caseF']
        for case in DictOfCases.keys():
            if case not in ['case00']:
                residuals[case], exitProperty = read_outfile(case+'-exit.RSAdat_residuals')
        write_summaryFile(residuals, exitProperty, DictOfCases, caseList)
    return 0

#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr Evaluate RSA"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
