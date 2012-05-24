"""
nenzfr_utils.py -- Small utility functions needed by the main program.
"""

import sys, os
import shlex, subprocess, string
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------

def prepare_input_script(substituteDict, jobName):
    """
    Prepare the actual input file for Eilmer3 from a template.
    """
    templateFileName = E3BIN+"/nenzfr_data_files/nozzle.input.template"
    scriptFileName = jobName + ".py"
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open(scriptFileName, 'w')
    fp.write(text)
    fp.close()
    return

def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    print "About to run cmd:", string.join(args)
    # p = subprocess.Popen(args)
    # wait until the subprocess is finished
    # stdoutData, stderrData = p.communicate()
    return subprocess.check_call(args)

def quote(str):
    """
    Put quotes around a string.
    """
    return '"' + str + '"'

#---------------------------------------------------------------
# Following functions are required by:
#   nenzfr_sensitivity.py
#   nenzfr_RSA.py

def read_case_summary(FileToRead=None):
    """
    Reads the file "perturbation_cases.dat" to determine which variables
    have been perturbed and values for the perturbed variables for each
    case consituting the sensitivity calculation.

    :returns: perturbedVariables - a list of variable names
            : DictOfCases - a dictionary with the case names as keys
                            and the values of each perturbed variable.
    """
    if FileToRead is None:
        fp = open('perturbation_cases.dat','r')
    else:
        fp = open(FileToRead,'r')

    varList = fp.readline().strip().split(" ")
    perturbedVariables = [k for k in varList if k!="#" and k!=""]
    fp.readline()
    DictOfCases = {}
    for line in fp.readlines():
        caseData = line.strip().split(" ")
        caseName = caseData[0]
        DictOfCases[caseName] = [float(k) for k in caseData \
                                 if k!=caseName and k!=""]
    fp.close()
    return perturbedVariables, DictOfCases

def read_nenzfr_outfile(FileToRead):
    """
    Reads the nenzfr out-file containing the statistics of the
    exit flow properties.

    :FileToRead: the name of the file to be read. Default should be
                 something like "nozzle-exit.stats"
    :returns: a list of all the exit flow property names  and a
              dictionary of the mean-values.
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
    fp.close()
    return exitDataDict, exitProperty

def read_estcj_outfile(FileToRead):
    """
    Read just the line in the estcj output file that contains
    the nozzle supply temperature and enthalpy.
    """
    supplyDict = {}
    lineToRead = 0

    fp = open(FileToRead,'r')
    for line in fp.readlines():
        # We only want to read one line in the file
        if lineToRead == 1:
            for data in line.split(','):
                # Data is a list of strings with one string for
                # each property. The form of the strings is
                # 'property: value unit'. We split the list and
                # create "values" which is a list of lists of the
                # form
                #     values = [[property],[value, unit]]
                # We then put this into a dictionary.
                values = [x.split() for x in data.split(':')]
                supplyDict[values[0][0]] = float(values[1][0])
            lineToRead = 0

        if line.strip() in \
            ['State 5s: equilibrium condition (relaxation to pe)',]:
            # We want to read the next line
            lineToRead = 1
    fp.close()
    return supplyDict

#---------------------------------------------------------------
# The following is required by:
#    nenzfr_perturbed.py
#    nenzfr_quasi_transient.py

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


