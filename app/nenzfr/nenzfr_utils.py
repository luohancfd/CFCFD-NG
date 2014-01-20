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
       
    # Flush before using subprocess to ensure output is in the right order.
    sys.stdout.flush()    
    
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

def read_gmodelFile_from_config(jobName):
    """
    Read the jobName.config file and return the name of 
    the gas model that was used
    """
    fp = open(jobName+'.config','r')
    for line in fp.readlines():
        if line.startswith('gas_model_file'):
            data = line.strip().split()
            gas_model_file = str(data[-1])
            break
    return gas_model_file

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

def read_nenzfr_outfile(FileToRead,inclStats=0):
    """
    Reads the nenzfr out-file containing the statistics of the
    exit flow properties.

    :FileToRead: the name of the file to be read. Default should be
                 something like "nozzle-exit.stats"
    :inclStats: specify whether the min, max, and standard deviation
                should also be returned
    :returns: a list of all the exit flow property names  and a
              dictionary of the mean-values.
    """
    fp = open(FileToRead,'r')
    fp.readline() # first line of file (CoreRadiusFraction which we don't use)
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
        if inclStats == 0: # Return just the mean value (default)
            exitDataDict[variable] = float(values[1])
        elif inclStats == 1: # Return a dictionary of all values
            exitDataDict[variable] = {'mean':float(values[1]),
                                      'min':float(values[2]),
                                      'max':float(values[3]),
                                      'std':float(values[4])}
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
   
def prepare_run_script(configDict, jobName, Cluster):
    """
    Prepare the actual run file for Nenzfr from a template.
    
    Also builds you a cool config file to go with it.
    """

    templateFileName = "run_template_"+Cluster +".sh"
    # Check that the templateFileName is in the current directory,
    # if it isn't copy it from "E3BIN/nenzfr_data_files/"
    if not os.path.exists(templateFileName):
        command_text = 'cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./'
        #print command_text
        run_command('cp '+E3BIN+'/nenzfr_data_files/'+templateFileName+' ./')
    scriptFileName = "run_" + jobName + ".sh"
    cfgFileName = jobName + ".cfg"
    substituteDict = {'config_filename': cfgFileName, 
                      'caseName': configDict['caseName']}
    
    #fill up the shell script
    fp = open(templateFileName, 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.safe_substitute(substituteDict)
    fp = open(scriptFileName, 'w')
    fp.write(text)
    fp.close()
    
    #now fill up the cfg file
    
    build_cfg_file(configDict, cfgFileName)
    
    return scriptFileName, cfgFileName

def build_cfg_file(cfg, filename):
    """
    Takes a nenzfr config dictionary and uses it to build a nenzfr config file.close

    Based on a program I built for Luke to convert old nenzfr optparse inputs
    to the new cfg file style.
    
    Chris James (09/05/13)
    
    """
    
    config_file = open(filename,"w")  #txt_output file creation
    print "Opening file '{0}'.".format(filename)
    
    print "Starting the process of adding inputs to the config file."

    facility_string = "facility = '{0}'".format(cfg['facility'])
    config_file.write(facility_string+'\n')
    
    empty_line_string = " \n"
    
    config_file.write(empty_line_string)
    
    gasName_string = "gasName = '{0}'".format(cfg['gasName'])
    config_file.write(gasName_string+'\n')
    
    config_file.write(empty_line_string)

    chemModel_string = "chemModel = '{0}'".format(cfg['chemModel'])
    config_file.write(chemModel_string+'\n')
    
    config_file.write(empty_line_string)    

    p1_string = "p1 = {0}".format(cfg['p1'])
    config_file.write(p1_string+'\n')    

    T1_string = "T1 = {0}".format(cfg['T1'])
    config_file.write(T1_string+'\n')     

    Vs_string = "Vs = {0}".format(cfg['Vs'])
    config_file.write(Vs_string+'\n')

    pe_string = "pe = {0}".format(cfg['pe'])
    config_file.write(pe_string+'\n')
    
    config_file.write(empty_line_string)   

    areaRatio_string = "areaRatio = {0}".format(cfg['areaRatio'])
    config_file.write(areaRatio_string+'\n')

    config_file.write(empty_line_string)
    
    jobName_string = "jobName = '{0}'".format(cfg['jobName'])
    config_file.write(jobName_string+'\n')

    config_file.write(empty_line_string)
    
    contourFileName_string = "contourFileName = '{0}'".format(cfg['contourFileName'])
    config_file.write(contourFileName_string+'\n')
    
    config_file.write(empty_line_string)

    gridFileName_string = "gridFileName = '{0}'".format(cfg['gridFileName'])
    config_file.write(gridFileName_string+'\n')      

    config_file.write(empty_line_string)

    exitSliceFileName_string = "exitSliceFileName = '{0}'".format(cfg['exitSliceFileName'])
    config_file.write(exitSliceFileName_string+'\n') 

    config_file.write(empty_line_string)

    justStats_string = "justStats = {0}".format(cfg['justStats'])
    config_file.write(justStats_string+'\n')
    
    config_file.write(empty_line_string)

    blockMarching_string = "blockMarching = {0}".format(cfg['blockMarching'])
    config_file.write(blockMarching_string+'\n')
    
    config_file.write(empty_line_string)

    nni_string = "nni = {0}".format(cfg['nni'])
    config_file.write(nni_string+'\n')
    
    nnj_string = "nnj = {0}".format(cfg['nnj'])
    config_file.write(nnj_string+'\n')

    config_file.write(empty_line_string)    

    nbi_string = "nbi = {0}".format(cfg['nbi'])
    config_file.write(nbi_string+'\n')
    
    nbj_string = "nbj = {0}".format(cfg['nbj'])
    config_file.write(nbj_string+'\n')
    
    config_file.write(empty_line_string)
    
    bx_string = "bx = {0}".format(cfg['bx'])
    config_file.write(bx_string+'\n')
    
    by_string = "by = {0}".format(cfg['by'])
    config_file.write(by_string+'\n')

    config_file.write(empty_line_string)
    
    max_time_string = "max_time = {0}".format(cfg['max_time'])
    config_file.write(max_time_string+'\n')  
    
    config_file.write(empty_line_string)
    
    max_step_string = "max_step = {0}".format(cfg['max_step'])
    config_file.write(max_step_string+'\n')
    
    config_file.write(empty_line_string)
    
    Tw_string = "Tw = {0}".format(cfg['Tw'])
    config_file.write(Tw_string+'\n')

    config_file.write(empty_line_string)
    
    BLTrans_string = "BLTrans = '{0}'".format(cfg['BLTrans'])
    config_file.write(BLTrans_string+'\n')   
    
    config_file.write(empty_line_string)
    
    TurbVisRatio_string = "TurbVisRatio = {0}".format(cfg['TurbVisRatio'])
    config_file.write(TurbVisRatio_string+'\n')  
    
    config_file.write(empty_line_string)
    
    TurbInten_string = "TurbInten = {0}".format(cfg['TurbInten'])
    config_file.write(TurbInten_string+'\n')  
    
    config_file.write(empty_line_string)
    
    coreRfraction_string = "coreRfraction = {0}".format(cfg['coreRfraction'])
    config_file.write(coreRfraction_string+'\n')  

    config_file.close()

    return                         

    


