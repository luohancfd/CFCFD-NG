#!/usr/bin/env python
# nenzfr_GCI.py
#
# This script calculates the Grid Convergence Index (GCI)
# for each nozzle exit flow property. The calculation 
# method follows that of:
#
#
#
# Luke Doherty
# School of Mechanical and Mining Engineering
# The University of Queensland

VERSION_STRING = "23-March-2013"

import string
import sys, os, gzip, glob, copy
import optparse

from nenzfr_utils import read_nenzfr_outfile, run_command
from nenzfr_stats import get_slice_data
from nenzfr_sensitivity import add_extra_variables
from e3_flow import StructuredGridFlow

import numpy as np
import math
#from matplotlib import pyplot as plt

from scipy.optimize import curve_fit
from scipy.stats import linregress

E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#---------------------------------------------------------------------
def calculate_representative_cell_size(dataFolder):
    """
    """
    fileList = glob.glob(dataFolder+"/flow/t9999/*.gz")
    GridFlowData = []
    blk_volume = []
    blk_cells = []
    for f in fileList:
        f1 = gzip.open(f,"r")
        GridFlowData.append(StructuredGridFlow())
        GridFlowData[-1].read(f1)
        blk_volume.append( sum(sum(GridFlowData[-1].data['volume'][:,:,:])) )
        blk_cells.append( GridFlowData[-1].ni * GridFlowData[-1].nj )
        f1.close()

    total_volume = sum(blk_volume)
    total_cells = sum(blk_cells)
    if dataFolder in ['case00_fine']:
        print "total_cells=",total_cells
        print "total_volume=",total_volume
        
    h = np.sqrt( total_volume/float(total_cells) )
    
    return h

def calculate_representative_cell_size_v2(dataFolder, jobName):
    """
    The easiest way to calculate the area is to use e3post to generate a
    slice across the entire nozzle flow-field. The resulting data file
    contains "volume" information for each cell. For 2D grids this is
    actually the area (method tip from Rowan G. 14-01-2013)
    """
    # Save the current working directory for later use 
    workingDir = os.getcwd()
    # Change into the simulation directory
    os.chdir(dataFolder)
    
    # We check first that the required slice file has not already been
    # created.
    if not os.path.exists(jobName+".data"):
        # e3post will require the gas model file so we search for it...
        if os.path.exists("gas-model.lua"):
            gmodelFile = 'gas-model.lua'
        else:
            gmodelFile = glob.glob("cea-lut-*.lua.gz")[0]
        run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' %
                          (jobName, gmodelFile))
                    +('--output-file=%s ' % (jobName+".data"))
                    +'--slice-list=":,:,:,0" ')
    
    # Now read in data from the slice file
    variable_list, data = get_slice_data(jobName+".data")
    # Calculate the total area (volume) and determine the 
    # representative cell size
    total_volume = sum(data['volume'])
    h = np.sqrt( total_volume/np.shape(data['volume'])[0] )
    # Change back into the working directory
    os.chdir(workingDir)
    return h

def read_caseList_file(FileToRead):
    """
    Short function that reads a file and returns
    a list ignoring lines that begin with a # or
    that are blank.
    """
    fp = open(FileToRead)
    data = fp.readlines()
    fp.close()
    
    case_list = []
    for line in data:
        if not line.startswith('#') and line.strip() not in [""]:
            case_list.append(line.strip())
    return case_list

#---------------------------------------------------------------------------
# The following functions were defined and used when I was trying to fit an
# aymptotic curve to the data. For the currently implemented models they are 
# not required.
def func_decreasing(x, a, b, c):
    """
    Function to fit when y decreases
    with x i.e. approaches asymptote 
    from above.
    """
    return a + b*np.exp(c*x)

def func_increasing(x, a, b, c):
    """
    Function to fit when y increases 
    with x i.e. approaches asymptote
    from below.
    """
    #return a + b*np.exp(-c/x)
    return a - b*np.exp(c*x)

def func(x, a, b):
    #return a - b*np.exp(c*x)
    #return a + b*np.power(x,c)
    return a + b*x

#---------------------------------------------------------------------------
def main():
    """
    Examine the command-line options, load the necessary data and then calculate 
    the Grid Convergence Index for each nozzle exit flow property.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    
    op.add_option('--run-defaults', dest='runDefaults', action='store_true',
                  default=True, help="calculate GCI using defaults for all other inputs.")

    op.add_option('--job',dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]") 
    
    op.add_option('--caseListFile',dest='caseListFile', default='case_list.txt',
                  help=("the name of a text file containing a list of relative pathnames "
                       "for the simulations whose grid error is to be calculated. Files "
                       "should be listed from finest grid to coarsest grid. Comment "
                       "lines should begin with a #. [default: %default]"))
    
    op.add_option('--casesToFitFile',dest='caseFitList', default='case_list.txt',
                   help=("the name of a text file containing a list of relative pathnames "
                         "for simulations that will be used to generate a line-of-best"
                         "-fit. Only applicable when --method='linear-fit'. Also, the "
                         "list MUST be a sub-set of the files given in --caseListFile."
                         "[default: %default]"))
 
    op.add_option('--safety-factor',dest='safetyFactor', type='float', default=1.0,
                  help=("specify a safety-factor to be applied to the calculated "
                        "convergence error [default: %default]"))
    
    op.add_option('--method', dest='method', choices=['linear-fit','ref-solution'],
                  default='linear-fit', help=("specify the calculation method. "
                  "Choices are (1) 'linear-fit' : fit a straight line and extrapolate "
                  "cell size = zero, or (2) 'ref-solution' : consider the solution on "
                  "the finest grid as the 'truth'. [default: %default]" ))
    
    opt, args = op.parse_args() 

    # Read the case list file
    case_list = read_caseList_file(opt.caseListFile)
    
    # If required read in the list of cases that will be used to 
    # generate a line-of-best fit.
    if opt.method in ['linear-fit',]:
        if opt.caseFitList not in ['case_list.txt',]:
            case_fit_list = read_caseList_file(opt.caseFitList)
            print "Linear model will be fitted using simulations given in ",opt.caseFitList
        else:
            case_fit_list = copy.copy(case_list) 
            print "Linear model will be fitted using simulations given in ",opt.caseListFile
        
        if len(case_fit_list) == 1:
            print "ERROR: case_fit_list contains only one entry."
            print "       This is not enough to generate a linear fit."
            print "       Check input --caseToFitFile"
            return -2
        if len(case_fit_list) == 2:
            print "WARNING: case_fit_list contains only two entries."
            print "         Fitted line may not be robust."
        
    else:
        print "Using ",case_list[0]," as the reference ('truth') solution"
    
    # Read in nozzle exit flow data
    print "Reading exit flow files..."
    FlowData = {}
    exitProperty = {}
    for case in case_list:
        data, props = read_nenzfr_outfile(case+opt.jobName+"-exit.stats")
        # Add dynamic pressure, unit Reynolds number, mass flux and 
        # static-to-dynamic pressure ratio
        data, props = add_extra_variables(data, props)
        FlowData[case] = data
        exitProperty[case] = props 
     
    # Calculate the cell size for each simulation and write a summary
    # file. If a summary file already exists we read that rather then
    # re-calculate the cell size.
    if not os.path.exists("cell_sizes.dat"):
        print "Calculating representative cell size..."
        # Initialise cellsize dictionary and x array...
        cellsize = {}
        x = np.array([])
        # Begin writing a summary file
        fp = open("cell_sizes.dat","w")
        fp.write("# Representative Cell Sizes\n")
        fp.write("# Columns are:\n")
        fp.write("#    case, cell-size (metre)\n")
        # Loop through each case, calculate cell size and its 
        # reciprocal (x = 1/h) and write data to the summary file
        for case in case_list:
            h = calculate_representative_cell_size_v2(case, opt.jobName)
            cellsize[case] = h
            #x = np.append(x, h)
            #fp.write("{0:s}\t{1:1.10e}\t{2:1.10e}\n".format(case,h,1.0/h))
            fp.write("{0:s}\t{1:1.10e}\n".format(case,h))
            print case, ": ",h
        fp.close()
    else:
        print "Reading cell size data from 'cell_sizes.dat'..."
        fp = open("cell_sizes.dat","r")
        contents = fp.readlines()
        fp.close()
        cellsize = {}
        x = np.array([])
        for line in contents:
            if not line.startswith("#"):
                lineData = line.strip().split("\t")
                cellsize[lineData[0]] = np.float(lineData[1])
                #x = np.append(x, np.float(lineData[2]))
                print lineData[0], ": ",lineData[1]
    

    # Now we set about calculating the grid convergence error for
    # each exit-flow property using the desired method
    print "Calculating grid errors for each flow property..."
    linear_fit_params = {}
    grid_errors = {}
    true_value = {}
    for exitVar in exitProperty[case_list[0]]:
        grid_errors[exitVar] = {}
        if opt.method in ['linear-fit',]:
            # Assemble "x" and "y" vectors
            y = np.array([])
            x = np.array([])
            for case in case_fit_list:
                y = np.append(y, FlowData[case][exitVar])
                x = np.append(x, cellsize[case])

            # Straight line fit through data
            slope, intercept, r_value, p_value, std_err = linregress(x,y)
            
            # Store the parameters of the fitterd curve so we may 
            # write them to a file later
            linear_fit_params[exitVar] = {'slope':slope,'intercept':intercept,\
                                      'r_value':r_value,'p_value':p_value,\
                                      'std_err':std_err}
            # Stort the "true" value for later use
            true_value[exitVar] = intercept
            
            # Error relative to the linear fit intercept
            for case in case_list:
                if intercept != 0.:
                    grid_errors[exitVar][case] = \
                               (FlowData[case][exitVar] - intercept)/intercept *\
                               100.0 * opt.safetyFactor
                else:
                    grid_errors[exitVar][case] = float('NaN')
            
        elif opt.method in ['ref-solution',]:
            finest = case_list[0]
            true_value[exitVar] = FlowData[finest][exitVar]
            # Error relative to the finest grid solution
            for case in case_list:
                if true_value[exitVar] != 0.:
                    grid_errors[exitVar][case] = \
                               (FlowData[case][exitVar] - true_value[exitVar])/\
                                true_value[exitVar] * 100.0 * opt.safetyFactor
                else:
                    grid_errors[exitVar][case] = float('NaN')

        elif opt.method in ['roacheGCI']:
            print "Calculation method roacheGCI has not yet been implemented"
            continue
            ## Check that the grids follow the recommendation
            #assert h3/h1 > 1.3
            #print "h1=",h1
            #print "h2=",h2
            #print "h3=",h3
            #print "h3/h1=",h3/h1
            ## Calculate ratios of grid cell sizes
            #r21 = h2/h1
            #r32 = h3/h1
            #
            #exitVar = 'p'
            #epsilon32 = coarseGridFlowData[exitVar] - mediumGridFlowData[exitVar]
            #epsilon21 = mediumGridFlowData[exitVar] - fineGridFlowData[exitVar]
            #
            #s = math.copysign(1, epsilon32/epsilon21)
            #
            #print "s=",s
            #print "epsilon32=",epsilon32
            #print "epsilon21=",epsilon21
            #
            #p_old = 1.0 #1/np.log(r21) * np.abs( 2*np.log( np.abs(epsilon32/epsilon21) ) )
            #
            #q_new = np.log( (r21**p_old-s)/(r32**p_old-s) )
            #p_new = 1/np.log(r21) * np.abs( np.log( np.abs(epsilon32/epsilon21) ) + q_new )
            #
            #count = 0
            ##while np.abs(p_new - p_old) > 1.0e-5: #and count < 100:
            ##    count += 1
            ##    p_old = p_new
            ##    q_new = np.log( (r21**p_old-s)/(r32**p_old-s) )
            ##    p_new = 1/np.log(r21) * np.abs( np.log( np.abs(epsilon32/epsilon21) ) + q_new )
            ##    #print np.abs(p_new-p_old)
            ##    #print count
            #
            ##def error_in_p(x, e32=epsilon32, e21=epsilon21, r32=r32, r21=r21, s=s):
            #
            ##p = secant(error_in_p, 1.1, 1.5, tol=1.0e-4,\
            ##           limits=[0,10])
            #
            #exitVar = 'p'
            #epsilon32 = coarseGridFlowData[exitVar] - mediumGridFlowData[exitVar]
            #epsilon21 = mediumGridFlowData[exitVar] - fineGridFlowData[exitVar]

    # Now write out some summary files
    print "Writing data files..."
    if opt.method in ['linear-fit']:
        # For the linear fit method we write a file containing
        # all the regression results
        fp = open('linear-fit-regression-params.dat','w')
        fp.write('# property:        slope:    intercept:   r_value:    p_value:      std_err\n')
        for exitVar in exitProperty[case_list[0]]:
            fp.write('{0:>12s}  '.format(exitVar))
            fp.write('{0:>12g}  '.format(linear_fit_params[exitVar]['slope']))
            fp.write('{0:>12g}  '.format(linear_fit_params[exitVar]['intercept']))
            fp.write('{0:>9g}  '.format(linear_fit_params[exitVar]['r_value']))
            fp.write('{0:>10g}  '.format(linear_fit_params[exitVar]['p_value']))
            fp.write('{0:>12g}\n'.format(linear_fit_params[exitVar]['std_err']))
        fp.close()
    
    # Summary of the errors
    fp = open("grid-convergence-errors-"+opt.method+".dat",'w')
    fp.write('# safety-factor: {0:1.5f}\n'.format(opt.safetyFactor))
    fp.write('#: property: "true-value": ')
    for case in case_list:
        if case not in [case_list[-1],]:
            fp.write('{0:s}: '.format(case))
        else:
            fp.write('{0:s}\n'.format(case))
    for exitVar in exitProperty[case_list[0]]:
        fp.write('{0:>12s}  '.format(exitVar))
        fp.write('{0:>12g}  '.format(true_value[exitVar]))
        for case in case_list:
            if case not in [case_list[-1],]:
                fp.write('{0:>9.4f}  '.format(grid_errors[exitVar][case]))
            else:
                fp.write('{0:>9.4f}\n'.format(grid_errors[exitVar][case]))
    fp.close()
    
     
    jnk

    

    x = np.array([1./300.,1./200.,1./100.,1./100,1./80.])
    
    xFit = np.linspace(0,np.max(x),1000)
    #for exitVar in exitProperty[case_list[0]]:
    exitVar = 'p'
    # Assemble "y" vectors
    y = np.array([])
    for case in case_list:
        y = np.append(y, FlowData[case][exitVar])
    print
    print x
    print
    print x/x[0]
    print
    print y
    
    # Calculate an initial guess for the equation coefficients
    ##a = y[0]/2
    #b = 1.0 #-( 0.5*(x[-1]+x[0])/(x[-1]-x[0])*y[0] - x[0]/(x[-1]-x[0])*y[-1] )
    #c = (y[-1] - y[0])/(x[0] - x[-1])
    #a = y[-1] + 1 + c*x[-1]
    ##b = 1.0
    ##c = (y[0] - y[-1])/(x[0] - x[-1])
    ##a = y[0] - 1 - c*x[0]
    
    #print "a=",a," b=",b," c=",c
    
    #yTestDec = func_decreasing(xFit,a,b,c)
    #yTestInc = func_increasing(xFit,a,b,c)
    #yTest = func(xFit,a,5*b,c/10.)
    
    #plt.plot(x,y,'.b',xFit,yTest,'-r') #yTestDec,'-r',xFit,yTestInc,'-k')
    ##plt.xscale('log')
    ##plt.yscale('log')
    #plt.ylim(600,700)
    #plt.show()
    
    #jnk
    
    #if y[0] > y[-1]:
    #    # Infer that the data is approaching the asymptote
    #    # from above.
    #    print "Decreasing function fitted for ",exitVar
    #    try: 
    #        popt, pcov = curve_fit(func_decreasing, x, y, p0=[y[-1],1.,5000.])
    #
    #        yFit = func_decreasing(xFit, popt[0], popt[1], popt[2])
    #    except RuntimeError:
    #        print "Error - curve_fit failed for ",exitVar
    #else:
    #    # Infer that the data is approaching the asymptote
    #    # from below.
    #    print "Increasing function fitted for ",exitVar
    #    #try:
    #    #    popt, pcov = curve_fit(func, x, y, p0=[a,5.*b,c/10.])
    #    #    yFit = func(xFit, popt[0], popt[1], popt[2])
    #    #    #popt, pcov = curve_fit(func_increasing, x, y, p0=None)
    #    # 
    #    #    #yFit = func_decreasing(xFit, popt[0], popt[1], popt[2])
    #    #except RuntimeError:
    #    #    print "Error - curve_fit failed for ",exitVar
    popt, pcov, infodict, errmsg, ier = \
          curve_fit(func, x, y, p0=[y[0],1.0], maxfev=5000, full_output=1)
    yFit = func(xFit, popt[0], popt[1])
    
    #print "a_opt=",popt[0]," b_opt=",popt[1]," c_opt=",popt[2]
    #print "a_var=",pcov[0]," b_var=",pcov[1]," c_var=",pcov[2]
    print "opt=",popt
    print "cov=",pcov
    print exitVar,"_true estimate=",func(0,popt[0],popt[1])
    print infodict
    
    plt.plot(x,y,'.b',xFit,yFit,'-r')
    plt.ylim(600,700)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.show() 
    
    
    
    
    return 0
    
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr GCI:\n Calculate Grid Convergence Index"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
