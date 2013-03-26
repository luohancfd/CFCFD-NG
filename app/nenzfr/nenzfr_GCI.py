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

VERSION_STRING = "13-Jan-2013"

import string
import sys, os, gzip, glob
import optparse
from nenzfr_utils import read_nenzfr_outfile, run_command
from nenzfr_stats import get_slice_data
from e3_flow import StructuredGridFlow
import numpy as np
import math
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

#def calculate_grid_area(contourFile):
#    """
#    
#    """
#    fp = open(contourFile)
#    
#    return area
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
    
    workingDir = os.getcwd()
    
    os.chdir(dataFolder)
    
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
    #print data['volume']
    total_volume = sum(data['volume'])
    h = np.sqrt( total_volume/np.shape(data['volume'])[0] )
    
    os.chdir(workingDir)
    return h

def main():
    """
    Examine the command-line options, load the necessary
    data and then calculate the Grid Convergence Index for
    each nozzle exit flow property.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    
    op.add_option('--job',dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]") 
    op.add_option('--coarse', dest='coarse', 
                  help=("(path)name of folder that contains "
                        "the nenzfr simulation on a coarse grid."))
    op.add_option('--medium', dest='medium',
                  help=("(path)name of folder that contains "
                        "the nenzfr simulation on a medium grid."))
    op.add_option('--fine', dest='fine',
                  help=("(path)name of folder that contains "
                        "the nenzfr simulation on a fine grid."))
    opt, args = op.parse_args()
    
    # Read in nozzle exit flow data
    coarseGridFlowData, exitProperty = \
                  read_nenzfr_outfile(opt.coarse+"/"+opt.jobName+"-exit.stats")
    mediumGridFlowData, exitProperty = \
                  read_nenzfr_outfile(opt.medium+"/"+opt.jobName+"-exit.stats")
    fineGridFlowData, exitProperty = \
                  read_nenzfr_outfile(opt.fine+"/"+opt.jobName+"-exit.stats") 
    
    # Calculate the cell size for each simulation
    #h3 = calculate_representative_cell_size(opt.coarse)
    #h2 = calculate_representative_cell_size(opt.medium)
    #h1 = calculate_representative_cell_size(opt.fine)
    #
    #print "h1=",h1
    #print "h2=",h2
    #print "h3=",h3
    #print "h3/h1=",h3/h1
    
    h3 = calculate_representative_cell_size_v2(opt.coarse, opt.jobName)
    h2 = calculate_representative_cell_size_v2(opt.medium, opt.jobName)
    h1 = calculate_representative_cell_size_v2(opt.fine, opt.jobName)
    
    # Check that the grids follow the recommendation
    assert h3/h1 > 1.3
    print "h1=",h1
    print "h2=",h2
    print "h3=",h3
    print "h3/h1=",h3/h1  
    # Calculate ratios of grid cell sizes
    r21 = h2/h1
    r32 = h3/h1
    
    var = 'p'
    epsilon32 = coarseGridFlowData[var] - mediumGridFlowData[var]
    epsilon21 = mediumGridFlowData[var] - fineGridFlowData[var]
    
    s = math.copysign(1, epsilon32/epsilon21)
    
    print "s=",s
    print "epsilon32=",epsilon32
    print "epsilon21=",epsilon21
    
    p_old = 1.0 #1/np.log(r21) * np.abs( 2*np.log( np.abs(epsilon32/epsilon21) ) )
    
    q_new = np.log( (r21**p_old-s)/(r32**p_old-s) )
    p_new = 1/np.log(r21) * np.abs( np.log( np.abs(epsilon32/epsilon21) ) + q_new )
    
    count = 0
    #while np.abs(p_new - p_old) > 1.0e-5: #and count < 100:
    #    count += 1
    #    p_old = p_new
    #    q_new = np.log( (r21**p_old-s)/(r32**p_old-s) )
    #    p_new = 1/np.log(r21) * np.abs( np.log( np.abs(epsilon32/epsilon21) ) + q_new )
    #    #print np.abs(p_new-p_old)
    #    #print count
    
    #def error_in_p(x, e32=epsilon32, e21=epsilon21, r32=r32, r21=r21, s=s):
        
    #p = secant(error_in_p, 1.1, 1.5, tol=1.0e-4,\
    #           limits=[0,10])
    
    ## The easiest way to calculate the area is to use e3post to generate a
    ## slice across the entire nozzle flow-field. The resulting data file 
    ## contains "volume" information for each cell. For 2D grids this is
    ## actually the area (method tip from Rowan G. 14-01-2013)
    #workingDir = os.getcwd()
    # 
    #os.chdir(opt.coarse)
    ## e3post will require the gas model file so we search for it...
    #if os.path.exists("gas-model.lua"):
    #    gmodelFile = 'gas-model.lua'
    #else:
    #    gmodelFile = glob.glob("cea-lut-*.lua.gz")[0]
    #run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' %
    #                  (opt.jobName, gmodelFile))
    #            +('--output-file=%s ' % (opt.jobName+".data"))
    #            +'--slice-list="0,:,:,0" ')
    #
    ## Now read in data from newly created slice file
    #variable_list, data = get_slice_data(opt.jobName+".data")
    ##print data['volume']
    #total_volume = sum(data['volume'])
    #h = np.sqrt( total_volume/np.shape(data['volume'])[0] )
    #print total_volume
    #print np.shape(data['volume'])
    #print h
    #os.chdir(workingDir)


    #area = calculate_grid_area()


    
    var = 'p'
    epsilon32 = coarseGridFlowData[var] - mediumGridFlowData[var]
    epsilon21 = mediumGridFlowData[var] - fineGridFlowData[var]
    
    
    #s = 
    
    return 0
    
if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr GCI:\n Calculate Grid Convergence Index"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
