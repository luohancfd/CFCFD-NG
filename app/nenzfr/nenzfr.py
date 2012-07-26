#!/usr/bin/env python
"""
nenzfr.py -- NENZF reloaded (Non-Equilibrium Nozzle Flow reloaded)

This script coordinates the running of the T4 nozzle calculation.
The intention is to provide a fairly quick calculation of 
the test flow conditions at the exit plane of the selected nozzle.

Behind the scene, estcj.py is used to get an estimate of 
the flow condition at the nozzle throat and then Eilmer3 is used
to compute the expanding flow in the divergent part of the nozzle.
Finally, a profile is examined at the downstream-end of the nozzle
to extract nominal flow condition data.

.. Authors: Peter Jacobs, Luke Doherty, Wilson Chan and Rainer Kirchhartz
   School of Mechanical and Mining Engineering
   The University of Queensland

.. TODO: Luke, we need to think of a good way to encode all of the options.
   I've added quite a few command-line options and more are needed.
   Maybe a dictionary of customised parameters for each case. Wilson.
"""

VERSION_STRING = "08-May-2012"

from string import upper
import sys, os
import optparse
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)

# We've put a most of the functions in other places so that this file
# becomes just the top-level coordinating function that wrangles
# the command-line options and works out what to do.
from nenzfr_utils import prepare_input_script, run_command, quote, \
                         read_gmodelFile_from_config
from nenzfr_stats import *
from nenzfr_parallel import run_in_block_marching_mode, read_block_dims

#---------------------------------------------------------------

def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by estcj and Eilmer3.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"
                        " [default: %default]"))
    op.add_option('--p1', dest='p1', type='float', default=None,
                  help=("shock tube fill pressure or static pressure, in Pa"))
    op.add_option('--T1', dest='T1', type='float', default=None,
                  help=("shock tube fill temperature, in degrees K"))
    op.add_option('--Vs', dest='Vs', type='float', default=None,
                  help=("incident shock speed, in m/s"))
    op.add_option('--pe', dest='pe', type='float', default=None,
                  help=("equilibrium pressure (after shock reflection), in Pa"))
    op.add_option('--chem', dest='chemModel', default='eq',
                  choices=['eq', 'neq', 'frz', 'frz2'], 
                  help=("chemistry model: " "eq=equilibrium; " 
                        "neq=non-equilibrium; " "frz=frozen " 
                        "[default: %default]"))
    op.add_option('--area', dest='areaRatio', default=1581.165,
                  help=("nozzle area ratio. only used for estcj calc. "
                        "use when --cfile(--gfile) are "
                        "specified. [default: %default]"))
    op.add_option('--job', dest='jobName', default='nozzle',
                  help="base name for Eilmer3 files [default: %default]")
    op.add_option('--cfile', dest='contourFileName', 
                  default='Bezier-control-pts-t4-m10.data',
                  help="file containing Bezier control points "
                       "for nozzle contour [default: %default]")
    op.add_option('--gfile', dest='gridFileName', default='None',
                  help="file containing nozzle grid. "
                  "overrides --cfile if both are given "
                  "[default: %default]")
    op.add_option('--exitfile', dest='exitSliceFileName', 
                  default='nozzle-exit.data',
                  help="file for holding the nozzle-exit data [default: %default]")
    op.add_option('--just-stats', dest='justStats', action='store_true', 
                  default=False,
                  help="skip the detailed calculations and "
                  "just retrieve exit-flow statistics")
    op.add_option('--block-marching', dest='blockMarching', action='store_true', 
                  default=False, help="run nenzfr in block-marching mode")
    # The following defaults suit Luke's Mach 10 calculations.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells [default: %default]"))
    op.add_option('--nnj', dest='nnj', type='int', default=100,
                  help=("number of radial cells [default: %default]"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk) "
                        "[default: %default]"))
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help=("number of radial blocks [default: %default]"))
    op.add_option('--bx', dest='bx', type='float', default=1.10,
                  help=("clustering in the axial direction [default: %default]"))
    op.add_option('--by', dest='by', type='float', default=1.002,
                  help=("clustering in the radial direction [default: %default]"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow [default: %default]"))
    op.add_option('--max-step', dest='max_step', type='int', default=800000,
                  help=("maximum simulation steps allowed [default: %default]"))
    #
    op.add_option('--Twall', dest='Tw', type='float', default=300.0,
                  help=("Nozzle wall temperature, in K "
                        "[default: %default]"))
    op.add_option('--BLTrans', dest='BLTrans', default="x_c[-1]*1.1",
                  help=("Transition location for the Boundary layer. Used "
                        "to define the turbulent portion of the nozzle. "
                        "[default: >nozzle length i.e. laminar nozzle]"))
    op.add_option('--TurbVisRatio', dest='TurbVisRatio', type='float',
                  default=100.0, help=("Turbulent to Laminar Viscosity Ratio "
                  "[default: %default]"))
    op.add_option('--TurbIntensity', dest='TurbInten', type='float', 
                  default=0.05, help=("Turbulence intensity at the throat "
                  "[default: %default]"))
    op.add_option('--CoreRadiusFraction', dest='coreRfraction', type='float',
                  default=2.0/3.0, help=("Radius of core flow as a fraction of "
                  "the nozzle exit radius [default: %default]"))
    opt, args = op.parse_args()
    #
    # Get the nozzle contour file into the current work directory.
    if not os.path.exists(opt.contourFileName):
        run_command('cp '+E3BIN+'/nenzfr_data_files/'+opt.contourFileName+' .')
    # Set up the equilibrium gas-model file as a look-up table.
    if opt.chemModel in ['eq',]:
        if opt.gasName in ['n2']:
            eqGasModelFile = 'cea-lut-'+upper(opt.gasName)+'.lua.gz'
        else:
            eqGasModelFile = 'cea-lut-'+opt.gasName+'.lua.gz'
        if not os.path.exists(eqGasModelFile):
            run_command('build-cea-lut.py --gas='+opt.gasName)
        gmodelFile = eqGasModelFile
    else:
        # We'll assume that the gas-model file of default name is set up.
        # TODO: Luke, this needs to be modified, I suspect.
        gmodelFile = 'gas-model.lua'
    #
    # If we have already run a calculation, it may be that we just want
    # to extract the exit-flow statistics again.
    if opt.justStats:
        gmodelFile = read_gmodelFile_from_config(opt.jobName)
        print_stats(opt.exitSliceFileName,opt.jobName,opt.coreRfraction,gmodelFile)
        return 0
    #
    # Go ahead with a new calculation.
    # First, make sure that we have the needed parameters.
    bad_input = False
    if opt.p1 is None:
        print "Need to supply a float value for p1."
        bad_input = True
    if opt.T1 is None:
        print "Need to supply a float value for T1."
        bad_input = True
    if opt.Vs is None:
        print "Need to supply a float value for Vs."
        bad_input = True
    if opt.pe is None:
        print "Need to supply a float value for pe."
        bad_input = True
    if bad_input:
        return -2
    #
    # Runs estcj to get the equilibrium shock-tube conditions up to the nozzle-supply region.
    command_text = E3BIN+('/estcj.py --gas=%s --T1=%g --p1=%g --Vs=%g --pe=%g --task=st --ofn=%s' % 
                          (opt.gasName, opt.T1, opt.p1, opt.Vs, opt.pe, opt.jobName))
    run_command(command_text)
    # Switch off block-sequencing flag if we have selected to run in MPI block-marching mode.
    if opt.blockMarching: 
        sequenceBlocksFlag = 0
    else: 
        sequenceBlocksFlag = 1
        opt.nbj = 1
    # Set up the input script for Eilmer3.
    paramDict = {'jobName':quote(opt.jobName), 'gasName':quote(opt.gasName),
                 'T1':opt.T1, 'p1':opt.p1, 'Vs':opt.Vs, 'pe':opt.pe,
                 'contourFileName':quote(opt.contourFileName),
                 'gridFileName':quote(opt.gridFileName), 
                 'chemModel':quote(opt.chemModel),
                 'areaRatio':opt.areaRatio, 'seq_blocks':sequenceBlocksFlag,
                 'nni':opt.nni, 'nnj':opt.nnj, 'nbi':opt.nbi, 'nbj':opt.nbj,
                 'bx':opt.bx, 'by':opt.by,
                 'max_time':opt.max_time, 'max_step':opt.max_step, 
                 'HOME':'$HOME', 'Tw':opt.Tw, 'TurbVisRatio':opt.TurbVisRatio,
                 'TurbInten':opt.TurbInten, 'BLTrans':opt.BLTrans}
    prepare_input_script(paramDict, opt.jobName)
    # Run the simulation code.
    if opt.blockMarching:
        # Run Eilmer3 either in the multi-processor block-marching mode
        # where the sets of blocks are set up by the functions in
        # nenzfr_parallel.py and a C++ simulation is done for each set.
        run_in_block_marching_mode(opt, gmodelFile)
        # Generate slice list for exit plane.
        exitPlaneSlice = '-' + str(opt.nbj) + ':-1,-2,:,0'
        # Generate slice list for centreline.
        blockDims, noOfBlks = read_block_dims(opt.jobName + '.config')
        centrelineSlice = '0,:,1,0'
        for blk in range(opt.nbj, noOfBlks, opt.nbj):
            centrelineSlice += ';' + str(blk) + ',:,1,0'
    else:
        # Run Eilmer3 either in the single processor mode
        # where the block-sequencing happens in the C++ code.
        run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
        run_command(E3BIN+('/e3shared.exe --job=%s --run' % (opt.jobName,)))
        # Generate slice list for exit plane and centreline.
        exitPlaneSlice = '-1,-2,:,0'
        centrelineSlice = ':,:,1,0'
    #
    # Exit plane slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s ' % (opt.exitSliceFileName,))
                +('--slice-list="%s" ' % exitPlaneSlice)
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Centerline slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s-centreline.data ' % (opt.jobName,))
                +('--slice-list="%s" ' % centrelineSlice)
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Prep files for plotting with Paraview            
    run_command(E3BIN+('/e3post.py --job=%s --vtk-xml --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Generate averaged exit flow properties
    print_stats(opt.exitSliceFileName,opt.jobName,opt.coreRfraction,gmodelFile)                
    # Compute viscous data at the nozzle wall
    run_command(E3BIN+'/nenzfr_compute_viscous_data.py --job=%s --nbj=%s' % (opt.jobName, opt.nbj))
    #
    if opt.contourFileName == "contour-t4-m10.data":
        # The following are additional commands specific to Luke D. and the Mach 10 nozzle.
        #
        # Extract a slice from the last block along jk index directions at the i-index that
        # is closest to the point x=1.642, y=0.0, z=0.0
        run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                           (opt.jobName, gmodelFile))
                   +('--output-file=%s2 --slice-at-point="-1,jk,1.642,0,0" ' % 
                     (opt.exitSliceFileName,))
                   +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
        #
        # 29-07-2011 Luke D. 
        # Copy the original exit stats file to a temporary file
        run_command(('cp %s-exit.stats %s-exit.stats_temp') % (opt.jobName, opt.jobName,))
        # (Re) Generate the exit stats using the slice-at-point data
        print_stats(opt.exitSliceFileName+'2',opt.jobName,opt.coreRfraction,gmodelFile)
        run_command('cp %s-exit.stats %s-exit.stats2' % (opt.jobName, opt.jobName,))
        # Now rename the temporary exit stats file back to its original name
        run_command('mv %s-exit.stats_temp %s-exit.stats' % (opt.jobName, opt.jobName,))
        # End specific commands
    #
    return 0
    
#---------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr: Calculate Shock Tunnel Test Flow Conditions"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
