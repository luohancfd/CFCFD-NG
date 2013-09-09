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

.. Authors: Peter Jacobs, Luke Doherty, Wilson Chan, Rainer Kirchhartz,
   and Chris James.
   School of Mechanical and Mining Engineering
   The University of Queensland

.. TODO: Luke, we need to think of a good way to encode all of the options.
   I've added quite a few command-line options and more are needed.
   Maybe a dictionary of customised parameters for each case. Wilson.
"""

VERSION_STRING = "13-May-2013"

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
from nenzfr_input_utils import input_checker

#---------------------------------------------------------------

def main(cfg={}):
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by estcj and Eilmer3.
    """

    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('-c', '--config_file', dest='config_file',
                  help=("filename for the config file"))
    opt, args = op.parse_args()
    config_file = opt.config_file
       
    if not cfg: #if the configuration dictionary has not been filled up already, load it from a file
        try: #from Rowan's onedval program
            execfile(config_file, globals(), cfg)
        except IOError:
            print "There was a problem reading the config file: '{0}'".format(config_file)
            print "Check that it conforms to Python syntax."
            print "Bailing out!"
            sys.exit(1)
               
    # First, check our input and assign any default values required:
    cfg = input_checker(cfg)
    #bail out here if there is an issue
    if cfg['bad_input']:
        return -2
    #
    # Get the nozzle contour file into the current work directory.
    if not os.path.exists(cfg['contourFileName']):
        run_command('cp '+E3BIN+'/nenzfr_data_files/'+cfg['contourFileName']+' .')
    # Set up the equilibrium gas-model file as a look-up table.
    if not cfg['justStats']:
        if cfg['chemModel'] in ['eq',]:
            if cfg['gasName'] in ['n2']:
                eqGasModelFile = 'cea-lut-'+upper(cfg['gasName'])+'.lua.gz'
            else:
                eqGasModelFile = 'cea-lut-'+cfg['gasName']+'.lua.gz'
            if not os.path.exists(eqGasModelFile):
                run_command('build-cea-lut.py --gas='+cfg['gasName'])
            gmodelFile = eqGasModelFile
        else:
            # We'll assume that the gas-model file of default name is set up.
            # TODO: Luke, this needs to be modified, I suspect.
            gmodelFile = 'gas-model.lua'
    #
    # If we have already run a calculation, it may be that we just want
    # to extract the exit-flow statistics again.
    if cfg['justStats']:
        gmodelFile = read_gmodelFile_from_config(cfg['jobName'])
        print_stats(cfg['exitSliceFileName'],cfg['jobName'],cfg['coreRfraction'],gmodelFile)
        return 0
    #
    # Here each facility type will be doing different things (no estcj for expansion tube)
    if cfg['facility'] == 'reflected-shock-tunnel':
        #
        # Runs estcj to get the equilibrium shock-tube conditions up to the nozzle-supply region.
        command_text = E3BIN+('/estcj.py --gas=%s --T1=%g --p1=%g --Vs=%g --pe=%g --task=st --ofn=%s' % 
                              (cfg['gasName'], cfg['T1'], cfg['p1'], cfg['Vs'], cfg['pe'], cfg['jobName']))
        run_command(command_text)
    # Switch off block-sequencing flag if we have selected to run in MPI block-marching mode.
    if cfg['blockMarching']: 
        sequenceBlocksFlag = 0
    else: 
        sequenceBlocksFlag = 1
        cfg['nbj'] = 1
    # Set up the input script for Eilmer3.
    # the majority of the inputs will already be in our cfg dictionary,
    # but we'll add any things that are still needed
    # add HOME and sequenceBlocksFlag
    cfg['HOME'] = '$HOME'; cfg['seq_blocks'] = sequenceBlocksFlag
    # now need to put in empty variables for the other facilities so the
    # template file works properly
    if cfg['facility'] == 'reflected-shock-tunnel':
        # need to set the expansion tube and gun tunnel parameters to nothing for the substitution
        cfg['T7'] = None; cfg['p7'] = None; cfg['V7'] = None
        cfg['T0'] = None; cfg['p0'] = None
        cfg['pitot_input_file'] = None
    elif cfg['facility'] == 'expansion-tube':
        # need to set the reflected shock tunnel and gun tunnel parameters to nothing for the substitution
        cfg['T1'] = None; cfg['p1'] = None; cfg['Vs'] = None 
        cfg['pe'] = None; cfg['T0'] = None; cfg['p0'] = None
        if 'pitot_input_file' in cfg:
            cfg['T7'] = None; cfg['p7'] = None; cfg['V7'] = None  
        if 'pitot_input_file' not in cfg:
             cfg['pitot_input_file'] = None
    elif cfg['facility'] == 'gun-tunnel':
        # need to set the reflected shock tunnel and expansion tube parameters to nothing for the substition
        cfg['T1'] = None; cfg['p1'] = None; cfg['Vs'] = None; cfg['pe'] = None
        cfg['T7'] = None; cfg['p7'] = None; cfg['V7'] = None
        cfg['pitot_input_file'] = None
        
    # need to put an extra set of quotes around some of the strings for when
    # they are put into the eilmer 3 template
    cfg['facility'] = quote(cfg['facility'])
    cfg['contourFileName'] = quote(cfg['contourFileName'])
    cfg['gridFileName'] = quote(cfg['gridFileName'])
    cfg['chemModel'] = quote(cfg['chemModel'])  
    cfg['gasName'] = quote(cfg['gasName'])
    if cfg['pitot_input_file']:
        cfg['pitot_input_file'] = quote(cfg['pitot_input_file'])
          
    prepare_input_script(cfg, cfg['jobName'])
    run_command(E3BIN+('/e3prep.py --job=%s --do-svg --clean-start' % (cfg['jobName'],)))
    # Run the simulation code.
    if cfg['blockMarching']:
        # Run Eilmer3 either in the multi-processor block-marching mode.
        run_command(E3BIN+('/e3march.py --job=%s --run --nbj=%d' % (cfg['jobName'],cfg['nbj'])))
        # Generate slice list for exit plane.
        exitPlaneSlice = '-' + str(cfg['nbj']) + ':-1,-2,:,0'
        # Generate slice list for centreline.
        noOfBlks = None
        fp = file(cfg['jobName'] + '.config')
        for line in fp:
            if 'nblock =' in line:
                tokens = line.split()
                noOfBlks = int(tokens[2])
        if noOfBlks is None:
            raise RuntimeError("Didn't manage to find the number of blocks in config file.")
        centrelineSlice = '0,:,1,0'
        for blk in range(cfg['nbj'], noOfBlks, cfg['nbj']):
            centrelineSlice += ';' + str(blk) + ',:,1,0'
    else:
        # Run Eilmer3 either in the single processor mode
        # where the block-sequencing happens in the C++ code.
        run_command(E3BIN+('/e3shared.exe --job=%s --run' % (cfg['jobName'],)))
        # Generate slice list for exit plane and centreline.
        exitPlaneSlice = '-1,-2,:,0'
        centrelineSlice = ':,:,1,0'
    #
    # Exit plane slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=0001 --gmodel-file=%s ' % 
                       (cfg['jobName'], gmodelFile))
                +('--output-file=%s ' % (cfg['exitSliceFileName'],))
                +('--slice-list="%s" ' % exitPlaneSlice)
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Centerline slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=0001 --gmodel-file=%s ' % 
                       (cfg['jobName'], gmodelFile))
                +('--output-file=%s-centreline.data ' % (cfg['jobName'],))
                +('--slice-list="%s" ' % centrelineSlice)
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Generate files for plotting with Paraview.          
    run_command(E3BIN+('/e3post.py --job=%s --vtk-xml --tindx=0001 --gmodel-file=%s ' % 
                       (cfg['jobName'], gmodelFile))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Compute viscous data at the nozzle wall
    run_command(E3BIN+'/nenzfr_compute_viscous_data.py --job=%s --nbj=%s' % (cfg['jobName'], cfg['nbj']))
    # Generate averaged exit flow properties
    print_stats(cfg['exitSliceFileName'],cfg['jobName'],cfg['coreRfraction'],gmodelFile)                
    #
    return 0
    
#---------------------------------------------------------------

#def main(cfg={}):
#    print "in main()."
#    print "leaving main()."

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print "NENZFr: Calculate Shock Tunnel Test Flow Conditions"
        print "   Version:", VERSION_STRING
        print "   To get some useful hints, invoke the program with option --help."
        sys.exit(0)
    return_flag = main()
    sys.exit(return_flag)
