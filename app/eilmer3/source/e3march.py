#! /usr/bin/env python
"""
Block-marching version of the Eilmer3 gas-flow simulation program.

This program started as a set of shell scripts written by Wilson 
to accelerate the calculation of his combustor flows by working 
on a few blocks at a time, marching along the axial-direction of 
the combustor.  This approach should work for many flows that are
essentially contained within a single duct or have a simple 
domain and a dominant supersonic flow in essentially one direction.

Usage
-----
Command line::

  e3march.py [options]

Options::

| e3march.py [--help] [--job=<jobFileName>]
|            [--zip-files|--no-zip-files]
|            [--run] [--restart=<runNumber>]
|            [--nbj=<nbj>] [--nbk=<nbk>]
|            [--max-dt=<dt>]
|            [--polish-time=<time,dt>]
|            [--verbosity=<int>]
|            [--max-wall-clock=<seconds>]

.. Authors: Peter Jacobs, Wilson Chan, Luke Doherty, Rainer Kirchhartz,
            Chris James and Rowan Gollan.
   School of Mechanical and Mining Engineering
   The University of Queensland

.. Versions:
   04-Sep-2013: ported the essential bits over from nenzfr.py and e3prep.py
                The future nenzfr.py really shouldn't have to dig into .config 
                files and do all the coordination of running the MPI tasks.
                We'll do all of that stuff here and leave it to the higher-level
                tasks of nozzle-flow simulation.
"""

# ----------------------------------------------------------------------
#
import sys, os, gzip
import shlex, subprocess, string
from getopt import getopt
from glob import glob
import traceback
from ConfigParser import SafeConfigParser
if sys.version_info < (2, 7):
    from OrderedDict import *
import time
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)
from libprep3 import *

verbosity_level = 0

shortOptions = ""
longOptions = ["help", "job=",
               "zip-files", "no-zip-files",
               "run", "restart=",
               "nbj=", "nbk=", 
               "max-dt=",
               "polish-time=",
               "verbosity=",
               "max-wall-clock="]

def printUsage():
    print ""
    print "Usage: e3march.py [--help] [--job=<jobFileName>]"
    print "       [--zip-files|--no-zip-files]"
    print "       [--run] [--restart=<runNumber>]"
    print "       [--nbj=<nbj>] [--nbk=<nbk>]"
    print "       [--max-dt=<dt>]"
    print "       [--polish-time=<time,dt]"
    print "       [--verbosity=<level>]"
    print "       [--max-wall-clock=<seconds>]"
    return

#----------------------------------------------------------------------

def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    global verbosity_level
    # Flush before using subprocess to ensure output is in the right order.
    sys.stdout.flush()    
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    if verbosity_level > 0: print "About to run cmd:", string.join(args)
    try:
        result = subprocess.check_call(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print "Subprocess command failed."
        print "    cmd=", args
        raise RuntimeError("Subprocess command failed.")
    return result
    
def quote(str):
    """
    Put quotes around a string.
    """
    return '"' + str + '"'

def get_value_from_ini_file(fileName, parameterName):
    """
    Rummage through an INI file to get a particular parameter value.
    """
    fp = open(fileName,'r'); lines = fp.readlines(); fp.close()
    stringValue = None
    for line in lines:
        if line.startswith(parameterName):
            data = line.strip().split()
            stringValue = data[-1]
            break
    if stringValue is None: 
        raise RuntimeError('Parameter %s not found in file %s' % (parameterName, fileName))
    return stringValue

#---------------------------------------------------------------

def run_in_block_marching_mode(cfgDict):
    """
    The core of the multi-block space-marching code.
    """
    global verbosity_level
    if verbosity_level > 0: print "Set a few overall parameters"
    startWallClockSeconds = time.time()
    jobName = cfgDict['jobName']
    blksPerSlice = cfgDict['nbj'] * cfgDict['nbk']
    max_time = cfgDict['max_time']
    max_dt = cfgDict['max_dt']
    gmodelFile = cfgDict['gmodelFile']
    restartFromRun = cfgDict['restartFromRun']
    # Set up mpirun parameters.
    MPI_PARAMS = "mpirun -np " + str(2 * blksPerSlice) + " "
    if restartFromRun == 0:
        blockDims, numberOfBlks = read_block_dims(jobName+'.config')
    else:
        blockDims, numberOfBlks = read_block_dims('./master/'+jobName+'.config')
    # Compute number of e3mpi runs needed.
    # The first e3mpi run starts with two columns and then 
    # every e3mpi run after that moves along by one column.
    if numberOfBlks % blksPerSlice == 0:
        numberOfEilmer3Runs = numberOfBlks/blksPerSlice - 1
    else:
        print "Oops, we were just computing the number of Eilmer3 runs needed"
        print "    when we found that there was a mismatch in the total"
        print "    number of blocks and the number of blocks per column."
        print "    numberOfBlks=", numberOfBlks
        print "    blksPerSlice=", blksPerSlice
        raise RuntimeError("Invalid data.")
    # Compute the maximum run time (in the simulation universe) for each Eilmer3 run.
    maxRunTime = max_time / numberOfEilmer3Runs
    #
    if restartFromRun == 0:
        if (os.path.exists(jobName+".config") and os.path.exists(jobName+".control") and
            os.path.exists("grid/t0000/") and os.path.exists("flow/t0000/")):
            pass
        else:
            raise RuntimeError("Oops, missing files from prep stage.")
        if os.path.exists("flow/t0001/"):
            raise RuntimeError("Oops, seem to have flow/t0001/ from a previous calculation.")
        if os.path.exists("master"):
            raise RuntimeError("Oops, still have master directory from a previous calculation.")
        run_command('mkdir master')
        run_command('cp %s.config master/%s.config' % (jobName, jobName,))
        run_command('cp %s.control master/%s.control' % (jobName, jobName,))
        run_command('cp block_labels.list master/block_labels.list')
        run_command('mv flow master/')
        run_command('mv grid master/')
        # We will accumulate the converged solution blocks in the master area.
        # Note that the tindx value will be 1 (not 9999).
        run_command('mkdir master/flow/t0001')
        #
        if verbosity_level > 0:
            print "Set up local config files in the usual places for Eilmer3."
        # We are going to assume that the configuration of block subsets
        # will not change run-to-run for the time-marching code.
        run_command('mkdir flow grid')
        run_command('mkdir flow/t0000 grid/t0000')
        create_config_file(blksPerSlice, "master/"+jobName+".config", jobName+".config")
        create_blkLabels_file(blksPerSlice)
        # Update control file to reflect the max_time for each Eilmer3 run.
        update_max_time(maxRunTime, jobName+".control")
    #
    # Start looping for the number of Eilmer3 runs.
    for run in range(restartFromRun,numberOfEilmer3Runs):
        print "-----------------------------------------------"
        print " Eilmer3 in block-marching mode - Run %d of %d. " % (run, numberOfEilmer3Runs-1)
        # For each run, there are a subset of blocks that are integrated in time.
        # This set of blocks is arranged in two columns A (upstream) and B (downstream).
        # Within the run, the blocks have to be labelled locally 0..2*blksPerSlice-1.
        # Globally, these same blocks are labelled firstBlk..lastBlk
        firstBlk = run * blksPerSlice
        lastBlk = firstBlk + (2 * blksPerSlice)
        #
        # Copy the grid and flow files needed for this run.
        # Upstream column (A) 
        for blk in range(firstBlk, firstBlk+blksPerSlice):
            localBlkId = blk - run*blksPerSlice
            if run == 0:
                # Upstream column comes from the master copy on the first run.
                src = 'master/flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
                dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                run_command(['cp', src, dest])
                src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
                dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                run_command(['cp', src, dest])
            elif run != restartFromRun:
                # On subsequent runs, the flow data comes from the downstream column
                # that has most recently been iterated.
                src = 'flow/t0001/'+jobName+'.flow.b'+str(localBlkId+blksPerSlice).zfill(4)+'.t0001.gz'
                dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                run_command(['cp', src, dest])
                # We also need to change the time in the header line back to 0.0 seconds in
                # each file. The files need to be unzipped and zipped again for this process.
                run_command(['gunzip'] + [dest])
                dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000'
                run_command(['sed'] + ['-i'] + ['1s/.*/ 0.000000000000e+00/'] + [dest])
                run_command(['gzip'] + [dest])
            else:
                # We are at the beginning of a restart.
                pass
            # Grid comes always from the master area.
            src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
        # Downstream column (B) always comes freshly from the master copy.
        for blk in range(firstBlk+blksPerSlice, lastBlk):
            localBlkId = blk - run*blksPerSlice
            src = 'master/flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
            src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
        # Also update .config file to reflect dimensions of each new block.
        for blk in range(firstBlk, lastBlk):
            localBlkId = blk - run*blksPerSlice
            update_block_dims(localBlkId, blockDims[blk], jobName+".config")
        if run == 1:
            # Update the inflow boundary conditions if we are in the second run.
            # The inflow condition should not change from this run onwards.
            update_inflowBC(jobName,blksPerSlice, jobName+".config")
        if run > 0:
            # Propagate the dt_global from the previous run. This will hopefully help 
            # achieve faster convergence, but we may have to limit it for stability
            # in case the wall geometry changes suddenly.
            update_dt_global(jobName, max_dt)
            # Propagate the inflow profile across all blocks to generate a starting solution
            # that will help achieve a faster convergence to the steady-state solution.
            # Propagate only the last profile slice of set A to the blocks in set B.
            for blk in range(blksPerSlice, 2*blksPerSlice):
                flowFileName = 'flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
                profileFileName = 'blk-'+str(blk)+'-slice-1.dat'
                propagate_data_west_to_east(flowFileName, profileFileName)
        #
        if verbosity_level > 0: print "Run Eilmer3 time integration process."
        run_command(MPI_PARAMS+E3BIN+('/e3mpi.exe --job=%s --run' % (jobName,)))
        #
        print "Post-process to get profiles for the inflow for the next run."
        # Extract flow profiles consisting of the the last 2 slices for each block 
        # in the upstream column, set A.  These will be used as inflow 
        # boundary conditions in the next stage of the march.
        profile_list_str = "0,east"
        for blk in range(1, blksPerSlice):
            profile_list_str += ';' + str(blk) + ',east'
        run_command(E3BIN+('/e3post.py --job=%s --tindx=0001 ' % (jobName,))
                    +('--static-flow-profile="'+profile_list_str+'" ')
                    +('--gmodel-file=%s' % gmodelFile))
        # Extract the last slice for each block in the downstream column, set B.
        # This will be used to propagate flow data across the block at the start of
        # the next stage of the march.
        for blk in range(blksPerSlice, 2*blksPerSlice):
            run_command(E3BIN+('/e3post.py --job=%s --tindx=0001 ' % (jobName,))
                        +('--output-file=blk-%d-slice-1.dat ' % blk)
                        +('--slice-list="%d,-1,:,:" ' % blk)
                        +('--gmodel-file=%s' % gmodelFile))
        #
        # Save the (presumed) converged blocks (A) in the upstream column back to the master area.
        for blk in range(firstBlk, firstBlk+blksPerSlice):
            localBlkId = blk - run*blksPerSlice
            src = 'flow/t0001/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0001.gz'
            dest = 'master/flow/t0001/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0001.gz'
            run_command(['cp', src, dest])
        if run == numberOfEilmer3Runs-1:
            # On the last run, also save the final column.
            for blk in range(firstBlk+blksPerSlice, lastBlk):
                localBlkId = blk - run*blksPerSlice
                src = 'flow/t0001/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0001.gz'
                dest = 'master/flow/t0001/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0001.gz'
                run_command(['cp', src, dest])
        # Throw away the local grid files.
        run_command(['rm'] + glob('grid/t0000/*'))
        #
        elapsedWallClockSeconds = time.time() - startWallClockSeconds
        print("Elapsed wall-clock seconds = %.1f" % elapsedWallClockSeconds)
        if ( cfgDict['max_wall_clock_seconds'] > 0 and
             elapsedWallClockSeconds > cfgDict['max_wall_clock_seconds'] ):
            print("Stopping part-way because we have exceeded wall-clock time limit")
            print("    finished run = %d" % run)
            print("    elapsedWallClockSeconds = %.1f" % elapsedWallClockSeconds)
            break
        #------------ end of run loop ----------------------------------
    #   
    # Clean up the temporary folders and files and bring the master files back to the
    # usual Eilmer3 working area.
    run_command('rm -r flow grid')
    run_command('mv master/flow flow')
    run_command('mv master/grid grid')
    run_command('mv master/%s.config %s.config' % (jobName, jobName,))
    run_command('mv master/%s.control %s.control' % (jobName, jobName,))
    run_command('mv master/block_labels.list block_labels.list')
    run_command('rm -rf master')
    run_command(['rm'] + glob('blk-*-slice*.dat'))
    run_command(['rm'] + glob('*.static-flow-profile'))
    return

def create_config_file(blksPerSlice, originalConfigFileName, newConfigFileName):
    """
    Create a modified config file to run time-marching simulation on subsets of blocks.
    We do this once because we are going to assume that every subset of blocks is alike, 
    with the same block-to-block connections and the same intra-block discretization.

    Updated to use the ConfigParser module, so that it's more robust. PJ, 07-Sep-2013.
    """
    # We always run with 2 slices, an upstream slice (A) and a downstream slice (B).
    blksPerRun = 2 * blksPerSlice 
    # We will omit the unwanted blocks from the new config file.
    unwantedBlock = 'block/' + str(blksPerRun)
    # Section names for the east faces of blocks in the downstream slice.
    sliceB_section_names = ['block/'+str(blk)+'/face/east' 
                            for blk in range(blksPerSlice, blksPerRun)]
    #
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(originalConfigFileName)
    outfile = file(newConfigFileName, 'w')
    parser.set('global_data', 'nblock', str(blksPerRun))
    for section_name in parser.sections():
        # Assume that sections are stored in the same order as they appear
        # in the original config file.
        if unwantedBlock in section_name: break
        # Blocks in east faces of downstream slice are no longer connected. 
        if section_name in sliceB_section_names:
            parser.set(section_name, 'bc', 'EXTRAPOLATE_OUT')
            parser.set(section_name, 'other_block', '-1')
            parser.set(section_name, 'other_face', 'none')
        outfile.write('['+section_name+']\n')
        for name, value in parser.items(section_name):
            outfile.write('%s = %s\n' % (name, value))
    outfile.write("# end file\n") # Same footer as in original file.
    outfile.close()
    return

def create_blkLabels_file(blksPerSlice):
    """
    Create block_labels.list file to run time-marching simulation on subsets of blocks.
    """
    outfile = file('block_labels.list', 'w')
    outfile.write('# indx label\n')
    for i in range(2 * blksPerSlice):
        outfile.write(str(i) + ' temp-' + str(i) + '\n')
    outfile.close()
    return

def read_block_dims(configFileName):
    """
    Returns the number blocks, and the dimensions for each block.
    """
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(configFileName)
    numberOfBlocks = parser.getint('global_data', 'nblock')
    blockDimsList = [(parser.getint('block/'+str(blk), 'nni'),
                      parser.getint('block/'+str(blk), 'nnj'),
                      parser.getint('block/'+str(blk), 'nnk'))
                     for blk in range(numberOfBlocks)]
    return blockDimsList, numberOfBlocks

def update_block_dims(targetBlock, blockDims, configFileName):
    """
    Update the block dimensions in .config file to reflect that of
    the block that is used in the current Eilmer3 run.
    """
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(configFileName)
    parser.set('block/'+str(targetBlock), 'nni', str(blockDims[0]))
    parser.set('block/'+str(targetBlock), 'nnj', str(blockDims[1]))
    parser.set('block/'+str(targetBlock), 'nnk', str(blockDims[2]))
    outfile = file(configFileName, 'w')
    parser.write(outfile)
    outfile.close()
    return

def update_max_time(maxTime, controlFileName):
    """
    Update the max_time in the control file to reflect the maximum
    run time for each Eilmer3 run.
    """
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(controlFileName)
    parser.set('control_data', 'max_time', str(maxTime))
    outfile = file(controlFileName, 'w')
    parser.write(outfile)
    outfile.close()
    return

def update_dt_global(jobName, max_dt):
    """
    Update dt value in the "jobName.control" file
    based on the last value in "jobName.times" file.

    max_dt is either a string (representing the float value) or None
    """
    # Grab the dt_global value from the last line of the .times file.
    f1 = open(jobName+'.times','r'); data = f1.readlines(); f1.close()
    final_dt_global = float(data[-1].split()[-1])
    if max_dt is None:
        new_dt_global = final_dt_global
    else:
        new_dt_global = min(final_dt_global, float(max_dt))
    # Write out a new .control file with the appropriate line updated.
    controlFileName = jobName+'.control'
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(controlFileName)
    parser.set('control_data', 'dt', str(new_dt_global))
    outfile = file(controlFileName, 'w')
    parser.write(outfile)
    outfile.close()
    return

def update_inflowBC(jobName, blksPerSlice, configFileName):
    """
    Update the inflow boundary condition for upstream-slice (A) blocks
    to use a StaticProfBC().

    The data for this staticProfile has come from a previous run of
    the time-stepping code.
    """
    if sys.version_info < (2, 7):
        parser = SafeConfigParser(dict_type=OrderedDict)
    else:
        parser = SafeConfigParser()
    parser.optionxform = str
    parser.read(configFileName)
    # Section names for the west faces of blocks in the upstream slice.
    sliceA_section_names = ['block/'+str(blk)+'/face/west' 
                            for blk in range(0, blksPerSlice)]
    for section_name in parser.sections():
        if section_name in sliceA_section_names:
            i = sliceA_section_names.index(section_name)
            parser.set(section_name, 'bc', 'STATIC_PROF')
            parser.set(section_name, 'filename', jobName+'-blk-'+str(i)+'-face-east.static-flow-profile')
            parser.set(section_name, 'n_profile', '2')
    outfile = file(configFileName, 'w')
    parser.write(outfile)
    outfile.close()
    return

def propagate_data_west_to_east(flowFileName, profileFileName):
    """
    Extract a profile from a given file and propagate this profile across
    the initial flow solution for the block.
    """
    profile = []
    # Read in file that contains the profile. 
    fi = open(profileFileName, "r")
    fi.readline() # Read away the header line that contains variable names.
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = buf.split()
        profile.append(tokens)
    fi.close()
    # Rewrite the flow solution data with the profile data.
    fi = gzip.open(flowFileName, "r")
    fo = gzip.open('flow/t0000/tmp.gz', "w")
    # Read and write the first three header lines, just as they are.
    buf = fi.readline(); fo.write(buf)
    buf = fi.readline(); fo.write(buf)
    buf = fi.readline(); fo.write(buf)
    # Extract dimensions from the third line.
    tokens = buf.strip().split()
    nni = int(tokens[0]); nnj = int(tokens[1]); nnk = int(tokens[2])
    # Check that the number of cells in the j-direction matches the given profile.
    if nnj*nnk != len(profile):
        print 'The number of cells in the slice ', nnj*nnk,\
              'does not match that of the given profile ', len(profile)
        sys.exit()
    # Read and replace flow data for every cell across the block.
    for k in range(nnk):
        for j in range(nnj):
            for i in range(nni):
                # Read in the current data for the cell.
                tokens = fi.readline().strip().split()
                # Write the pos.x, pos.y, pos.z and vol variables to the new data file.
                fo.write(' ' + ' '.join(tokens[0:4]))
                # Fill out values for the other variables.
                fo.write(' ' + ' '.join(profile[k*nnj+j][4:]) + '\n')
    fi.close(); fo.close()
    # Change the output file name to the inflow file name.
    run_command('mv flow/t0000/tmp.gz ' + flowFileName)
    return

#----------------------------------------------------------------------

def main(uoDict):
    """
    Top-level function for the block-marching application.
    """
    global verbosity_level
    jobName = uoDict.get("--job", "test")
    rootName, ext = os.path.splitext(jobName)
    if os.path.exists(jobName):
        jobFileName = jobName
    else:
        jobFileName = rootName + ".py"
    print "Job file: %s" % jobFileName
    zipFiles = 1  # Default: use zip file format for grid and flow data files.
    if uoDict.has_key("--zip-files"): zipFiles = 1
    if uoDict.has_key("--no-zip-files"): zipFiles = 0
    cfgDict = {'jobName': jobName, 'zipFiles': zipFiles}
    cfgDict['restartFromRun'] = int(uoDict.get("--restart", "0"))
    cfgDict['nbj'] = int(uoDict.get("--nbj", "1"))
    cfgDict['nbk'] = int(uoDict.get("--nbk", "1"))
    cfgDict['max_dt'] = uoDict.get("--max-dt", None) # either a string or None
    cfgDict['max_wall_clock_seconds'] = int(uoDict.get('--max-wall-clock', '-1'))
    # Get some parameters from the files previously written by e3prep.py
    cfgDict['gmodelFile'] = get_value_from_ini_file(jobName+'.config', 'gas_model_file')
    cfgDict['max_time'] = float(get_value_from_ini_file(jobName+'.control', 'max_time'))
    # Do the real work...
    run_in_block_marching_mode(cfgDict)
    # At this time, the full-domain flow solution should be sitting in ./flow/t0001/
    if uoDict.has_key("--polish-time"):
        print "Polish the flow solution by running the whole domain a short time."
        print "    The final solution should be found into ./flow/t0002/"
        print "    Need to set max_time, dt."
        print "    FIX-ME need to finish this code..."
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    verbosity_level = int(uoDict.get("--verbosity", "0"))
    if verbosity_level > 0 or len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        print "Begin e3march.py..."
        print "Source code revision string: ", get_revision_string()
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(1) # abnormal exit; didn't know what to do
    else:
        try:
            main(uoDict)
        except:
            print "This run of e3march.py has gone bad."
            traceback.print_exc(file=sys.stdout)
            sys.exit(2) # abnormal exit; we tried and failed
    if verbosity_level > 0: print "Done."
    sys.exit(0) # Finally, a normal exit.

