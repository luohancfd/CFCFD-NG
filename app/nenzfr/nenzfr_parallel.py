"""
nenzfr_parallel.py -- Functions to coordinate the running of e3mpi for nenzfr.

.. Authors: Wilson Chan, Luke Doherty and Peter Jacobs
"""

import sys, os
E3BIN = os.path.expandvars("$HOME/e3bin")
sys.path.append(E3BIN)
from glob import glob
import gzip
from nenzfr_utils import prepare_input_script, run_command, quote

#---------------------------------------------------------------

def run_in_block_marching_mode(opt, gmodelFile):
    """
    Run nenzfr in multi-block space-marching mode.
    """

    print "-------------------------------------------"
    print " Running nenzfr in block-marching mode ... " 
    print "-------------------------------------------"

    # Assume that the number of blocks per set is the number of radial blocks.
    blksPerSet = opt.nbj

    # Set up mpirun parameters
    MPI_PARAMS = "mpirun -np " + str(2 * blksPerSet) + " "

    # For the multi-block space-marching mode, we initialise the starting
    # solution with the inflow conditions.
    run_command("sed -i 's/fill_condition=initial/fill_condition=inflow/g' %s.py" % opt.jobName)

    # Clean up any remnant files from previous runs 
    run_command(['rm', '-r', 'flow', 'grid', 'temp'] + glob('*.master'))

    # Run e3prep
    run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
    
    # Make a copy of the flow and grid folders, and the config, control
    # and block_labels.list files, and rename to *.master files
    run_command('cp %s.config %s.config.master' % (opt.jobName, opt.jobName,))
    run_command('cp %s.control %s.control.master' % (opt.jobName, opt.jobName,))
    run_command('cp -r block_labels.list block_labels.list.master')
    run_command('mv flow flow.master')
    run_command('mv grid grid.master')
    run_command('cp -r flow.master/t0000 flow.master/intermediate')
    run_command('mkdir flow grid flow/t0000 grid/t0000 flow.master/t9999')
    
    # Modify .config and block_labels.list files, and reading in the
    # number of blocks and the respective dimensions for each block.
    inputFileName = opt.jobName + ".config.master"
    outputFileName = opt.jobName + ".config"
    mod_cfg_file(blksPerSet, inputFileName, outputFileName)
    mod_blkLabels_file(blksPerSet)
    blockDims, noOfBlks = read_block_dims(inputFileName)
    
    # Compute number of runs needed.
    if noOfBlks % blksPerSet == 0:
        noOfEilmer3Runs = noOfBlks / blksPerSet - 1
    else:
        raise RunTimeError("Mismatch in total number of blocks and number of blocks per set.")

    # Read timing file (file that specifies the max_time for each Eilmer3 run).
    inputFileName = opt.jobName + ".timing"
    # Creates timing file if it does not exist.
    if not os.path.exists(inputFileName):
        timingFile = file(inputFileName, 'w')
        time_slice = opt.max_time / noOfEilmer3Runs
        for t in range(noOfEilmer3Runs):
            timingFile.write('%.6e\n' % time_slice)
        timingFile.close()
    # Read timing file.
    maxTimeList = read_timing_file(inputFileName)
    # Ensure that the number of time slices is equal to the number of runs.
    if len(maxTimeList) != noOfEilmer3Runs:
        raise RunTimeError("Number of max_time specified must match the total number of Eilmer3 runs.")

    # Create temporary folder to hold flow and grid files which are to be renamed for each run.
    run_command('mkdir temp')

    # Start looping for the number of Eilmer3 runs.
    firstBlk = 0
    for run in range(noOfEilmer3Runs):
        print "-----------------------------------------------"
        print " nenzfr in block-marching mode - Run %d of %d. " % (run, noOfEilmer3Runs)
        # Identify the last block of this run.
        lastBlk = firstBlk + (2 * blksPerSet)
        # Copy the grid and flow files needed for this run from the 
        # flow.master/intermediate folder.
        for blk in range(firstBlk, lastBlk):
            # Copy to files to temporary folder first.
            targetCmd = ['cp'] + glob('flow.master/intermediate/*b' + str(blk).zfill(4) + '*') + ['temp/'] 
            run_command(targetCmd)
            targetCmd = ['cp'] + glob('grid.master/t0000/*b' + str(blk).zfill(4) + '*') + ['temp/']
            run_command(targetCmd) 
            # From the second Eilmer3 run onwards, ...
            if run != 0:                
                # Change the block numbers of each solution file. Simultaneously
                # change the root name to 'tmp' to avoid over-writing files.
                targetCmd = 'mv temp/' + str(opt.jobName) + '.flow.b' +\
                            str(blk).zfill(4) + '.t0000.gz temp/tmp.flow.b' +\
                            str(blk - run * blksPerSet).zfill(4) + '.t0000.gz'
                run_command(targetCmd)
                targetCmd = 'mv temp/' + str(opt.jobName) + '.grid.b' +\
                            str(blk).zfill(4) + '.t0000.gz temp/tmp.grid.b' +\
                            str(blk - run * blksPerSet).zfill(4) + '.t0000.gz'
                run_command(targetCmd)
            # Also update config file to reflect dimensions of each new block.
            inputFileName = opt.jobName + ".config"
            targetBlock = blk - run * blksPerSet
            update_block_dims(targetBlock, blockDims[blk], inputFileName)
        # ------ PJ FIX-ME--------------------------------------
        # The following rename and mv commands appear malformed.
        # Need to ask Wilson for advice.
        # ------------------------------------------------------
        # Restore the original root name for the files.
        run_command(['rename'] + ['tmp'] + [str(opt.jobName)] + glob('temp/*'))
        # Copy the files from temp folder to the flow & grid folders for this Eilmer3 tun.
        run_command(['mv'] + glob('temp/*flow*') + ['flow/t0000/'])
        run_command(['mv'] + glob('temp/*grid*') + ['grid/t0000/'])
        # Update control file to reflect the max_time for the current run.
        inputFileName = opt.jobName + ".control"
        update_max_time(maxTimeList[run], inputFileName)
        # Update the inflow boundary conditions if we are in the second run.
        # The inflow condition should not change from the second run onwards.
        if run == 1:
            inputFileName = opt.jobName + ".config"
            update_inflowBC(blksPerSet, inputFileName)
        # Propagate the inflow profile across all blocks to generate a starting solution
        # that will help achieve a faster convergence to the steady-state solution.
        if run != 0:
            # Propagate only the last profile slice of set A to the blocks in set B.
            for blk in range(blksPerSet, 2*blksPerSet):
                flowFileName = 'flow/t0000/' + str(opt.jobName) + '.flow.b' +\
                               str(blk).zfill(4) + '.t0000.gz'
                profileFileName = 'blk-' + str(blk - blksPerSet) + \
                                  '-slices-1-and-2.dat'
                propagate_data_west_to_east(flowFileName, profileFileName)
        # Run Eilmer3.
        run_command(MPI_PARAMS+E3BIN+('/e3mpi.exe --job=%s --run' % (opt.jobName,)))
        # Post-process to get profiles for the inflow for the next run.
        for blk in range(0, blksPerSet):
            # Extract the last 2 slices for each block in set A.
            run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 ' % (opt.jobName,))
                        +('--output-file=blk-%d-slices-1-and-2.dat ' % blk)
                        +('--slice-list="%d,-1,:,0;%d,-2,:,0" ' % (blk, blk))
                        +('--gmodel-file=%s' % gmodelFile))
        # Clean up the grid folder.
        run_command(['rm'] + glob('grid/t0000/*'))
        # Copy flow solutions to temporary folder first.
        run_command(['mv'] + glob('flow/t9999/*') + ['temp/'])
        # Then restore block numbers for each solution file to its original block number.
        for blk in range(firstBlk, lastBlk):
            if run != 0:
                # Change the block numbers of each solution file. Simultaneously  
                # change the root name to 'tmp' to avoid over-writing files.
                targetCmd = 'mv temp/' + str(opt.jobName) + '.flow.b' +\
                            str(blk - run * blksPerSet).zfill(4) +\
                            '.t9999.gz temp/tmp.flow.b' + str(blk).zfill(4) + '.t9999.gz'
                run_command(targetCmd)
        # Restore the original root name for the files.
        run_command(['rename'] + ['tmp'] + [str(opt.jobName)] + glob('temp/*'))
        # Copy files from the temp folder back to the flow.master t9999 folder.
        run_command(['cp'] + glob('temp/*flow*') + ['flow.master/t9999/'])
        # Now to use the flow solutions from this run as the starting solution for
        # the next Eilmer3 run, we have to convert the files from t9999 to t0000.
        run_command(['rename'] + ['t9999'] + ['t0000'] + glob('temp/*'))
        # We also need to change the time in the header line back to 0.0 seconds in
        # each file. The files need to be unzipped and zipped again for this process.
        run_command(['gunzip'] + glob('temp/*'))
        run_command(['sed'] + ['-i'] + ['1s/.*/ 0.000000000000e+00/'] + glob('temp/*'))
        run_command(['gzip'] + glob('temp/*'))
        # Move files from the temp folder back to the flow.master intermediate folder.
        run_command(['mv'] + glob('temp/*flow*') + ['flow.master/intermediate/'])
        # Update the index of the first block for the next run.
        firstBlk = firstBlk + blksPerSet
        
    # Clean up the temporary folders and files.
    run_command('rm -r flow grid temp flow.master/intermediate')
    run_command('mv flow.master flow')
    run_command('mv grid.master grid')
    run_command('mv %s.config.master %s.config' % (opt.jobName, opt.jobName,))
    run_command('mv %s.control.master %s.control' % (opt.jobName, opt.jobName,))
    run_command('mv block_labels.list.master block_labels.list')
    return

def mod_cfg_file(blksPerSet, inputFileName, outputFileName):
    """
    Modify config file to allow nenzfr to run in block-marching mode.
    """
    # Initialise some variables and lists.
    lineBuffer = []; Line = []; setBLines = []; lineIndex = 0
    # Read input file
    f = file(inputFileName)
    # We always run with 2 sets, so the number of blocks per run
    # is always twice that of the number of blocks per set.
    blksPerRun = 2 * blksPerSet 
    # Delete the unwanted blocks in the config file and 
    # find the line numbers for the blocks in Sets A and B.
    unwantedBlock = '[block/' + str(blksPerRun) + ']'
    for line in f:
        # Stop appending once we have hit the unwanted blocks.
        if unwantedBlock in line: break
        # Edit to reflect new number of blocks.
        if 'nblock =' in line:
            line = "nblock = " + str(blksPerRun) + "\n"
        # Append line numbers of the target string for Set B. The index for
        # the blocks in Set B always starts from blksPerSet to blksPerRun.
        for blk in range(blksPerSet, blksPerRun):
            targetString = '[block/' + str(blk) + '/face/east]'
            if targetString in line:
                # Append target line numbers.
                setBLines.append(lineIndex)
        # Append lines from the input file to a buffer.
        lineBuffer.append(line)
        # Increase line index counter.
        lineIndex += 1
    # Append the footer that is found in Eilmer3 config files.
    lineBuffer.append("# end file")
    f.close()
    # Start editing the lines in lineBuffer.
    # For Set B, we modify the EAST boundary to reflect an ExtrapolateOutBC.
    for i in range(len(setBLines)):
        lineBuffer[setBLines[i]+2] = "bc = 2\n"
        # "Disconnect" edges of the blocks after changing the boundary condition.
        lineBuffer[setBLines[i]+15] = "other_block = -1\n"
        lineBuffer[setBLines[i]+16] = "other_face = none\n"
    # Write the new config file.
    outfile = file(outputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def mod_blkLabels_file(blksPerSet):
    """
    Modify block_labels.list file to allow nenzfr to run in block-marching mode.
    """
    # Generate list of blocks.
    lineBuffer = []
    for i in range(2 * blksPerSet):
        lineBuffer.append(str(i) + ' temp-' + str(i) + '\n')
    # Write the new block_labels.list file.
    outfile = file('block_labels.list', 'w')
    outfile.write('# indx label\n')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def read_block_dims(inputFileName):
    """
    Find the number blocks, and also read the
    nni, nnj and nnk dimensions for each block.
    """
    blockDims = []; blockIndex = 0
    f = file(inputFileName)
    for line in f:
        # Find number of blocks for the full simulation.
        if 'nblock =' in line:
            tokens = line.split()
            noOfBlks = int(tokens[2])
        # Find the nni, nnj, nnk for each block.
        if 'nni =' in line:
            blockDims.append([])
            tokens = line.split()
            blockDims[blockIndex].append(tokens[2])
        if 'nnj =' in line:
            tokens = line.split()
            blockDims[blockIndex].append(tokens[2])
        if 'nnk =' in line:
            tokens = line.split()
            blockDims[blockIndex].append(tokens[2])
            blockIndex += 1
    f.close()
    return blockDims, noOfBlks

def read_timing_file(inputFileName):
    """
    Read the timing file that contains the max_time for each 
    Eilmer3 run, and return a list of the max_time for each run.
    """
    maxTimeList = []
    f = file(inputFileName)
    for line in f:
        maxTimeList.append(line)
    f.close()
    return maxTimeList

def update_max_time(maxTime, inputFileName):
    """
    Update the max_time in the control file to reflect the max_time 
    that is specified for the current Eilmer3 run.
    """
    lineBuffer = []
    # Read in .control file and update to reflect the new max_time.
    f = file(inputFileName)
    for line in f:
        if 'max_time =' in line:
            line = "max_time = " + maxTime 
        # Append lines from the input file to a buffer.
        lineBuffer.append(line)
    f.close()
    # Write the new config file.
    outfile = file(inputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def update_block_dims(targetBlock, blockDims, inputFileName):
    """
    Update the block dimensions in .control file to reflect that of
    the block that is used in the current Eilmer3 run.
    """
    lineBuffer = []; lineIndex = 0
    targetString = '[block/' + str(targetBlock) + ']'
    # Read in file
    f = file(inputFileName)
    for line in f:
        # Record line index if target string is matched.
        if targetString in line:
            targetLine = lineIndex
        # Append lines from the input file to a buffer.
        lineBuffer.append(line)
        # Increase the line index counter.
        lineIndex += 1
    f.close()
    # Update the nni, nnj and nnk of the target block.
    lineBuffer[targetLine+3] = "nni = " + blockDims[0] + "\n"
    lineBuffer[targetLine+4] = "nnj = " + blockDims[1] + "\n"
    lineBuffer[targetLine+5] = "nnk = " + blockDims[2] + "\n"
    # Write the new config file.
    outfile = file(inputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def update_inflowBC(blksPerSet, inputFileName):
    """
    Update the inflow boundary condition for set A blocks
    to reflect a StaticProfBC().
    """
    lineBuffer = []; setALines = []; lineIndex = 0
    # We always run with 2 sets, so the number of blocks per run
    # is always twice that of the number of blocks per set.
    blksPerRun = 2 * blksPerSet
    # Read in input file.
    f = file(inputFileName)
    for line in f:
        for blk in range(0, blksPerSet):
            targetString = '[block/' + str(blk) + '/face/west]'
            if targetString in line:
                # Append target line numbers.
                setALines.append(lineIndex)
        # Append lines from the input file to a buffer.
        lineBuffer.append(line)
        # Increment line index counter.
        lineIndex += 1
    f.close()
    # For Set A, we modify the WEST boundary to reflect a StaticProfBC.
    for i in range(len(setALines)):
        lineBuffer[setALines[i]+2] = "bc = 10\n"
        lineBuffer[setALines[i]+4] = "filename = blk-" + str(i) + \
                                     "-slices-1-and-2.dat\n"
        lineBuffer[setALines[i]+5] = "n_profile = 2\n"
    # Write the new config file.
    outfile = file(inputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def propagate_data_west_to_east(flowFileName, profileFileName):
    """
    Extract a profile from a given file and propagate this profile across
    the initial flow solution for the specified block.
    """
    profile = []
    # Read in file that contains profile. Note that the input file contains
    # contains 2 profiles (One for the inner ghost cells and the other for 
    # the outer ghost cells). We want only that for the inner ghost cells.
    fi = open(profileFileName, "r")
    fi.readline() # Read away the header line that contains variable names.
    while True:
        buf = fi.readline().strip()
        if len(buf) == 0: break
        tokens = [float(word) for word in buf.split()]
        profile.append(tokens)
    fi.close()
    # Keep only the first profile slice.
    ncells_for_profile = len(profile) / 2
    del profile[ncells_for_profile:]
    # Propagate the profile across the flow solutions.
    fi = gzip.open(flowFileName, "r")
    fo = gzip.open('flow/t0000/tmp.gz', "w")
    # Read and write the first three header lines.
    buf = fi.readline(); fo.write(buf)
    buf = fi.readline(); fo.write(buf)
    buf = fi.readline(); fo.write(buf)
    # Extract dimensions from the third line.
    buf.strip(); tokens = buf.split()
    nni = int(tokens[0]); nnj = int(tokens[1])
    # Check that the number of cells in the j-direction matches the given profile.
    if nnj != len(profile):
        print 'The number of cells in the j-direction ', nnj,\
              'does not match that of the given profile ', len(profile)
        sys.exit()
    # Read and replace flow data.
    for j in range(nnj):
        for i in range(nni):
            # Read in the pos.x, pos.y, pos.z and vol variables first.
            buf = fi.readline().strip(); tokens = buf.split()
            # Write the pos.x, pos.y, pos.z and vol variables to the new data file.
            newline = ' ' + tokens[0] + ' ' + tokens[1] +\
                      ' ' + tokens[2] + ' ' + tokens[3] + ' '
            fo.write(newline)
            # Fill out values for the other variables.
            for variable in range(4, len(profile[0])):
                fo.write("%.12e " % profile[j][variable])
            fo.write("\n")
    fi.close(); fo.close()
    # Change the output file name to the inflow file name.
    run_command('mv flow/t0000/tmp.gz ' + flowFileName)
    return
