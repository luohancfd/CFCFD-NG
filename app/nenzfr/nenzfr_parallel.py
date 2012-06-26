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

restart = "no"
restartFromRun = 126
if restart not in ["yes"]:
    restartFromRun = 0

#---------------------------------------------------------------

def run_in_block_marching_mode(opt, gmodelFile):
    """
    Run nenzfr in multi-block space-marching mode.
    """
    print "-------------------------------------------"
    print " Running nenzfr in block-marching mode ... " 
    print "-------------------------------------------"
    #
    print "Set a few overall parameters"
    jobName = opt.jobName
    # Assume that the number of blocks per set is the number of radial blocks.
    blksPerColumn = opt.nbj
    # Set up mpirun parameters.
    MPI_PARAMS = "mpirun -np " + str(2 * blksPerColumn) + " "
    # For the multi-block space-marching mode, we initialise the starting
    # solution with the inflow conditions in the nenzfr-generated input script.
    run_command("sed -i 's/fill_condition=initial/fill_condition=inflow/g' %s.py" % jobName)
    #
    if restart not in ["yes"]:
        # Clean up any remnant files from previous runs
        run_command('rm -rf flow grid master')
        #
        print "Run e3prep to generate a collection of blocks for the full nozzle."
        run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % jobName))
        blockDims, noOfBlks = read_block_dims(jobName+'.config')
    elif restart in ["yes"]:
        blockDims, noOfBlks = read_block_dims('./master/'+jobName+'.config')
    # Compute number of e3mpi runs needed.
    # The first e3mpi run starts with two columns and then 
    # every e3mpi run after that moves along by one column.
    if noOfBlks % blksPerColumn == 0:
        noOfEilmer3Runs = noOfBlks/blksPerColumn - 1
    else:
        raise RuntimeError("Mismatch in total number of blocks and number of blocks per column.")
    # Compute the maximum run time for each Eilmer3 run.
    maxRunTime = opt.max_time / noOfEilmer3Runs
    #
    if restart not in ["yes"]:
        print "Set up the master copy of the blocks and files."
        run_command('mkdir master')
        run_command('cp %s.config master/%s.config' % (jobName, jobName,))
        run_command('cp %s.control master/%s.control' % (jobName, jobName,))
        run_command('cp -r block_labels.list master/block_labels.list')
        run_command('mv flow master/')
        run_command('mv grid master/')
        # We will accumulate the converged solution blocks in the master area.
        run_command('mkdir master/flow/t9999')
        #
        print "Set up current-run config files in the usual places for Eilmer3."
        run_command('mkdir flow grid')
        run_command('mkdir flow/t0000 grid/t0000')
        # Modify .config and block_labels.list files for the run of a subset of blocks.
        mod_cfg_file(blksPerColumn, "master/"+jobName+".config", jobName+".config")
        mod_blkLabels_file(blksPerColumn)
        # Update control file to reflect the max_time for each Eilmer3 run.
        update_max_time(maxRunTime, jobName+".control")

    # Start looping for the number of Eilmer3 runs.
    for run in range(restartFromRun,noOfEilmer3Runs):
        print "-----------------------------------------------"
        print " nenzfr in block-marching mode - Run %d of %d. " % (run, noOfEilmer3Runs-1)
        # For each run, there are a subset of blocks that are integrated in time.
        # This set of blocks is arranged in two columns A (upstream) and B (downstream).
        # Within the run, the blocks have to be labelled locally 0..2*blksPerColumn-1.
        # Globally, these same blocks are labelled firstBlk..lastBlk
        firstBlk = run * blksPerColumn
        lastBlk = firstBlk + (2 * blksPerColumn)
        #
        # Copy the grid and flow files needed for this run.
        #
        for blk in range(firstBlk, firstBlk+blksPerColumn):
            localBlkId = blk - run*blksPerColumn
            if run == 0:
                # Upstream column comes from the master copy on the first run.
                src = 'master/flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
                dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                run_command(['cp', src, dest])
                src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
                dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                run_command(['cp', src, dest])
            else:
                if run != restartFromRun:
                    # On subsequent runs, the flow data comes from the downstream column
                    # that has most recently been iterated.
                    src = 'flow/t9999/'+jobName+'.flow.b'+str(localBlkId+blksPerColumn).zfill(4)+'.t9999.gz'
                    dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
                    run_command(['cp', src, dest])
                    # We also need to change the time in the header line back to 0.0 seconds in
                    # each file. The files need to be unzipped and zipped again for this process.
                    run_command(['gunzip'] + [dest])
                    dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000'
                    run_command(['sed'] + ['-i'] + ['1s/.*/ 0.000000000000e+00/'] + [dest])
                    run_command(['gzip'] + [dest])
            # Grid comes always from the master area.
            src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
        # Downstream column always comes freshly from the master copy.
        for blk in range(firstBlk+blksPerColumn, lastBlk):
            localBlkId = blk - run*blksPerColumn
            src = 'master/flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'flow/t0000/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
            src = 'master/grid/t0000/'+jobName+'.grid.b'+str(blk).zfill(4)+'.t0000.gz'
            dest = 'grid/t0000/'+jobName+'.grid.b'+str(localBlkId).zfill(4)+'.t0000.gz'
            run_command(['cp', src, dest])
        # Also update .config file to reflect dimensions of each new block.
        for blk in range(firstBlk, lastBlk):
            localBlkId = blk - run*blksPerColumn
            update_block_dims(localBlkId, blockDims[blk], jobName+".config")
        if run == 1:
            # Update the inflow boundary conditions if we are in the second run.
            # The inflow condition should not change from this run onwards.
            update_inflowBC(blksPerColumn, jobName+".config")
        if run > 0:
            # Propagate the dt_global from the previous run. This will hopefully help 
            # achieve faster convergence.
            update_dt_global(jobName)
            # Propagate the inflow profile across all blocks to generate a starting solution
            # that will help achieve a faster convergence to the steady-state solution.
            # Propagate only the last profile slice of set A to the blocks in set B.
            for blk in range(blksPerColumn, 2*blksPerColumn):
                flowFileName = 'flow/t0000/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t0000.gz'
                profileFileName = 'blk-'+str(blk-blksPerColumn)+'-slices-1-and-2.dat'
                propagate_data_west_to_east(flowFileName, profileFileName)
        #
        print "Run Eilmer3."
        run_command(MPI_PARAMS+E3BIN+('/e3mpi.exe --job=%s --run' % (jobName,)))
        #
        print "Post-process to get profiles for the inflow for the next run."
        # Extract the last 2 slices for each block in the upstream column.
        for blk in range(0, blksPerColumn):
            run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 ' % (jobName,))
                        +('--output-file=blk-%d-slices-1-and-2.dat ' % blk)
                        +('--slice-list="%d,-1,:,0;%d,-2,:,0" ' % (blk, blk))
                        +('--gmodel-file=%s' % gmodelFile))
        #
        # Save the (presumed) converged blocks in the upstream column back to the master area.
        for blk in range(firstBlk, firstBlk+blksPerColumn):
            localBlkId = blk - run*blksPerColumn
            src = 'flow/t9999/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t9999.gz'
            dest = 'master/flow/t9999/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t9999.gz'
            run_command(['cp', src, dest])
        if run == noOfEilmer3Runs-1:
            # On the last run, also save the final column.
            for blk in range(firstBlk+blksPerColumn, lastBlk):
                localBlkId = blk - run*blksPerColumn
                src = 'flow/t9999/'+jobName+'.flow.b'+str(localBlkId).zfill(4)+'.t9999.gz'
                dest = 'master/flow/t9999/'+jobName+'.flow.b'+str(blk).zfill(4)+'.t9999.gz'
                run_command(['cp', src, dest])
        # Throw away the local grid files.
        run_command(['rm'] + glob('grid/t0000/*'))
        #------------ end of run loop ----------------------------------
        
    # Clean up the temporary folders and files and bring the master files back to the
    # usual Eilmer3 working area.
    run_command('rm -r flow grid')
    run_command('mv master/flow flow')
    run_command('mv master/grid grid')
    run_command('mv master/%s.config %s.config' % (jobName, jobName,))
    run_command('mv master/%s.control %s.control' % (jobName, jobName,))
    run_command('mv master/block_labels.list block_labels.list')
    run_command('rm -rf master')
    return

def mod_cfg_file(blksPerColumn, inputFileName, outputFileName):
    """
    Modify config file to allow nenzfr to run in block-marching mode.
    """
    # Initialise some variables and lists.
    lineBuffer = []; Line = []; setBLines = []; lineIndex = 0
    # Read input file
    f = file(inputFileName)
    # We always run with 2 sets, so the number of blocks per run
    # is always twice that of the number of blocks per set.
    blksPerRun = 2 * blksPerColumn 
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
        # the blocks in Set B always starts from blksPerColumn to blksPerRun.
        for blk in range(blksPerColumn, blksPerRun):
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

def mod_blkLabels_file(blksPerColumn):
    """
    Modify block_labels.list file to allow nenzfr to run in block-marching mode.
    """
    # Generate list of blocks.
    lineBuffer = []
    for i in range(2 * blksPerColumn):
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

def update_max_time(maxTime, inputFileName):
    """
    Update the max_time in the control file to reflect the maximum
    run time for each Eilmer3 run.
    """
    lineBuffer = []
    # Read in .control file and update to reflect the new max_time.
    f = file(inputFileName)
    for line in f:
        if 'max_time =' in line:
            line = "max_time = " + str(maxTime) + "\n"
        # Append lines from the input file to a buffer.
        lineBuffer.append(line)
    f.close()
    # Write the new .control file.
    outfile = file(inputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def update_dt_global(jobName):
    """
    Update dt value in the "jobName.control" file
    based on the last value in "jobName.times" file.
    """

    # Read in data from .times file
    f1 = open(jobName+'.times','r')
    data = f1.readlines()
    f1.close()

    # Grab the dt_global value from the last line of the file.
    # We leave it as a string for convenience
    final_dt_global = data[-1].split()[-1]

    # Now read in the data from the .control file
    f2 = open(jobName+'.control','r')
    data2 = f2.readlines()
    f2.close()

    # Write out a new .control file with the appropriate
    # line updated
    fout = open(jobName+'.control','w')
    for line in data2:
        if line.split()[0] == 'dt':
            linedata = line.split()
            linedata[-1] = final_dt_global
            newline = ''.join([' '.join(linedata),'\n'])
            fout.write(newline)
        else:
            fout.write(line)
    fout.close()

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

def update_inflowBC(blksPerColumn, inputFileName):
    """
    Update the inflow boundary condition for set A blocks
    to reflect a StaticProfBC().
    """
    lineBuffer = []; setALines = []; lineIndex = 0
    # We always run with 2 sets, so the number of blocks per run
    # is always twice that of the number of blocks per set.
    blksPerRun = 2 * blksPerColumn
    # Read in input file.
    f = file(inputFileName)
    for line in f:
        for blk in range(0, blksPerColumn):
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
