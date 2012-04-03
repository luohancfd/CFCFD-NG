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
   School of Mechancial and Mining Engineering
   The University of Queensland
"""

VERSION_STRING = "01-Mar-2012"

import shlex, subprocess, string
from subprocess import PIPE
from string import upper
import sys, os, gzip
import optparse
from glob import glob
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
    print "About to run cmd:", cmdText
    if (type(cmdText) is list):
        args = cmdText
    else:
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

def print_stats(sliceFileName,jobName):
    """
    Display statistics of flow properties at the nozzle exit.
    """
    print "Nozzle-exit statistics:"
    fp = open(sliceFileName, 'r')
    # Keep a list of variables in order of appearance.
    varLine = fp.readline().strip()
    items = varLine.split()
    if items[0] == '#': del items[0]
    if items[0] == 'Variables:': del items[0]
    variable_list = [item.split(':')[1] for item in items]
    # print "variable_list=", variable_list
    # Store the data in lists against these names.
    data = {}
    for var in variable_list:
        data[var] = []
    for line in fp.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        assert len(items) == len(variable_list)
        for i in range(len(items)):
            data[variable_list[i]].append(float(items[i]))
    fp.close()
    #
    # Identify edge of core flow.
    ys = data['pos.y']
    y_edge = ys[-1] * 2.0/3.0  # somewhat arbitrary...
    #
    # Compute and print area-weighted-average core flow values.
    exclude_list = ['pos.x', 'pos.y', 'pos.z', 'volume', 'vel.z', 'S']
    #
    fout = open(jobName+'-exit.stats','w')
    fout.write('%10s  %12s   %10s  %10s %10s\n' % \
                   ("variable","mean-value","minus","plus","std-dev"))
    fout.write(60*'-')
    fout.write('\n')
    #
    print "%10s  %12s    %10s %10s %10s" % \
        ("variable", "mean-value", "minus", "plus","std-dev")
    print 60*'-'
    for var in variable_list:
        if var in exclude_list: continue
        A = 0.0; F = 0.0;
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            if j == 0:
                y0 = 0.0
            else:
                y0 = 0.5*(ys[j-1]+ys[j])
            y1 = 0.5*(ys[j]+ys[j+1])
            dA = y1**2 - y0**2
            F += data[var][j] * dA
            A += dA
        mean = F/A
        # Identify low and high values.
        diff_minus = 0.0
        diff_plus = 0.0
        count = 0.0
        stddev = 0.0
        for j in range(len(ys)):
            if ys[j] > y_edge: break
            diff = data[var][j] - mean
            diff_minus = min(diff, diff_minus)
            diff_plus = max(diff, diff_plus)
            count += 1
            stddev += diff**2
        # Calculate the sample standard deviation
        stddev = (stddev/(count-1))**0.5
        print "%10s  %12.4g    %10.3g %10.3g %10.3g" % \
              (var, mean, diff_minus, diff_plus, stddev)
        fout.write('%10s  %12.4g    %10.3g %10.3g %10.3g\n' % \
              (var, mean, diff_minus, diff_plus, stddev))
    #    
    print 60*'-'
    #
    fout.write(60*'-')
    fout.close()
    return

def run_in_block_marching_mode(opt, gmodelFile):
    """
    Run Eilmer3 in multi-block space-marching mode.
    """

    print "-------------------------------------------"
    print " Running nenzfr in block-marching mode ... " 
    print "-------------------------------------------"

    # Assume that the number of blocks per set is the number of radial blocks.
    blksPerSet = opt.nbj

    # Set up mpirun parameters
    MPI_PARAMS = "mpirun -np " + str(2 * blksPerSet) + " "

    # Step 1: Run e3prep
    run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
    
    # Step 2: Make a copy of the flow and grid folders and 
    #         the config and control files, and rename to *master
    run_command('cp %s.config %s.config.master' % (opt.jobName, opt.jobName,))
    run_command('cp %s.control %s.control.master' % (opt.jobName, opt.jobName,))
    run_command('cp -r block_labels.list block_labels.list.master')
    run_command('mv flow flow.master')
    run_command('mv grid grid.master')
    run_command('cp -r flow.master/t0000 flow.master/intermediate')
    run_command('mkdir flow grid flow/t0000 grid/t0000 flow.master/t9999')
    
    # Step 3: Modify .config and block_labels.list files, and reading in the
    #         number of blocks and the respective dimensions for each block.
    inputFileName = opt.jobName + ".config.master"
    outputFileName = opt.jobName + ".config"
    mod_cfg_file(blksPerSet, inputFileName, outputFileName)
    mod_blkLabels_file(blksPerSet)
    blockDims, noOfBlks = read_block_dims(inputFileName)
    
    # Step 4: Compute number of runs needed.
    if noOfBlks % blksPerSet == 0:
        # The number of runs is the number of sets minus 1.
        noOfEilmer3Runs = noOfBlks / blksPerSet - 1
    else:
        # Abort the calculation if there's a mismatch in the total number 
        # of blocks and the number of blocks per set.
        print "There's a mismatch in the total number of blocks \
               and the number of blocks per set."
        return -2

    # Step 5: Read timing file (file that states the max_time for each eilmer3 run)
    inputFileName = opt.jobName + ".timing"
    maxTimeList = read_timing_file(inputFileName)
    # TODO -> Unfortunately, because the number of blocks in the throat region is
    #         not known when setting up the nenzfr calculation (it is determined
    #         automatically), the user might not get the number of lines in the
    #         timing file correct in his/her first go.
    if len(maxTimeList) != noOfEilmer3Runs:
        print "Error: The number of max_time specified must match \
               the total number of Eilmer3 runs."
        return -2

    # Step 6: Start looping for the number of Eilmer3 runs.
    firstBlk = 0
    # Create temporary folder to hold flow and grid files which are to be renamed for each run.
    run_command('mkdir temp')
    for run in range(noOfEilmer3Runs):
        print "-----------------------------------------------"
        print " nenzfr in block-marching mode - Run %d of %d. " % (run, noOfEilmer3Runs)
        # Identify last block of this run.
        lastBlk = firstBlk + (2 * blksPerSet)
        # Copy the grid and flow files needed for this run from the flow.master folder.
        for blk in range(firstBlk, lastBlk):
            # Copy to temporary folder first.
            targetCmd = ['cp'] + glob('flow.master/intermediate/*b' + str(blk).zfill(4) + '*')\
                        + ['temp/'] 
            run_command(targetCmd)
            targetCmd = ['cp'] + glob('grid.master/t0000/*b' + str(blk).zfill(4) + '*')\
                        + ['temp/']
            run_command(targetCmd)
            # then modify the file names.
            if run != 0:                
                # Change the block numbers of each solution file and move them to
                # temporary files first to avoid over-writing files.
                targetCmd = 'mv temp/' + str(opt.jobName) + '.flow.b' +\
                            str(blk).zfill(4) + '.t0000.gz temp/tmp.flow.b' +\
                            str(blk - run * blksPerSet).zfill(4) + '.t0000.gz'
                run_command(targetCmd)
                targetCmd = 'mv temp/' + str(opt.jobName) + '.grid.b' +\
                            str(blk).zfill(4) + '.t0000.gz temp/tmp.grid.b' +\
                            str(blk - run * blksPerSet).zfill(4) + '.t0000.gz'
                run_command(targetCmd)
            # Update config file to reflect dimensions of the new blocks.
            inputFileName = opt.jobName + ".config"
            targetBlock = blk - run * blksPerSet
            update_block_dims(targetBlock, blockDims[blk], inputFileName)
        # Rename temporary files back to original job name.
        targetCmd = ['rename'] + ['tmp'] + [str(opt.jobName)] + glob('temp/*')
        run_command(targetCmd)
        # Copy the files to the flow and grid folders for this Eilmer3 tun.
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
        # Run Eilmer3.
        run_command(MPI_PARAMS+E3BIN+('/e3mpi.exe --job=%s --run' % (opt.jobName,)))
        # Post-process to get profiles for the inflow for the next run.
        for blk in range(0, blksPerSet):
            # Take the last 2 slices for each block in set A.
            run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 ' % (opt.jobName,))
                        +('--output-file=blk-%d-slice-1.dat ' % blk)
                        +('--slice-list="%d,-1,:,0" --gmodel-file=%s' % (blk, gmodelFile)))
            run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 ' % (opt.jobName,))
                        +('--output-file=blk-%d-slice-2.dat ' % blk)
                        +('--slice-list="%d,-2,:,0" --gmodel-file=%s' % (blk, gmodelFile)))
        # Copy flow solutions to temporary folder first.
        run_command(['mv'] + glob('flow/t9999/*') + ['temp/'])
        # Clean up the grid folder.
        run_command(['rm'] + glob('grid/t0000/*'))
        # Restore flow files to its original name and copy them back to the master folder.
        for blk in range(firstBlk, lastBlk):
            # Restore original file names
            if run != 0:
                # Change the block numbers of each solution file and move them to 
                # temporary files first to avoid over-writing files.
                targetCmd = 'mv temp/' + str(opt.jobName) + '.flow.b' +\
                            str(blk - run * blksPerSet).zfill(4) +\
                            '.t9999.gz temp/tmp.flow.b' + str(blk).zfill(4) + '.t9999.gz'
                run_command(targetCmd)
        # Rename temporary files back to original job name.
        targetCmd = ['rename'] + ['tmp'] + [str(opt.jobName)] + glob('temp/*')
        run_command(targetCmd)
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
        # Move files from the temp folder back to the flow.master intermediate folder
        run_command(['mv'] + glob('temp/*flow*') + ['flow.master/intermediate/'])
        # Update the index of the first block for the next run
        firstBlk = firstBlk + blksPerSet
        #
    # Step 7: Clean up the temporary folders and files.
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
    # find the line numbers for the blocks in sets A and B
    unwantedBlock = '[block/' + str(blksPerRun) + ']'
    for line in f:
        # Stop appending once we have hit the unwanted blocks
        if unwantedBlock in line: break
        # Edit to reflect new number of blocks
        if 'nblock =' in line:
            line = "nblock = " + str(blksPerRun) + "\n"
        # Append line numbers of the target string for set B. The index for the
        # blocks in set B always starts from blksPerSet to blksPerRun.
        for blk in range(blksPerSet, blksPerRun):
            targetString = '[block/' + str(blk) + '/face/east]'
            if targetString in line:
                # Append target line numbers
                setBLines.append(lineIndex)
        # Append lines from the input file to a buffer
        lineBuffer.append(line)
        # Increase line index counter
        lineIndex += 1
    # Append the footer that is found in eilmer3 config files
    lineBuffer.append("# end file")
    f.close()
    # Start editing the lines in lineBuffer.
    # For set B, we modify the EAST boundary to reflect an ExtrapolateOutBC    
    for i in range(len(setBLines)):
        lineBuffer[setBLines[i]+2] = "bc = 2\n"
        # "Disconnect" edges of the blocks after changing the boundary condition.
        lineBuffer[setBLines[i]+13] = "other_block = -1\n"
        lineBuffer[setBLines[i]+14] = "other_face = none\n"
    # Write the new config file
    outfile = file(outputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return

def mod_blkLabels_file(blksPerSet):
    """
    Modify block_labels.list file to allow nenzfr to run in block-marching mode.
    """
    # Generate list of blocks
    lineBuffer = []
    for i in range(2 * blksPerSet):
        lineBuffer.append(str(i) + ' temp-' + str(i) + '\n')
    # Write the new block_labels.list file
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
        # Find number of blocks for the full simulation
        if 'nblock =' in line:
            tokens = line.split()
            noOfBlks = int(tokens[2])
        # Find the nni, nnj, nnk for each block
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
    # Read in .control file and update to reflect the new max_time
    f = file(inputFileName)
    for line in f:
        if 'max_time =' in line:
            line = "max_time = " + maxTime + "\n"
        # Append lines from the input file to a buffer
        lineBuffer.append(line)
    f.close()
    # Write the new config file
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
        # Record line index if target string is matched
        if targetString in line:
            targetLine = lineIndex
        # Append lines from the input file to a buffer
        lineBuffer.append(line)
        # Increase the line index counter
        lineIndex += 1
    f.close()
    # Update the nni, nnj and nnk of the target block
    lineBuffer[targetLine+3] = "nni = " + blockDims[0] + "\n"
    lineBuffer[targetLine+4] = "nnj = " + blockDims[1] + "\n"
    lineBuffer[targetLine+5] = "nnk = " + blockDims[2] + "\n"
    # Write the new config file
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
    # is always twice that of the number of blocks per set
    blksPerRun = 2 * blksPerSet
    # Read in input file
    f = file(inputFileName)
    for line in f:
        for blk in range(0, blksPerSet):
            targetString = '[block/' + str(blk) + '/face/west]'
            if targetString in line:
                # Append target line numbers
                setALines.append(lineIndex)
        # Append lines from the input file to a buffer
        lineBuffer.append(line)
        # Increment line index counter
        lineIndex += 1
    f.close()
    # For set A, we modify the WEST boundary to reflect a StaticProfBC
    # TODO -> Need to change this to accept new BC class that accepts 2
    #         slices of input profiles.
    for i in range(len(setALines)):
        lineBuffer[setALines[i]+2] = "bc = 10\n"
        lineBuffer[setALines[i]+4] = "filename = blk-" + str(i) + \
                                     "-slice-1.dat\n"
    # Write the new config file
    outfile = file(inputFileName, 'w')
    outfile.writelines(lineBuffer)
    outfile.close()
    return


def main():
    """
    Examine the command-line options to decide the what to do
    and then coordinate the calculations done by Eilmer3.
    """
    op = optparse.OptionParser(version=VERSION_STRING)
    op.add_option('--gas', dest='gasName', default='air',
                  choices=['air', 'air5species', 'n2', 'co2', 'h2ne'],
                  help=("name of gas model: "
                        "air; " "air5species; " "n2; " "co2; " "h2ne"))
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
                  default='contour-t4-m10.data',
                  help="file containing nozzle contour [default: %default]")
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
    # The following defaults suit Like's Mach 10 calculations.
    op.add_option('--nni', dest='nni', type='int', default=1800,
                  help=("number of axial cells"))
    op.add_option('--nnj', dest='nnj', type='int', default=100,
                  help=("number of radial cells"))
    op.add_option('--nbi', dest='nbi', type='int', default=180,
                  help=("number of axial blocks for the divergence section (nozzle_blk)"))
    op.add_option('--nbj', dest='nbj', type='int', default=1,
                  help=("number of radial blocks"))
    op.add_option('--bx', dest='bx', type='float', default=1.05,
                  help=("clustering in the axial direction"))
    op.add_option('--by', dest='by', type='float', default=1.002,
                  help=("clustering in the radial direction"))
    op.add_option('--max-time', dest='max_time', type='float', default=6.0e-3,
                  help=("overall simulation time for nozzle flow"))
    op.add_option('--max-step', dest='max_step', type='int', default=80000,
                  help=("maximum simulation steps allowed"))
    
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
    opt, args = op.parse_args()
    #
    # If we have already run a calculation, it may be that we just want
    # to extract the exit-flow statistics again.
    if opt.justStats:
        print_stats(opt.exitSliceFileName,opt.jobName)
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
    # Get the nozzle contour file into the current work directory.
    run_command('cp '+E3BIN+'/nenzfr_data_files/'+opt.contourFileName+' .')
    # Set up the equilibrium gas-model file as a look-up table.
    if opt.chemModel in ['eq',]:
        if opt.gasName in ['n2']:
            eqGasModelFile = 'cea-lut-'+upper(opt.gasName)+'-no-ions.lua.gz'
        else:
            eqGasModelFile = 'cea-lut-'+opt.gasName+'-no-ions.lua.gz'
        if not os.path.exists(eqGasModelFile):
            run_command('build-cea-lut --case='+opt.gasName+'-no-ions --extrapolate')
        gmodelFile = eqGasModelFile
    else:
        # We'll assume that the gas-model file of default name is set up.
        # TODO: Luke, this needs to be modified, I suspect.
        gmodelFile = 'gas-model.lua'
    # Runs estcj to get the equilibrium shock-tube conditions up to the nozzle-supply region.
    command_text = E3BIN+('/estcj.py --gas=%s --T1=%g --p1=%g --Vs=%g --pe=%g --task=st --ofn=%s' % 
                          (opt.gasName, opt.T1, opt.p1, opt.Vs, opt.pe, opt.jobName))
    run_command(command_text)
    # Switch off block-sequencing flag if we have selected to run in MPI block-marching mode.
    if opt.blockMarching: 
        sequenceBlocksFlag = 0
    else: 
        sequenceBlocksFlag = 1
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
    #
    # Run Eilmer3 either in the single processor mode, or the
    # multi-processor block-marching mode
    if opt.blockMarching:
        # Assume that the number of blocks per set to 
        # be the same as the number of radial blocks
        run_in_block_marching_mode(opt, gmodelFile)
    else:
        run_command(E3BIN+('/e3prep.py --job=%s --do-svg' % (opt.jobName,)))
        run_command(E3BIN+('/e3shared.exe --job=%s --run' % (opt.jobName,)))
    #
    # Exit plane slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s --slice-list="-1,-2,:,0" ' % 
                  (opt.exitSliceFileName,))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Centerline slice
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +('--output-file=%s-centreline.data --slice-list=":,:,1,0" ' % 
                  (opt.jobName,))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Prep files for plotting with Paraview            
    run_command(E3BIN+('/e3post.py --job=%s --vtk-xml --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
                +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
    # Generate averaged exit flow properties
    print_stats(opt.exitSliceFileName,opt.jobName)                
    # Compute viscous data at the nozzle wall
    run_command(E3BIN+'/nenzfr_compute_viscous_data.py --job=%s' % (opt.jobName,))
     
    #
    # The following are additional commands specific to Luke D. and the Mach 10 nozzle.
    # TODO: Luke, we need to think of a good way to encode all of the options.
    #       I've added quite a few command-line options and more are needed.
    #       Maybe a dictionary of customised parameters for each case.
    #
    # Extract a slice from the last block along jk index directions at the i-index that
    # is closest to the point x=1.642, y=0.0, z=0.0
    run_command(E3BIN+('/e3post.py --job=%s --tindx=9999 --gmodel-file=%s ' % 
                       (opt.jobName, gmodelFile))
               +('--output-file=%s2 --slice-at-point="-1,jk,1.642,0,0" ' % 
                 (opt.exitSliceFileName,))
               +'--add-mach --add-pitot --add-total-enthalpy --add-total-p')
     
    # 29-07-2011 Luke D. 
    # Copy the original exit stats file to a temporary file
    run_command(('cp %s-exit.stats %s-exit.stats_temp') % (opt.jobName, opt.jobName,))
    # (Re) Generate the exit stats using the slice-at-point data
    print_stats(opt.exitSliceFileName+'2',opt.jobName)
    run_command('cp %s-exit.stats %s-exit.stats2' % (opt.jobName, opt.jobName,))
    # Now rename the temporary exit stats file back to its original name
    run_command('mv %s-exit.stats_temp %s-exit.stats' % (opt.jobName, opt.jobName,))
    #
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
