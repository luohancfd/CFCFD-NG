#! /usr/bin/env python
"""
Python program to pick up the data after a simulation.

e3post.py is the principal post-processor for slicing and dicing your flow data.
There are so many things that you might want to do with your data from
a simulation that this program is structured as a library of functions 
to pick up the data, maybe add a few items, and then write the data to some 
other format.  Invoking it with the --help option will present the summary of options.

Usage
----- 

Command line::

  e3post.py [options]

Summary of options::

| e3post.py [--help] [--job=<jobFileName>] [--tindx=<index|all>]
|           [--zip-files|--no-zip-files]
|           [--moving-grid]
|           [--omegaz="[omegaz0,omegaz1,...]"]
| 
|           [--add-pitot-p] [--add-total-p] [--add-mach] [--add-total-enthalpy]
|           [--add-molef --gmodel-file="gas-model.lua"]
|           [--add-transport-coeffs --gmodel-file="gas-model.lua"]
|           [--add-user-computed-vars="user-script.py"]
| 
|           [--vtk-xml] [--binary-format] [--tecplot] [--plot3d] [--OpenFoam]
| 
|           [--output-file=<profile-data-file>]
|           [--slice-list="blk-range,i-range,j-range,k-range;..."]
|           [--slice-at-point="blk-range,index-pair,x,y,z;..."]
|           [--slice-along-line="x0,y0,z0,x1,y1,z1,N"]
|           [--surface-list="blk,surface-name;..."]
|           [--local-surface-list="blk,surface-name;..."]
|           [--static-flow-profile="blk,face-name"]
| 
|           [--heat-flux-list="blk-range,surf-range,i-range,j-range,k-range;..."]
|           [--tangent-slab-list="blk-range,i-range,j-range,k-range;..."]
| 
|           [--probe="x,y,z;..."]
| 
|           [--report-norms]
|           [--per-block-norm-list="jb,var-name,norm-name;..."
|           [--global-norm-list="var-name,norm-name;..."
|           [--ref-function=<python-script>]
|           [--compare-job=<jobFileName> [--compare-tindx=<index>]]
| 
|           [--prepare-restart] [--prepare-fstc-restart]
|           [--put-into-folders]
|
|           [--verbosity=<int>]

Examples
--------

* Extract the final solution frame and write VTK files for Paraview::

    e3post.py --job=cone20 --vtk-xml 

* Extract a particular solution frame::

    e3post.py --job=n90 --tindx=5 --vtk-xml --binary-format

* Extract a slice of a particular solution frame::

    e3post.py --job=n90 --output-file=n90_100_iy1.data --tindx=5 \\
        --slice-list="0,:,1,0"

* Compare one solution frame with another::

    e3post.py --job=euler_manufactured --tindx=6 \\
        --compare-job=euler_manufactured --compare-tindx=20

* Compare a solution frame with data provided by a function and write
  particular norms::

    e3post.py --job=euler_manufactured --tindx=20 \\
        --ref-function=euler_wrapper.py \\
        --per-block-norm-list="0,rho,L2;0,rho,L1" \\
        --global-norm-list="rho,L2" \\


Notes
-----

* slice-list:
  Several slices (separated by semicolons) may be specified 
  in the one string.  Each slice specification consists of 
  4 indices or index ranges separated by commas.  
  An index is a single integer value.  An index range may be 
  a colon-separated pair of integers, a colon and one limit 
  or just a colon by itself (to indicate the full range).
  The range limits are inclusive.
* slice-at-point:
  The index-pair is one of ij, jk or ki.  
  The program sets these indices to zero and searches along 
  the remaining index to find the cell nearest the specified 
  (x,y,z) point.  Once found, the slice over the index pair 
  is selected for output (by adding it to the slice-list.
  Beware that, for each block selected, slice-at-point will 
  always select a slice to output, even if it is not very close. 
* Note that you must use double-quotes to prevent the 
  command shell from pulling the string apart.
* add-pitot-p, add-total-p and add-mach work for writing profile files,
  surface (VTK) files, TECPLOT files, VTK-XML files and Plot3D files.
* When choosing Plot3D output, two grid files are generated.  The first,
  with .grd extension, is the true grid as used by the simulation with
  mesh location at the nodes.  The second, with extension .g, has
  cell-centred values and accompanies the cell-centred values in 
  the .f file.
* The angular velocities of the rotating blocks are written as 
  a list with Python syntax.
* The addition of mole-fractions needs to be done in the context
  of a valid gas model (because we need species and molecular masses).

.. Author: P.Jacobs and many others

.. Versions:
   19-March-2008 initial code
   Look at the revision control system logs for the many additions. 
"""

# ----------------------------------------------------------------------
#
import sys
import os
import ConfigParser
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
sys.path.append("") # so that we can find user's scripts in current directory
from getopt import getopt, GetoptError
from glob import glob
from gzip import GzipFile
from libprep3 import *
from e3_defs import *
from e3_grid import *
from e3_block import *
from e3_flow import *

verbosity_level = 0

shortOptions = ""
longOptions = ["help", "job=", "zip-files", "no-zip-files", "vtk-xml", "binary-format", "tecplot", 
               "prepare-restart", "put-into-folders", "tindx=", "output-file=", 
               "slice-list=", "slice-at-point=", "slice-along-line=", "surface-list=",
               "static-flow-profile=", "local-surface-list=",
               "ref-function=", "compare-job=", "compare-tindx=",
               "report-norms", "per-block-norm-list=", "global-norm-list=",
               "probe=", "add-pitot-p", "add-total-p",
               "add-molef", "gmodel-file=",
               "add-user-computed-vars=",
               "add-total-enthalpy", "add-mach", "heat-flux-list=", "vertex-velocity-list=", 
               "plot3d", "omegaz=", "tangent-slab-list=", "prepare-fstc-restart", "moving-grid",
               "add-transport-coeffs", "verbosity=", "bc-surface-list=", "OpenFoam"]

def printUsage():
    print ""
    print "Usage:"
    print "e3post.py [--help] [--job=<jobFileName>] [--tindx=<index|9999|last|all|xxxx>]"
    print "          [--zip-files|--no-zip-files]"
    print "          [--moving-grid]"
    print "          [--omegaz=\"[omegaz0,omegaz1,...]\"]"
    print ""
    print "          [--add-pitot-p] [--add-total-p] [--add-mach] [--add-total-enthalpy]"
    print "          [--add-molef --gmodel-file=\"gas-model.lua\"]"
    print "          [--add-transport-coeffs --gmodel-file=\"gas-model.lua\"]"
    print "          [--add-user-computed-vars=\"user-script.py\"]"
    print ""
    print "          [--vtk-xml] [--binary-format] [--tecplot] [--plot3d] [--OpenFoam]"
    print ""
    print "          [--output-file=<profile-data-file>]"
    print "          [--slice-list=\"blk-range,i-range,j-range,k-range;...\"]"
    print "          [--slice-at-point=\"blk-range,index-pair,x,y,z;...\"]"
    print "          [--slice-along-line=\"x0,y0,z0,x1,y1,z1,N\"]"
    print "          [--surface-list=\"blk,surface-name;...\"]"
    print "          [--local-surface-list=\"blk,surface-name;...\"]"
    print ""
    print "          [--static-flow-profile=\"blk,face-name;...\"]"
    print ""
    print "          [--heat-flux-list=\"blk-range,surf-range,i-range,j-range,k-range;...\"]"
    print "          [--bc-surface-list=\"blk-range,surf-range,i-range,j-range,k-range;...\"]"
    print "          [--tangent-slab-list=\"blk-range,i-range,j-range,k-range;...\"]"
    print ""
    print "          [--probe=\"x,y,z;...\"]"
    print ""
    print "          [--report-norms]"
    print "          [--per-block-norm-list=\"jb,var-name,norm-name;...\""
    print "          [--global-norm-list=\"var-name,norm-name;...\""
    print "          [--ref-function=<python-script>]"
    print "          [--compare-job=<jobFileName> [--compare-tindx=<index>]]"
    print ""
    print "          [--prepare-restart] [--prepare-fstc-restart]"
    print "          [--put-into-folders]"
    print ""
    print "          [--verbosity=<int>]"
    print ""
    print "For further information, see the online documentation, the Eilmer3 User Guide"
    print "and the source code."
    return

#----------------------------------------------------------------------

def read_time_indices():
    """
    The job.times should contain the map from solution index to simulation time.

    The file format is one entry per line, possibly with comment lines 
    starting with the # character.

    Returns the final tindx and a dictionary of (tindx,t)-pairs
    """
    fileName = rootName + ".times"
    try:
        fp = open(fileName, "r")
        # For unkown reason, sometimes a line has several
        # leading null-characters. So remove them to avoid
        # an error while reading the file. (Stefan Hess)
        buf = fp.readline().strip().strip("\0")
        times_dict = {}
        while len(buf) > 0:
            if buf[0] == '#': 
                buf = fp.readline().strip().strip("\0")
                continue
            tokens = buf.split()
            tindx = int(tokens[0]); t = float(tokens[1])
            times_dict[tindx] = t
            buf = fp.readline().strip().strip("\0")
        fp.close()
    except Exception, e:
        # For some reason, we could not read the job.times file
        # so let's assume that there should be a t=0 data set.
        # This will allow us to handle cases where we have run
        # the preparation program but not the main simulation.
        print "Problem reading .times file:", e
        tindx = 0; times_dict = { 0:0.0 }
    return tindx, times_dict


def read_block_labels():
    """
    The block_labels file contains the (string) names of each of the blocks,
    one per line.

    Comment lines start with a # character.
    """
    fileName = "block_labels.list"
    fp = open(fileName, "r")
    buf = fp.readline()
    block_labels = []
    while len(buf) > 0:
        if buf[0] == '#': 
            buf = fp.readline()
            continue
        tokens = buf.split()
        block_labels.append(tokens[1])
        buf = fp.readline()
    return len(block_labels), block_labels


def prepare_for_restart(rootName, tindx):
    """
    Do the house-keeping of changing the .times file and renaming the 9999
    flow data files so that the simulation may be restarted cleanly.

    This function essentially documents the process of setting up for a restart.
    .. Versions: PJ, 19-Feb-2009 for Rainer.
       24-Jan-2010 added some suggestions from Fabs

    This function is always verbose because we want the user to check things.
    """
    print "Preparing a restart from tindx=", tindx
    fileName = rootName + ".times"
    fp = open(fileName, "r")
    content = fp.read()
    fp.close()
    backupFileName = fileName+".original"
    if os.access(backupFileName, os.F_OK):
        print "It seems that you have tried to prepare for a restart before."
        print "   We will not do any work this time.  However, it you really."
        print "   want to proceed, just delete the file", backupFileName
        print "   and try the --prepare-restart command again."
        return
    fp = open(backupFileName, "w")
    fp.write(content)
    fp.close()
    lines = content.split('\n')
    last_tindx = 0
    # Overwrite the times file with the selected indices
    # and (possibly) the last index renamed.
    fp = open(fileName, "w")
    for line in lines:
        if len(line.strip()) == 0:
            continue
        tokens = line.split()
        first_item = tokens[0]
        if first_item == "#":
            fp.write("%s\n" % line)
            continue
        if int(first_item) < tindx:
            fp.write("%s\n" % line)
            last_tindx = int(first_item)
        elif int(first_item) == 9999:
            # We have reached the end of the solution set (tindx=9999)
            # so we need to relabel this solution and the flow files.
            tindx = last_tindx + 1
            tindx_str = "%04d" % tindx
            print "Relabel the 9999 time index to", tindx_str
            t = float(tokens[1])
            dt = float(tokens[2])
            fp.write("%04d %e %e\n" % (tindx, t, dt))
            pattern = os.path.join("flow", "t9999", rootName+".flow.b????.t9999*")
            solution_file_list = glob(pattern)
            for srcFileName in solution_file_list:
                destFileName = srcFileName.replace("9999", tindx_str)
                print "   Rename", srcFileName, "-->", destFileName
                os.renames(srcFileName, destFileName)
            break
        elif int(first_item) == tindx:
            fp.write("%s\n" % line)
            # Truncate the file at this line.
            # Leave the solution files as they are.
            break
    fp.close()
    #
    print "-------------------------------------------------------"
    print "The restart process in semi-automatic and may require"
    print "some manual adjustment of the restart parameters."
    print ""
    print "To confirm, here is the content of your new times file:"
    print ""
    fp = open(fileName, "r"); content = fp.read(); fp.close()
    print content
    print "-------------------------------------------------------"
    controlFileName = rootName + ".control"
    print "You may need to edit the time-stepping control file:", controlFileName
    cp = ConfigParser.ConfigParser()
    cp.read(controlFileName)
    print "Current values in that file are:"
    print "    dt =", cp.get("control_data", "dt")
    print "    max_time =", cp.get("control_data", "max_time")
    print "    max_step =", cp.get("control_data", "max_step")
    finishFileName = rootName + ".finish"
    print "But, in the file", finishFileName, "we see that" 
    cp2 = ConfigParser.ConfigParser()
    cp2.read(finishFileName)
    print "    dt =", cp2.get("simulation_end", "dt")
    print "--------------------------------------------------------"
    print "Done."
    return

def prepare_fstc_restart(rootName, tindx):
    """
    This function documents the process of setting up for a coupled computation restart.

    Do the house-keeping of changing the .times file and renaming the 9999
    flow data files so that the simulation may be restarted cleanly.
    Automatically replace in the .control file the old dt from the .finish file.
    """
    global verbosity_level
    if verbosity_level > 0: print "Preparing a fstc-restart from tindx=", tindx
    fileName = rootName + ".times"
    fp = open(fileName, "r")
    content = fp.read()
    fp.close()
    lines = content.split('\n')
    last_tindx = 0
    # Overwrite the times file with the selected indices
    # and (possibly) the last index renamed.
    fp = open(fileName, "w")
    for line in lines:
        if len(line.strip()) == 0:
            continue
        tokens = line.split()
        first_item = tokens[0]
        if first_item == "#":
            fp.write("%s\n" % line)
            continue
        if int(first_item) < tindx:
            fp.write("%s\n" % line)
            last_tindx = int(first_item)
        elif int(first_item) == 9999:
            # We have reached the end of the solution set (tindx=9999)
            # so we need to relabel this solution and the flow files.
            tindx = last_tindx + 1
            tindx_str = "%04d" % tindx
            if verbosity_level > 0: print "Relabel the 9999 time index to", tindx_str
            t = float(tokens[1])
            dt = float(tokens[2])
            fp.write("%04d %e %e\n" % (tindx, t, dt))
            pattern = os.path.join("flow", "t9999", rootName+".flow.b????.t9999*")
            solution_file_list = glob(pattern)
            for srcFileName in solution_file_list:
                destFileName = srcFileName.replace("9999", tindx_str)
                if verbosity_level > 0: print "   Rename", srcFileName, "-->", destFileName
                os.renames(srcFileName, destFileName)
            break
        elif int(first_item) == tindx:
            fp.write("%s\n" % line)
            # Truncate the file at this line.
            # Leave the solution files as they are.
            break
    fp.close()
    #
    controlFileName = rootName + ".control"
    finishFileName = rootName + ".finish"
    cp = ConfigParser.ConfigParser()
    cp.read(controlFileName)
    cp2 = ConfigParser.ConfigParser()
    cp2.read(finishFileName)
    replace_dt = cp2.get("simulation_end", "dt")
    cp.set("control_data", "dt", replace_dt)
    fp = open(controlFileName,'w')
    cp.write(fp)
    fp.close()
    return

def put_into_folders(rootName, times_dict):
    """
    Put the solution and grid files into folders.
    
    This should be needed only for old simulations where solution files
    were all written into the script-level directory.  Since that made
    the directory very busy and tedious to clean up, the newer arrangement
    writes the grid and solution files to their own directores.
    """
    global verbosity_level
    work_dir = os.getcwd()
    tindx_list = times_dict.keys()
    if verbosity_level > 0:
        print "Put the solution and grid files into folders."
        print "(This should be needed for old simulations, only.)"
        print "Working in directory:", work_dir
        print "time indices:", tindx_list
        print "Move grid files."
    grid_file_list = glob(rootName+".grid.b????.t????*")
    for old_name in grid_file_list:
        new_name = os.path.join("grid", "t0000", old_name)
        if verbosity_level > 0: print "   ", old_name, "-->", new_name
        os.renames(old_name, new_name)
    print "Move flow files."
    for tindx in tindx_list:
        tindx_str = "t%04d" % tindx
        if verbosity_level > 0: print "tindx:", tindx
        flow_file_list = glob(rootName+".flow.b????."+tindx_str+"*")
        for old_name in flow_file_list:
            new_name = os.path.join("flow", tindx_str, old_name)
            if verbosity_level > 0: print "   ", old_name, "-->", new_name
            os.renames(old_name, new_name)
    if verbosity_level > 0: print "Move history files."
    hist_file_list = glob(rootName+".hist.b????*")
    for old_name in hist_file_list:
        new_name = os.path.join("hist", old_name)
        if verbosity_level > 0: print "   ", old_name, "-->", new_name
        os.renames(old_name, new_name)
    if verbosity_level > 0: print "Move plot files."
    plot_file_list = glob(rootName+"*vtu*")
    for old_name in plot_file_list:
        new_name = os.path.join("plot", old_name)
        if verbosity_level > 0: print "   ", old_name, "-->", new_name
        os.renames(old_name, new_name)
    plot_file_list = glob(rootName+"*visit*")
    for old_name in plot_file_list:
        new_name = os.path.join("plot", old_name)
        if verbosity_level > 0: print "   ", old_name, "-->", new_name
        os.renames(old_name, new_name)
    if verbosity_level > 0: print "Done."
    return

#----------------------------------------------------------------------------
# The following function depends on definitions in e3_grid.py, e3_flow.py
# and e3_block.py

def select_surface_from_block(block_grid, block_flow, which_surface):
    """
    Selects one of the bounding surfaces from the block and constructs grid
    and flow objects that have one dimension less.  These objects can then be
    used when writing data to the plotting files so that we can have 
    corresponding VTK files for the surfaces of blocks.

    Note that this operation really only makes sense for 3D blocks.
    """
    which_surface = faceDict[which_surface] # string or integer ---> integer
    if which_surface == BOTTOM or which_surface == TOP:
        niv = block_grid.ni; njv = block_grid.nj
        surface_grid = StructuredGrid((niv,njv,1))
        if which_surface == BOTTOM:
            k = 0
        else:
            k = block_grid.nk-1
        for i in range(niv):
            for j in range(njv):
                surface_grid.x[i,j,0] = block_grid.x[i,j,k]
                surface_grid.y[i,j,0] = block_grid.y[i,j,k]
                surface_grid.z[i,j,0] = block_grid.z[i,j,k]
        nic = niv - 1; njc = njv - 1
        if which_surface == BOTTOM:
            k = 0
        else:
            k = block_flow.nk-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list
        surface_flow.ni = nic; surface_flow.nj = njc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((nic,njc,1),'d')
        for j in range(njc):
            for i in range(nic):
                for var in var_list:
                    surface_flow.data[var][i,j,0] = block_flow.data[var][i,j,k]
    if which_surface == NORTH or which_surface == SOUTH:
        niv = block_grid.ni; nkv = block_grid.nk
        surface_grid = StructuredGrid((niv,nkv,1))
        if which_surface == SOUTH:
            j = 0
        else:
            j = block_grid.nj-1
        for i in range(niv):
            for k in range(nkv):
                surface_grid.x[i,k,0] = block_grid.x[i,j,k]
                surface_grid.y[i,k,0] = block_grid.y[i,j,k]
                surface_grid.z[i,k,0] = block_grid.z[i,j,k]
        nic = niv - 1; nkc = nkv - 1
        if which_surface == SOUTH:
            j = 0
        else:
            j = block_flow.nj-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list;
        surface_flow.ni = nic; surface_flow.nj = nkc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((nic,nkc,1),'d')
        for k in range(nkc):
            for i in range(nic):
                for var in var_list:
                    surface_flow.data[var][i,k,0] = block_flow.data[var][i,j,k]
    if which_surface == WEST or which_surface == EAST:
        njv = block_grid.nj; nkv = block_grid.nk
        surface_grid = StructuredGrid((njv,nkv,1))
        if which_surface == WEST: 
            i = 0
        else:
            i = block_grid.ni-1
        for k in range(nkv):
            for j in range(njv):
                surface_grid.x[j,k,0] = block_grid.x[i,j,k]
                surface_grid.y[j,k,0] = block_grid.y[i,j,k]
                surface_grid.z[j,k,0] = block_grid.z[i,j,k]
        njc = njv - 1; nkc = nkv - 1
        if which_surface == WEST:
            i = 0
        else:
            i = block_flow.ni-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list;
        surface_flow.ni = njc; surface_flow.nj = nkc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((njc,nkc,1),'d')
        for k in range(nkc):
            for j in range(njc):
                for var in var_list:
                    surface_flow.data[var][j,k,0] = block_flow.data[var][i,j,k]
    return surface_grid, surface_flow

# Write the surface data at local frame of reference
#
def select_surface_from_block_local(block_grid, block_flow, which_surface):
    """
    Selects one of the bounding surfaces from the block and constructs grid
    and flow objects that have one dimension less.  These objects can then be
    used when writing data to the plotting files so that we can have 
    corresponding VTK files for the surfaces of blocks at local frame
    of reference.

    Note that this operation really only makes sense for 3D blocks.
    """
    which_surface = faceDict[which_surface] # string or integer ---> integer
    if which_surface == BOTTOM or which_surface == TOP:
        niv = block_grid.ni; njv = block_grid.nj
        surface_grid = StructuredGrid((niv,njv,1))
        if which_surface == BOTTOM:
            k = 0
        else:
            k = block_grid.nk-1
        for i in range(niv):
            for j in range(njv):
                surface_grid.x[i,j,0] = block_grid.x[i,j,k]
                surface_grid.y[i,j,0] = block_grid.y[i,j,k]
                surface_grid.z[i,j,0] = block_grid.z[i,j,k]
        nic = niv - 1; njc = njv - 1
        if which_surface == BOTTOM:
            k = 0
        else:
            k = block_flow.nk-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list
        surface_flow.ni = nic; surface_flow.nj = njc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((nic,njc,1),'d')
        for j in range(njc):
            for i in range(nic):
                for var in var_list:
                    surface_flow.data[var][i,j,0] = block_flow.data[var][i,j,k]
        # Change the velocity vector into local frame of reference
        for j in range(njc):
            for i in range(nic):
                    rotation_angle = math.atan2(block_flow.data["pos.y"][i,j,k], block_flow.data["pos.x"][i,j,k])
                    rotation_angle = math.atan2(block_flow.data["pos.y"][i,j,k], block_flow.data["pos.x"][i,j,k])
                    surface_flow_vel_x_local = block_flow.data["vel.x"][i,j,k]*math.cos(rotation_angle) \
                                               + block_flow.data["vel.y"][i,j,k]*math.sin(rotation_angle)
                    surface_flow_vel_y_local = block_flow.data["vel.y"][i,j,k]*math.cos(rotation_angle) \
                                               - block_flow.data["vel.x"][i,j,k]*math.sin(rotation_angle)
                    surface_flow.data["vel.x"][i,j,0] = surface_flow_vel_x_local
                    surface_flow.data["vel.y"][i,j,0] = surface_flow_vel_y_local
    if which_surface == NORTH or which_surface == SOUTH:
        niv = block_grid.ni; nkv = block_grid.nk
        surface_grid = StructuredGrid((niv,nkv,1))
        if which_surface == SOUTH:
            j = 0
        else:
            j = block_grid.nj-1
        for i in range(niv):
            for k in range(nkv):
                surface_grid.x[i,k,0] = block_grid.x[i,j,k]
                surface_grid.y[i,k,0] = block_grid.y[i,j,k]
                surface_grid.z[i,k,0] = block_grid.z[i,j,k]
        nic = niv - 1; nkc = nkv - 1
        if which_surface == SOUTH:
            j = 0
        else:
            j = block_flow.nj-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list;
        surface_flow.ni = nic; surface_flow.nj = nkc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((nic,nkc,1),'d')
        for k in range(nkc):
            for i in range(nic):
                for var in var_list:
                    surface_flow.data[var][i,k,0] = block_flow.data[var][i,j,k]
        # Change the velocity vector into local frame of reference
        for k in range(nkc):
            for i in range(nic):
                    rotation_angle = math.atan2(block_flow.data["pos.y"][i,j,k], block_flow.data["pos.x"][i,j,k])
                    surface_flow_vel_x_local = block_flow.data["vel.x"][i,j,k]*math.cos(rotation_angle) \
                                               + block_flow.data["vel.y"][i,j,k]*math.sin(rotation_angle)
                    surface_flow_vel_y_local = block_flow.data["vel.y"][i,j,k]*math.cos(rotation_angle) \
                                               - block_flow.data["vel.x"][i,j,k]*math.sin(rotation_angle)
                    surface_flow.data["vel.x"][i,k,0] = surface_flow_vel_x_local
                    surface_flow.data["vel.y"][i,k,0] = surface_flow_vel_y_local
    if which_surface == WEST or which_surface == EAST:
        njv = block_grid.nj; nkv = block_grid.nk
        surface_grid = StructuredGrid((njv,nkv,1))
        if which_surface == WEST: 
            i = 0
        else:
            i = block_grid.ni-1
        for k in range(nkv):
            for j in range(njv):
                surface_grid.x[j,k,0] = block_grid.x[i,j,k]
                surface_grid.y[j,k,0] = block_grid.y[i,j,k]
                surface_grid.z[j,k,0] = block_grid.z[i,j,k]
        njc = njv - 1; nkc = nkv - 1
        if which_surface == WEST:
            i = 0
        else:
            i = block_flow.ni-1
        var_list = copy(block_flow.vars)
        surface_flow = StructuredGridFlow()
        surface_flow.vars = var_list;
        surface_flow.ni = njc; surface_flow.nj = nkc; surface_flow.nk = 1
        for variable in var_list:
            surface_flow.data[variable] = zeros((njc,nkc,1),'d')
        for k in range(nkc):
            for j in range(njc):
                for var in var_list:
                    surface_flow.data[var][j,k,0] = block_flow.data[var][i,j,k]
        # Change the velocity vector into local frame of reference
        for k in range(nkc):
            for j in range(njc):
                    rotation_angle = math.atan2(block_flow.data["pos.y"][i,j,k], block_flow.data["pos.x"][i,j,k])
                    surface_flow_vel_x_local = block_flow.data["vel.x"][i,j,k]*math.cos(rotation_angle) \
                                               + block_flow.data["vel.y"][i,j,k]*math.sin(rotation_angle)
                    surface_flow_vel_y_local = block_flow.data["vel.y"][i,j,k]*math.cos(rotation_angle) \
                                               - block_flow.data["vel.x"][i,j,k]*math.sin(rotation_angle)
                    surface_flow.data["vel.x"][j,k,0] = surface_flow_vel_x_local
                    surface_flow.data["vel.y"][j,k,0] = surface_flow_vel_y_local
    return surface_grid, surface_flow

# --------------------------------------------------------------------
# The following functions and classes are for post-processing the heat-flux data
# into a human readable format
class HeatFluxData(object):
    """
    Python class to store heat-flux data for a single cell interface
    """
    def __init__(self,i=0,j=0,k=0,x=0.0,y=0.0,z=0.0,qc=0.0,qd=0.0,qr=0.0,Twall=0.0,
    	         Tcell=0.0,rho_cell=0.0,un_cell=0.0,Re_wall=0.0):
	"""
	Create a HeatFluxData object from provided data
	"""
	self.i = i
	self.j = j
	self.k = k
	self.x = x
	self.y = y
	self.z = z
	self.qc = qc
	self.qd = qd
	self.qr = qr
	self.Twall = Twall
	self.Tcell = Tcell
	self.rho_cell = rho_cell
	self.un_cell = un_cell
	self.Re_wall = Re_wall

class BoundaryHeatFluxData(object):
    """
    Python class to store heat-flux data for a block surface/boundary
    """
    def __init__(self):
	self.iface = []
	self.irange = [1000000000,-1]
	self.jrange = [1000000000,-1]
	self.krange = [1000000000,-1]
    
    def read(self,fp=None):
	"""
	Read in heat-flux data from file.
	"""
	if fp == None:
	    print "HeatFluxData.read(): no file was provided!"
	    sys.exit()
	buf = fp.readline() # time
	time_stamp = float(buf)
	buf = fp.readline() # variable-name list
	var_list = []
	for token in buf.split():
	    var_list.append(token.strip('"')) # just keep the name
	buf = fp.readline() # surface dimensions
	tks = buf.split()
	ni = int(tks[0]); nj = int(tks[1]); nk = int(tks[2])
	dim = ni*nj*nk
	
	for line in range(dim):
	    buf = fp.readline() # heat-flux data
	    # if len(buf)==0: break
	    tks = buf.split()
	    if len(tks) < 14:
	    	print "BoundaryHeatFluxData::read()"
	    	print "    The heat flux data files may have been created from an older"
	    	print "    version of eilmer3 and therefore Re_wall cannot be extracted."
	    	print "    Continuing on without Re_wall."
	    	Re_wall = 0.0
	    else:
	    	Re_wall = float(tks[13])
	    self.iface.append( HeatFluxData(i=int(tks[0]),j=int(tks[1]),k=int(tks[2]),
	                               x=float(tks[3]),y=float(tks[4]),z=float(tks[5]),
	                               qc=float(tks[6]),qd=float(tks[7]),qr=float(tks[8]),
	                               Twall=float(tks[9]),Tcell=float(tks[10]),
	                               rho_cell=float(tks[11]),un_cell=float(tks[12]),
                                       Re_wall=Re_wall ) )
	    # set ranges
	    # lower
	    if self.iface[-1].i < self.irange[0]: self.irange[0] = self.iface[-1].i
	    if self.iface[-1].j < self.jrange[0]: self.jrange[0] = self.iface[-1].j
	    if self.iface[-1].k < self.krange[0]: self.krange[0] = self.iface[-1].k
	    # upper
	    if self.iface[-1].i > self.irange[1]: self.irange[1] = self.iface[-1].i
	    if self.iface[-1].j > self.jrange[1]: self.jrange[1] = self.iface[-1].j
	    if self.iface[-1].k > self.krange[1]: self.krange[1] = self.iface[-1].k
        return
	    
    def get_iface(self,i,j,k):
	# check if indices are in-range
	if i < self.irange[0] or i > self.irange[1]:
	    print "i = %d is out-of-range"
	    sys.exit()
	if j < self.jrange[0] or j > self.jrange[1]:
	    print "j = %d is out-of-range"
	    sys.exit()
	if k < self.krange[0] or k > self.krange[1]:
	    print "k = %d is out-of-range"
	    sys.exit()
	# search for this data point
	for ihfd in self.iface:
	    if ihfd.i == i and ihfd.j == j and ihfd.k == k:
		return ihfd
	#
	print "BoundaryHeatFluxData.get_iface(): search failed!"
	sys.exit()

def read_all_heat_flux_data(rootName, nblock, tindx, zipFiles=0):
    """
    Returns all heat-flux data for a single flow solution.
    """
    global verbosity_level
    heat = []
    for jb in range(nblock):
	heat.append([])
	for js in range(6):
	    fileName = rootName+".heat"+(".b%04d.s%04d.t%04d" % (jb, js, tindx))
	    fileName = os.path.join("heat", "t%04d" % ( tindx ), fileName)
	    # test if this file exists (required due to disparate 2D/3D boundaries)
	    if os.path.isfile(fileName) == False \
	    	and os.path.isfile(fileName+".gz") == False: 
		# print "file %s does not exist" % ( fileName )
		break
	    if verbosity_level > 0: print "Read heat-flux data from", fileName
	    heat[-1].append(BoundaryHeatFluxData())
	    if zipFiles: 
                fp = GzipFile(fileName+".gz", "rb")
            else:
                fp = open(fileName, "r")
            heat[jb][-1].read(fp)
            fp.close()
    return heat
    
def write_heat_flux_profile(outputFileName, heat_flux_list_str, tindx, nblock, hf_data ):
    """
    Extracts and writes to file a profile of heat-flux data from a collection
    of surface/boundary slices.
    """
    global verbosity_level
    fp = open(outputFileName, "w")
    # write header
    fp.write("# Filename: %s\n" % outputFileName)
    fp.write("# Column 1: Distance along surface\n")
    fp.write("# Column 2: Conductive heat flux, q_cond (W/m**2)\n")
    fp.write("# Column 3: Diffusive heat flux, q_diff (W/m**2)\n")
    fp.write("# Column 4: Radiative heat flux, q_rad (W/m**2)\n")
    fp.write("# Column 5: Wall temperature, T_wall (K)\n")
    fp.write("# Column 6: Cell temperature, T_cell (K)\n")
    fp.write("# Column 7: Cell density, rho_cell (kg/m**3)\n")
    fp.write("# Column 8: Cell normal velocity, un_cell (m/s)\n")
    fp.write("# Column 9: Wall Reynolds number, Re_wall (ND)\n")
    fp.write("# Column 10: pos.x (m)\n")
    fp.write("# Column 11: pos.y (m)\n")
    fp.write("# Column 12: pos.z (m)\n")
    heat_flux_lists = heat_flux_list_str.split(';')
    if verbosity_level > 0: print "heat_flux_lists = ", heat_flux_lists
    first = True
    L = 0.0
    for heat_flux_str in heat_flux_lists:
	bstr,sstr,istr,jstr,kstr = heat_flux_str.split(',')
	bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)
        print bfirst, blast
	for jb in range(bfirst,blast+1):
	    sfirst,slast = decode_range_from_string(sstr, 0, len(hf_data[jb])-1)
	    for js in range(sfirst,slast+1):
		kfirst,klast = decode_range_from_string(kstr, hf_data[jb][js].krange[0], 
		    					hf_data[jb][js].krange[1])
		jfirst,jlast = decode_range_from_string(jstr, hf_data[jb][js].jrange[0], 
		    					hf_data[jb][js].jrange[1])
		ifirst,ilast = decode_range_from_string(istr, hf_data[jb][js].irange[0], 
		    					hf_data[jb][js].irange[1])
		if verbosity_level > 0:
                    print ("slice jb=%d js=%d i=%d:%d, j=%d:%d, k=%d:%d" %
                           (jb,js,ifirst,ilast,jfirst,jlast,kfirst,klast))
		for k in range(kfirst,klast+1):
		    for j in range(jfirst,jlast+1):
			for i in range(ifirst,ilast+1):
			    iface_data = hf_data[jb][js].get_iface(i,j,k)
			    pos = Vector3(iface_data.x,iface_data.y,iface_data.z)
			    if first:
				pos_prev = Vector3(iface_data.x,iface_data.y,
				    			iface_data.z)
				first = False
			    L += vabs(pos-pos_prev)
			    pos_prev = pos
			    fp.write("%e %e %e %e %e %e %e %e %e %e %e %e\n" % \
                                         ( L, iface_data.qc, iface_data.qd, iface_data.qr,
                                           iface_data.Twall, iface_data.Tcell,
                                           iface_data.rho_cell, iface_data.un_cell, 
                                           iface_data.Re_wall, 
                                           iface_data.x, iface_data.y, iface_data.z) )
    #
    return 0
    
# --------------------------------------------------------------------
# The following functions and classes are for post-processing the vertex velocity data
# into a human readable format
class VertexVelocityData(object):
    """
    Python class to store velocity data for a single cell vertex
    """
    def __init__(self,i=0,j=0,k=0,x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0,v=0.0):
	"""
	Create a VertexVelocityData object from provided data
	"""
	self.i = i
	self.j = j
	self.k = k
	self.x = x
	self.y = y
	self.z = z
	self.vx = vx
	self.vy = vy
	self.vz = vz
	self.v = v	


class BoundaryVertexVelocityData(object):
    """
    Python class to store vertex velocity data for a block surface/boundary
    """
    def __init__(self):
	self.vtx = []
	self.irange = [1000000000,-1]
	self.jrange = [1000000000,-1]
	self.krange = [1000000000,-1]
    
    def read(self,fp=None):
	"""
	Read in vertex velocity data from file.
	"""
	if fp == None:
	    print "BoundaryVertexVelocityData.read(): no file was provided!"
	    sys.exit()
	buf = fp.readline() # time
	time_stamp = float(buf)
	buf = fp.readline() # variable-name list
	var_list = []
	for token in buf.split():
	    var_list.append(token.strip('"')) # just keep the name
	buf = fp.readline() # surface dimensions
	tks = buf.split()
	ni = int(tks[0]); nj = int(tks[1]); nk = int(tks[2])
	dim = ni*nj*nk
	#
	for line in range(dim):
	    buf = fp.readline() # vertex velocity data
	    #if len(buf)==0: break
	    tks = buf.split()
	    v = ( float(tks[6])**2 + float(tks[7])**2 + float(tks[8])**2 )**0.5
	    self.vtx.append( VertexVelocityData(i=int(tks[0]),j=int(tks[1]),k=int(tks[2]),
	                               x=float(tks[3]),y=float(tks[4]),z=float(tks[5]),
	                               vx=float(tks[6]),vy=float(tks[7]),vz=float(tks[8]),
	                               v=v) )
	    # set ranges
	    # lower
	    if self.vtx[-1].i < self.irange[0]: self.irange[0] = self.vtx[-1].i
	    if self.vtx[-1].j < self.jrange[0]: self.jrange[0] = self.vtx[-1].j
	    if self.vtx[-1].k < self.krange[0]: self.krange[0] = self.vtx[-1].k
	    # upper
	    if self.vtx[-1].i > self.irange[1]: self.irange[1] = self.vtx[-1].i
	    if self.vtx[-1].j > self.jrange[1]: self.jrange[1] = self.vtx[-1].j
	    if self.vtx[-1].k > self.krange[1]: self.krange[1] = self.vtx[-1].k
        return
	    
    def get_vtx(self,i,j,k):
	# check if indices are in-range
	if i < self.irange[0] or i > self.irange[1]:
	    print "i = %d is out-of-range"
	    sys.exit()
	if j < self.jrange[0] or j > self.jrange[1]:
	    print "j = %d is out-of-range"
	    sys.exit()
	if k < self.krange[0] or k > self.krange[1]:
	    print "k = %d is out-of-range"
	    sys.exit()
	# search for this data point
	for ihfd in self.vtx:
	    if ihfd.i == i and ihfd.j == j and ihfd.k == k:
		return ihfd
	#
	print "BoundaryVertexVelocityData.get_vtx(): search failed!"
	sys.exit()
  
def read_all_vertex_velocity_data(rootName, nblock, tindx, zipFiles=0):
    """
    Returns all heat-flux data for a single flow solution.
    """
    global verbosity_level
    velocity = []
    for jb in range(nblock):
	velocity.append([])
	for js in range(6):
	    fileName = rootName+".vel"+(".b%04d.s%04d.t%04d" % (jb, js, tindx))
	    fileName = os.path.join("vel", "t%04d" % ( tindx ), fileName)
	    # test if this file exists (required due to disparate 2D/3D boundaries)
	    if os.path.isfile(fileName) == False \
	    	and os.path.isfile(fileName+".gz") == False: 
		# print "file %s does not exist" % ( fileName )
		break
            if verbosity_level > 0:
                print "Read vertex velocity data from", fileName
	    velocity[-1].append(BoundaryVertexVelocityData())
	    if zipFiles: 
                fp = GzipFile(fileName+".gz", "rb")
            else:
                fp = open(fileName, "r")
            velocity[jb][-1].read(fp)
            fp.close()
    return velocity

def write_vertex_velocity_profile(outputFileName, vertex_velocity_list_str, tindx, nblock, vtxv_data ):
    """
    Extracts and writes to file a profile of vertex velocity data from a collection
    of surface/boundary slices.
    """
    global verbosity_level
    fp = open(outputFileName, "w")
    # write header
    fp.write("# Filename: %s\n" % outputFileName)
    fp.write("# Column 1: Distance along surface\n")
    fp.write("# Column 2: Vertex velocity, magnitude (m/s)\n")    
    fp.write("# Column 3: Vertex velocity, x-direction (m/s)\n")
    fp.write("# Column 4: Vertex velocity, y-direction (m/s)\n")
    fp.write("# Column 5: Vertex velocity, z-direction (m/s)\n")
    fp.write("# Column 6: Vertex position, x (m)\n")
    fp.write("# Column 7: Vertex position, y (m)\n")
    fp.write("# Column 8: Vertex position, z (m)\n")

    vertex_velocity_lists = vertex_velocity_list_str.split(';')
    if verbosity_level > 0:
        print "vertex_velocity_lists = ", vertex_velocity_lists
    first = True
    L = 0.0
    for vertex_velocity_str in vertex_velocity_lists:
	bstr,sstr,istr,jstr,kstr = vertex_velocity_str.split(',')
	bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)

	for jb in range(bfirst,blast+1):
	    sfirst,slast = decode_range_from_string(sstr, 0, len(vtxv_data[jb])-1)
	    for js in range(sfirst,slast+1):
		kfirst,klast = decode_range_from_string(kstr, vtxv_data[jb][js].krange[0], 
		    					vtxv_data[jb][js].krange[1])
		jfirst,jlast = decode_range_from_string(jstr, vtxv_data[jb][js].jrange[0], 
		    					vtxv_data[jb][js].jrange[1])
		ifirst,ilast = decode_range_from_string(istr, vtxv_data[jb][js].irange[0], 
		    					vtxv_data[jb][js].irange[1])
		print ("slice jb=%d js=%d i=%d:%d, j=%d:%d, k=%d:%d" %
                   (jb,js,ifirst,ilast,jfirst,jlast,kfirst,klast))
		for k in range(kfirst,klast+1):
		    for j in range(jfirst,jlast+1):
			for i in range(ifirst,ilast+1):
			    vtx_data = vtxv_data[jb][js].get_vtx(i,j,k)
			    pos = Vector3(vtx_data.x,vtx_data.y,vtx_data.z)
			    if first:
				pos_prev = Vector3(vtx_data.x,vtx_data.y,
				    			vtx_data.z)
				first = False
			    L += vabs(pos-pos_prev)
			    pos_prev = pos
			    fp.write("%e %e %e %e %e %e %e %e \n" % ( L,
						vtx_data.v,
						vtx_data.vx,
						vtx_data.vy,
						vtx_data.vz,
						vtx_data.x,
						vtx_data.y,
						vtx_data.z ) )
    
    return 0

# --------------------------------------------------------------------
# The following functions and classes are for post-processing the surface data
# into a human readable format
class SurfaceData(object):
    """
    Python class to store surface data for a single cell interface
    """
    def __init__(self,i=0,j=0,k=0,x=0.0,y=0.0,z=0.0,T_surface=0.0,u_x=0.0,u_y=0.0,u_z=0.0,
    	         tke_surface=0.0,omega_surface=0.0,mass_flux=0.0):
	"""
	Create a SurfaceData object from provided data
	"""
	self.i = i
	self.j = j
	self.k = k
	self.x = x
	self.y = y
	self.z = z
	self.T_surface = T_surface
	self.u_x = u_x
	self.u_y = u_y
	self.u_z = u_z
	self.tke_surface = tke_surface
	self.omega_surface = omega_surface
	self.mass_flux = mass_flux

class BoundarySurfaceData(object):
    """
    Python class to store surface data for a block surface/boundary
    """
    def __init__(self):
	self.iface = []
	self.irange = [1000000000,-1]
	self.jrange = [1000000000,-1]
	self.krange = [1000000000,-1]
    
    def read(self,fp=None):
	"""
	Read in surface data from file.
	"""
	if fp == None:
	    print "SurfaceData.read(): no file was provided!"
	    sys.exit()
	buf = fp.readline() # time
	time_stamp = float(buf)
	buf = fp.readline() # variable-name list
	var_list = []
	for token in buf.split():
	    var_list.append(token.strip('"')) # just keep the name
	buf = fp.readline() # surface dimensions
	tks = buf.split()
	ni = int(tks[0]); nj = int(tks[1]); nk = int(tks[2])
	dim = ni*nj*nk
	
	for line in range(dim):
	    buf = fp.readline() 
	    # if len(buf)==0: break
	    tks = buf.split()
	    self.iface.append( SurfaceData(i=int(tks[0]),j=int(tks[1]),k=int(tks[2]),
	                               x=float(tks[3]),y=float(tks[4]),z=float(tks[5]),
	                               T_surface=float(tks[6]),u_x=float(tks[7]),u_y=float(tks[8]),
	                               u_z=float(tks[9]),tke_surface=float(tks[10]),
	                               omega_surface=float(tks[11]),mass_flux=float(tks[12]) ) )
	    # set ranges
	    # lower
	    if self.iface[-1].i < self.irange[0]: self.irange[0] = self.iface[-1].i
	    if self.iface[-1].j < self.jrange[0]: self.jrange[0] = self.iface[-1].j
	    if self.iface[-1].k < self.krange[0]: self.krange[0] = self.iface[-1].k
	    # upper
	    if self.iface[-1].i > self.irange[1]: self.irange[1] = self.iface[-1].i
	    if self.iface[-1].j > self.jrange[1]: self.jrange[1] = self.iface[-1].j
	    if self.iface[-1].k > self.krange[1]: self.krange[1] = self.iface[-1].k
        return
	    
    def get_iface(self,i,j,k):
	# check if indices are in-range
	if i < self.irange[0] or i > self.irange[1]:
	    print "i = %d is out-of-range"
	    sys.exit()
	if j < self.jrange[0] or j > self.jrange[1]:
	    print "j = %d is out-of-range"
	    sys.exit()
	if k < self.krange[0] or k > self.krange[1]:
	    print "k = %d is out-of-range"
	    sys.exit()
	# search for this data point
	for ihfd in self.iface:
	    if ihfd.i == i and ihfd.j == j and ihfd.k == k:
		return ihfd
	#
	print "BoundarySurfaceData.get_iface(): search failed!"
	sys.exit()

def read_all_surface_data(rootName, nblock, tindx, zipFiles=0):
    """
    Returns all surface data for a single flow solution.
    """
    global verbosity_level
    surf = []
    for jb in range(nblock):
	surf.append([])
	for js in range(6):
	    fileName = rootName+".surf"+(".b%04d.s%04d.t%04d" % (jb, js, tindx))
	    fileName = os.path.join("surf", "t%04d" % ( tindx ), fileName)
	    # test if this file exists (required due to disparate 2D/3D boundaries)
	    if os.path.isfile(fileName) == False \
	    	and os.path.isfile(fileName+".gz") == False: 
		# print "file %s does not exist" % ( fileName )
		break
	    if verbosity_level > 0: print "Read surface data from", fileName
	    surf[-1].append(BoundarySurfaceData())
	    if zipFiles: 
                fp = GzipFile(fileName+".gz", "rb")
            else:
                fp = open(fileName, "r")
            surf[jb][-1].read(fp)
            fp.close()
    return surf
    
def write_surface_profile(outputFileName, bc_surface_list_str, tindx, nblock, sf_data ):
    """
    Extracts and writes to file a profile of surface data from a collection
    of surface/boundary slices.
    """
    global verbosity_level
    fp = open(outputFileName, "w")
    # write header
    fp.write("# Filename: %s\n" % outputFileName)
    fp.write("# Column 1: Surface temperature, T_surface (K)\n")
    fp.write("# Column 2: x-axis velocity, u_x (m/s)\n")
    fp.write("# Column 3: y-axis velocity, u_y (m/s)\n")
    fp.write("# Column 4: z-axis velocity, u_z (m/s)\n")
    fp.write("# Column 5: surface turbulence kinetic energy, tke_surface (m**2/s**2)\n")
    fp.write("# Column 6: surface dissipation rate, omega_surface (1/s)\n")
    fp.write("# Column 7: mass flux, mass_flux (kg/m**3/s)\n")
    fp.write("# Column 8: pos.x (m)\n")
    fp.write("# Column 9: pos.y (m)\n")
    fp.write("# Column 10: pos.z (m)\n")
    surface_lists = bc_surface_list_str.split(';')
    if verbosity_level > 0: print "surface_lists = ", surface_lists
    first = True
    L = 0.0
    for surface_str in surface_lists:
	bstr,sstr,istr,jstr,kstr = surface_str.split(',')
	bfirst,blast = decode_range_from_string(bstr, 0, nblock-1)
        print bfirst, blast
	for jb in range(bfirst,blast+1):
	    sfirst,slast = decode_range_from_string(sstr, 0, len(sf_data[jb])-1)
	    for js in range(sfirst,slast+1):
		kfirst,klast = decode_range_from_string(kstr, sf_data[jb][js].krange[0], 
		    					sf_data[jb][js].krange[1])
		jfirst,jlast = decode_range_from_string(jstr, sf_data[jb][js].jrange[0], 
		    					sf_data[jb][js].jrange[1])
		ifirst,ilast = decode_range_from_string(istr, sf_data[jb][js].irange[0], 
		    					sf_data[jb][js].irange[1])
		if verbosity_level > 0:
                    print ("slice jb=%d js=%d i=%d:%d, j=%d:%d, k=%d:%d" %
                           (jb,js,ifirst,ilast,jfirst,jlast,kfirst,klast))
		for k in range(kfirst,klast+1):
		    for j in range(jfirst,jlast+1):
			for i in range(ifirst,ilast+1):
			    iface_data = sf_data[jb][js].get_iface(i,j,k)
			    pos = Vector3(iface_data.x,iface_data.y,iface_data.z)
			    if first:
				pos_prev = Vector3(iface_data.x,iface_data.y,
				    			iface_data.z)
				first = False
			    L += vabs(pos-pos_prev)
			    pos_prev = pos
			    fp.write("%e %e %e %e %e %e %e %e %e %e\n" % \
                                         ( iface_data.T_surface, iface_data.u_x, iface_data.u_y,
                                           iface_data.u_z, iface_data.tke_surface, iface_data.omega_surface,
                                           iface_data.mass_flux,
                                           iface_data.x, iface_data.y, iface_data.z) )
    #
    return 0
    
def flatten(L):
    """
    Flatten a list of items, which may or may not be lists.

    This function lifted from
    http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    """
    out = []
    for item in L:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out

def parse_user_script(fname):
    execfile(fname, globals())
    # Test that the user supplied a list of variable names
    try:
        if not isinstance(var_names, list):
            print "In user script: ", fname
            print "the variable names should be supplied as a list called 'var_names'"
            print "Cannot continue without this list."
            print "Bailing out!"
            sys.exit(1)
    except NameError:
        print "In user script: ", fname
        print "It appears that the list of variable names was not set."
        print "A list called 'var_names' should be supplied."
        print "Cannot continue without this list."
        print "Bailing out!"
        sys.exit(1)
    # And test that the user supplied a function to call
    try:
        if not callable(compute_vars):
            print "In user script: ", fname
            print "the variable 'compute_vars' is not defined as a function."
            print "Cannot continue without this function defined."
            print "Bailing out!"
            sys.exit(1)
    except NameError:
        print "In user script: ", fname
        print "the 'compute_vars' function is not defined."
        print "Cannont continue without this function defined."
        print "Bailing out!"
        sys.exit(1)

    return var_names, compute_vars

if __name__ == '__main__':
    try:
        userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    except GetoptError, e:
        print "One (or more) of your command-line options was no good."
        print "    ", e
        printUsage()
        sys.exit(1)
    uoDict = dict(userOptions[0])
    verbosity_level = int(uoDict.get("--verbosity", "0"))
    if verbosity_level > 0 or len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        print "Begin e3post.py..."
        print "Source code revision string: ", get_revision_string()
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
        sys.exit(0)
    #
    jobName = uoDict.get("--job", "test")
    rootName, ext = os.path.splitext(jobName)
    zipFiles = True  # Default: use zip file format for grid and flow data files.
    if uoDict.has_key("--no-zip-files"): zipFiles = False
    if uoDict.has_key("--zip-files"): zipFiles = True
    movingGrid = uoDict.has_key("--moving-grid")
    if uoDict.has_key("--omegaz"):
        mystr = uoDict.get("--omegaz", None)
        print "mystr=", mystr
        # We may have to flatten the list if it has been constructed
        # as a list containing other lists.
        # e.g. "[0.0,2100.0,3*[0.0,]]"
        omegaz = flatten(eval(mystr))
        print "omegaz=", omegaz, "len(omegaz)=", len(omegaz)
    else:
        omegaz = None
    #
    final_time_indx, times_dict = read_time_indices()
    nblock, block_labels = read_block_labels()
    tindx_str = uoDict.get("--tindx", str(final_time_indx))
    #
    if uoDict.has_key("--put-into-folders"):
        put_into_folders(rootName, times_dict)
        sys.exit(0)
    #
    # Restarts have to have a particular tindx integer value 
    # from which to start the continuing simulation.
    if uoDict.has_key("--prepare-restart"):
        tindx = int(tindx_str)
        prepare_for_restart(rootName, tindx)
        sys.exit(0)
    #
    if uoDict.has_key("--prepare-fstc-restart"):
        tindx = int(tindx_str)
        prepare_fstc_restart(rootName, tindx)
        sys.exit(0)    
    #
    # If we get to here, we are not preparing for a restart,
    # so it may make sense to be dealing with multiple values for tindx.
    if tindx_str == "all":
        tindx_list = times_dict.keys()
        tindx_list.sort()
    elif tindx_str in ["last", "9999"]:
        tindx_list = times_dict.keys()
        tindx_list.sort()
        tindx_list = [tindx_list[-1],]
    elif tindx_str in ["xxxx", "XXXX"]:
        tindx_list = ["xxxx",]
        # Now, we need to dip into the txxxx solution file,
        # extract the simulation time value from the first line and
        # add it to times_dict.
        times_dict["xxxx"] = read_time_from_flow_file(rootName, "xxxx", zipFiles)
    else:
        try:
            tindx_list = [int(tindx_str,10),]
        except:
            raise ValueError("Do not know what to do with option tindx= %s" % tindx_str)
    #
    if uoDict.has_key("--vtk-xml") or uoDict.has_key("--ref-function") or \
            uoDict.has_key("--compare-job") or uoDict.has_key("--surface-list") or \
            uoDict.has_key("--local-surface-list"): 
        # At this point, the tindx_list may have several entries and 
        # the Visit and PVD files accumulate information about each entry.
        # This allows the construction of animations built from multiple
        # solution frames.
        begin_Visit_file(rootName, nblock)
        begin_PVD_file(rootName)

    aux_var_names = None
    compute_vars = None
    if uoDict.has_key("--add-user-computed-vars"):
        uname = uoDict.get("--add-user-computed-vars", "dummy-script-name")
        aux_var_names, compute_vars = parse_user_script(uname)
    #
    # For the times that have been specified, do something...
    if verbosity_level > 0: print "About to process the following solutions:", tindx_list
    for tindx in tindx_list:
        if uoDict.has_key("--vtk-xml"):
            if verbosity_level > 0:
                print "Assemble VTK-XML files for t=", times_dict[tindx]
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles, movingGrid)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            write_VTK_XML_files(rootName, tindx, nblock, grid, flow, times_dict[tindx],
                                uoDict.has_key("--binary-format"))
        #
        if uoDict.has_key("--OpenFoam"):
            if verbosity_level > 0:
                print "writing OpenFoam grid"
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles, movingGrid)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            write_OpenFoam_files(rootName, nblock, grid, flow)
        #
        if uoDict.has_key("--tecplot"):
            if verbosity_level > 0:
                print "Assemble Tecplot file for t=", times_dict[tindx]
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles, movingGrid)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            write_Tecplot_file(rootName, tindx, nblock, grid, flow, times_dict[tindx])
        #
        if uoDict.has_key("--plot3d"):
            if verbosity_level > 0:
                print "Write out Plot3d grid for t=", times_dict[tindx]
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles, movingGrid)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            # tindx may be an integer, or already a string such as "xxxx"
            if type(tindx) is int:
                tindx_str = "%04d" % tindx
            elif type(tindx) is string:
                tindx_str = tindx
            else:
                raise RuntimeError("WTF: tindx is neither an int nor string.")
            fname = rootName+(".t%04s" % tindx_str)+".grd"
            plotPath = "plot"
            if not os.access(plotPath, os.F_OK):
                os.makedirs(plotPath)
            fname = os.path.join(plotPath, fname)
            if verbosity_level > 0:
                print "Writing true (node-centred) grid for t=", times_dict[tindx]
            write_plot3d_grid(fname, grid)
            if verbosity_level > 0:
                print "Writing cell-centred grid and function files for t=", times_dict[tindx]
            write_plot3d_files(rootName, tindx, nblock, grid, flow, times_dict[tindx])
        #
        if uoDict.has_key("--slice-list") or uoDict.has_key("--slice-at-point"):
            outputFileName = uoDict.get("--output-file", "profile.data")
            if verbosity_level > 0:
                print "Extract slices of data for t=", times_dict[tindx]
                print "    outputFileName=", outputFileName
            slice_list_str = uoDict.get("--slice-list", "")
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            if uoDict.has_key("--slice-at-point"):
                slice_list_str += convert_string(uoDict.get("--slice-at-point", ""),
                                                 nblock, grid, flow)
            if len(slice_list_str) > 0:
                if slice_list_str[0] == ';':
                    # Drop the leading semicolon if it is present.
                    # It may have been added by convert_string()
                    slice_list_str = slice_list_str[1:]
                write_profile_data(outputFileName, slice_list_str, tindx, nblock, grid, flow)
            else:
                print "Probably and error; no slices written."
        #
        if uoDict.has_key("--slice-along-line"):
            outputFileName = uoDict.get("--output-file", "profile.data")
            if verbosity_level > 0:
                print "Extract slices of data for t=", times_dict[tindx]
                print "    outputFileName=", outputFileName
            slice_line_str = uoDict.get("--slice-along-line", "")
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            if len(slice_line_str) > 0:
                write_profile_along_line(outputFileName, slice_line_str, tindx, nblock, 
                                         grid, flow, dimensions)
            else:
                print "Probably and error; no slices written."
        #
        if uoDict.has_key("--ref-function"):
            aScriptName = uoDict.get("--ref-function", "")
            if len(aScriptName) > 0:
                if verbosity_level > 0:
                    print "Compare with reference function for t=", times_dict[tindx]
                    print "   ref_function is from script:", aScriptName
                execfile(aScriptName)
                grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
                add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
                compute_difference_in_flow_data(ref_function, nblock, grid, flow, 
                                                times_dict[tindx])
                write_VTK_XML_files(rootName, tindx, nblock, grid, flow, times_dict[tindx],
                                    uoDict.has_key("--binary-format"))
                norms = compute_volume_weighted_norms(nblock, grid, flow)
                pretty_print_norms(norms,
                                   uoDict.get("--per-block-norm-list", ""),
                                   uoDict.get("--global-norm-list", ""))
        #
        if uoDict.has_key("--compare-job"):
            compareJobName = uoDict.get("--compare-job", jobName)
            compareRootName, compareExt = os.path.splitext(compareJobName)
            compareTindx = uoDict.get("--compare-tindx", "9999")
            if compareTindx in ["last", "9999"]:
                compareTindx = tindx_list[-1]
            else:
                compareTindx = int(compareTindx)
            if verbosity_level > 0:
                print "Compare with reference solution for t=", times_dict[tindx]
                print "   Comparison solution:", compareRootName, " compareTindx=", compareTindx
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            grid2, flow2, dimensions = read_all_blocks(compareRootName, nblock, compareTindx, zipFiles)
            add_auxiliary_variables(nblock, flow2, uoDict, omegaz, aux_var_names, compute_vars)
            compute_difference_in_flow_data2(nblock, grid, flow, grid2, flow2, times_dict[tindx])
            write_VTK_XML_files(rootName, tindx, nblock, grid, flow, times_dict[tindx],
                                uoDict.has_key("--binary-format"))
            norms = compute_volume_weighted_norms(nblock, grid, flow)
            pretty_print_norms(norms,
                               uoDict.get("--per-block-norm-list", ""),
                               uoDict.get("--global-norm-list", ""))
        #
        if uoDict.has_key("--report-norms"):
            if verbosity_level > 0:
                print "Report norms for t=", times_dict[tindx]
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            norms = compute_volume_weighted_norms(nblock, grid, flow)
            pretty_print_norms(norms, 
                               uoDict.get("--per-block-norm-list", ""),
                               uoDict.get("--global-norm-list", ""))
        #
        if uoDict.has_key("--surface-list"):
            if verbosity_level > 0:
                print "Extract a set of surfaces for t=", times_dict[tindx], "and write as VTK files."
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            surface_list_str = uoDict.get("--surface-list", "")
            surface_list_str = surface_list_str.lower().strip()
            if verbosity_level > 0:
                print "surface_list_str=", surface_list_str
            surface_list = surface_list_str.split(";")
            surface_grid_list = []
            surface_flow_list = []
            for surf in surface_list:
                items = surf.split(",")
                blk_id = int(items[0])
                which_surface = faceDict[items[1]]
                sgrid, sflow = select_surface_from_block(grid[blk_id], flow[blk_id], which_surface)
                surface_grid_list.append(sgrid)
                surface_flow_list.append(sflow)
            outputName = uoDict.get("--output-file", rootName)
            write_VTK_XML_files(outputName+".surface", tindx, 
                                len(surface_grid_list), 
                                surface_grid_list, surface_flow_list, times_dict[tindx],
                                uoDict.has_key("--binary-format"))
        #
        if uoDict.has_key("--local-surface-list"):
            if verbosity_level > 0:
                print "Extract a set of surfaces for t=", times_dict[tindx], "at local frame of reference and write as VTK files."
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            surface_list_str = uoDict.get("--local-surface-list", "")
            surface_list_str = surface_list_str.lower().strip()
            if verbosity_level > 0:
                print "local_surface_list_str=", surface_list_str
            surface_list = surface_list_str.split(";")
            surface_grid_list = []
            surface_flow_list = []
            for surf in surface_list:
                items = surf.split(",")
                blk_id = int(items[0])
                which_surface = faceDict[items[1]]
                sgrid, sflow = select_surface_from_block_local(grid[blk_id], flow[blk_id], which_surface)
                surface_grid_list.append(sgrid)
                surface_flow_list.append(sflow)
            outputName = uoDict.get("--output-file", rootName)
            write_VTK_XML_files(outputName+".local-surface", tindx, 
                                len(surface_grid_list), 
                                surface_grid_list, surface_flow_list, times_dict[tindx],
                                uoDict.has_key("--binary-format"))
        #
        if uoDict.has_key("--static-flow-profile"):
            if verbosity_level > 0:
                print "Write a set of static-flow profiles for t=", times_dict[tindx]
                print "(to be used for simulation inflow)"
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            blkface_list_str = uoDict.get("--static-flow-profile", "")
            blkface_list_str = blkface_list_str.lower().strip()
            if verbosity_level > 0:
                print "blkface_list_str=", blkface_list_str
            blkface_list = blkface_list_str.split(";")
            blkface_grid_list = []
            blkface_flow_list = []
            for surf in blkface_list:
                items = surf.split(",")
                blk_id = int(items[0])
                which_surface = faceDict[items[1]] # converts string or number to int
                outputName = rootName + "-blk-" + str(blk_id) + "-face-" \
                    + faceName[which_surface] + ".static-flow-profile"
                write_static_flow_profile_from_block(outputName, flow[blk_id], which_surface)
        #
        if uoDict.has_key("--probe"):
            if verbosity_level > 0:
                print "Probe data for t=", times_dict[tindx]
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            add_auxiliary_variables(nblock, flow, uoDict, omegaz, aux_var_names, compute_vars)
            # Pull apart coordinates list
            coordinate_list_str = uoDict.get("--probe", "")
            for xyz_str in coordinate_list_str.split(';'):
                coords = xyz_str.split(',')
                x = float(coords[0]); y = 0.0; z = 0.0
                if len(coords) > 1: y = float(coords[1])
                if len(coords) > 2: z = float(coords[2])
                jb, i, j, k = locate_cell_and_block(grid, flow, dimensions, 0, 0, 0, 0, x, y, z)
                if verbosity_level > 0:
                    print "    coords=(%g,%g,%g)" % (x, y, z)
                    print "    jb=", jb, "ijk=", i, j, k
                flow[jb].write_gnuplot_header(sys.stdout)
                flow[jb].write_gnuplot_data_for_cell(sys.stdout, i, j, k)
        #
        if uoDict.has_key("--heat-flux-list"):
            outputFileName = uoDict.get("--output-file", "hf_profile.data")
            if verbosity_level > 0:
                print "Extract heat flux for t=", times_dict[tindx], "and write as text files."
                print "    outputFileName=", outputFileName
            hf_data = read_all_heat_flux_data(rootName, nblock, tindx, zipFiles)
            heat_flux_list_str = uoDict.get("--heat-flux-list", "")
            if len(heat_flux_list_str) > 0:
		write_heat_flux_profile(outputFileName, heat_flux_list_str, tindx, nblock, hf_data )
            else:
                print "Probably an error; no heat flux profile written."
        #
        if uoDict.has_key("--vertex-velocity-list"):
            outputFileName = uoDict.get("--output-file", "vtxv_profile.data")
            if verbosity_level > 0:
                print "Extract vertex velocities for t=", times_dict[tindx], "and write as text files."
                print "    outputFileName=", outputFileName
            vtxv_data = read_all_vertex_velocity_data(rootName, nblock, tindx, zipFiles)
            vertex_velocity_list_str = uoDict.get("--vertex-velocity-list", "")
            if len(vertex_velocity_list_str) > 0:
		write_vertex_velocity_profile(outputFileName, vertex_velocity_list_str, tindx, nblock, vtxv_data )
            else:
                print "Probably and error; no vertex velocity list written."
        #
        if uoDict.has_key("--tangent-slab-list"):
            outputFileName = uoDict.get("--output-file", "TS-profile.data")
            if verbosity_level > 0:
                print "Perform tangent-slab calculation on slices of data for t=", times_dict[tindx]
                print "    outputFileName=", outputFileName
            slice_list_str = uoDict.get("--tangent-slab-list", "")
            grid, flow, dimensions = read_all_blocks(rootName, nblock, tindx, zipFiles)
            if len(slice_list_str) > 0:
                if slice_list_str[0] == ';':
                    # Drop the leading semicolon if it is present.
                    # It may have been added by convert_string()
                    slice_list_str = slice_list_str[1:]
                tangent_slab_along_slice(outputFileName, slice_list_str, tindx, nblock, grid, flow)
            else:
                print "Probably an error; no tangent-slab calculations performed."

        #
        if uoDict.has_key("--bc-surface-list"):
            outputFileName = uoDict.get("--output-file", "sf_profile.data")
            if verbosity_level > 0:
                print "Extract surface data for t=", times_dict[tindx], "and write as text files."
                print "    outputFileName=", outputFileName
            sf_data = read_all_surface_data(rootName, nblock, tindx, zipFiles)
            bc_surface_list_str = uoDict.get("--bc-surface-list", "")
            if len(bc_surface_list_str) > 0:
		write_surface_profile(outputFileName, bc_surface_list_str, tindx, nblock, sf_data )
            else:
                print "Probably an error; no surface profile written."
    #
    if uoDict.has_key("--vtk-xml") or uoDict.has_key("--ref-function") or \
            uoDict.has_key("--compare-job") or uoDict.has_key("--surface-list") or \
            uoDict.has_key("--local-surface-list"): 
        finish_PVD_file(rootName)

    if verbosity_level > 0: print "End e3post.py."
    sys.exit(0)
