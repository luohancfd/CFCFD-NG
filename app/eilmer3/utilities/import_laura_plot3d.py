#!/usr/bin/env python
# Author: Rowan J. Gollan
# Date: 26-Apr-2010, 02-Sep-2010 integrated into Eilmer3 (PJ).
# Place: NASA Langley, Hampton, Virginia, USA
#
# This script will convert a Laura solution (in 
# Plot3D fortran unformatted form) to a VTK file for 
# viewing in Paraview.  
# It can handle solutions with multiple blocks and 
# should be able to handle solutions with different 
# numbers of of flow variables included.
#

from math import sqrt
from numpy import zeros, shape
import re
import sys
import os
from getopt import getopt
from gzip import GzipFile
sys.path.append("/sw/lib/python2.3/site-packages/Numeric") # for Tim's MacOSX
sys.path.append(os.path.expandvars("$HOME/e3bin")) # installation directory
from cfpylib.util.FortranFile import FortranFile
from libprep3 import *
from e3_grid import StructuredGrid, read_plot3d_grid_binary_file
from e3_flow import uflowz, quoted_string

shortOptions = ""
longOptions = ["help", "job="]

def printUsage():
    print ""
    print "Usage: import_laura_plot3d.py [--help] [--job=<jobName>]"
    return

#------------------------------------------------------------------
# Functions to read Laura's plot3d name, grid and data files.
# The mesh is block-structured and the flow data is vertex-based.

def create_var_list(nfile):
    nf = open(nfile, "r")
    var_names = []
    while 1:
        line = nf.readline()
        if not line:
            break
        var_names.append(line[1:-1])
    for i in range(len(var_names)):
        # Strip out xml-style structures and internal spaces
        var_names[i] = re.sub('<sub>', '_', var_names[i])
        var_names[i] = re.sub('</sub>', '', var_names[i])
        var_names[i] = re.sub('<greek>', '', var_names[i])
        var_names[i] = re.sub('</greek>', '', var_names[i])
        var_names[i] = re.sub(' ', '', var_names[i])
    print "found names=", var_names
    return var_names

def read_laura_function_file(fname, var_list):
    """
    Reads the multiblock plot3d function file and
    returns lists of dimensions and flow data at cell vertices.

    When indexing particular data items, data[m][n][i,j,k]
    is the data item for block m, variable n, vertex [i,j,k].
    """
    f = FortranFile(fname)
    nblocks = f.readInts()[0]
    ndim = f.readInts()
    ni = []; nj = []; nk = []; nvar = [];
    for m in range(nblocks):
        ni.append(ndim[4*m])
        nj.append(ndim[4*m+1])
        nk.append(ndim[4*m+2])
        nvar.append(ndim[4*m+3])
    #
    d = f.readReals('d')
    data = []
    pos = 0
    for m in range(nblocks):
        data.append({})
        assert nvar[m] == len(var_list)
        for name in var_list:
            data[m][name] = zeros((ni[m], nj[m], nk[m]))
            for k in range(nk[m]):
                for j in range(nj[m]):
                    for i in range(ni[m]):
                        data[m][name][i,j,k] = d[pos]
                        pos = pos + 1
    #
    return (ni, nj, nk, data)

def node_data_to_cell_data(ni, nj, nk, grids, data, var_list):
    # Do on a block-by-block basis
    cells = []
    for m in range(len(grids)):
        cells.append({})
        nic = ni[m] - 1
        njc = nj[m] - 1
        nkc = nk[m] - 1
        for n in var_list:
            cells[m][n] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.x'] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.y'] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.z'] = zeros((nic, njc, nkc), 'd')
        cells[m]['vol'] = zeros((nic, njc, nkc), 'd')
        #
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    p0 = Vector(grids[m].x[i,j,k],       grids[m].y[i,j,k],       grids[m].z[i,j,k])
                    p1 = Vector(grids[m].x[i+1,j,k],     grids[m].y[i+1,j,k],     grids[m].z[i+1,j,k])
                    p2 = Vector(grids[m].x[i+1,j+1,k],   grids[m].y[i+1,j+1,k],   grids[m].z[i+1,j+1,k])
                    p3 = Vector(grids[m].x[i,j+1,k],     grids[m].y[i,j+1,k],     grids[m].z[i,j+1,k])
                    p4 = Vector(grids[m].x[i,j,k+1],     grids[m].y[i,j,k+1],     grids[m].z[i,j,k+1])
                    p5 = Vector(grids[m].x[i+1,j,k+1],   grids[m].y[i+1,j,k+1],   grids[m].z[i+1,j,k+1])
                    p6 = Vector(grids[m].x[i+1,j+1,k+1], grids[m].y[i+1,j+1,k+1], grids[m].z[i+1,j+1,k+1])
                    p7 = Vector(grids[m].x[i,j+1,k+1],   grids[m].y[i,j+1,k+1],   grids[m].z[i,j+1,k+1])
                    centre = hexahedron_centroid(p0, p1, p2, p3, p4, p5, p6, p7)
                    vol = hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7)
                    #
                    try:
                        a = quad_centroid(p0, p1, p2, p3)
                        b = quad_centroid(p1, p5, p6, p2)
                        c = quad_centroid(p4, p5, p6, p7)
                        d = quad_centroid(p0, p4, p7, p3)
                        e = quad_centroid(p0, p1, p5, p4)
                        f = quad_centroid(p3, p2, p6, p7)
                        #
                        p01 = 0.5*(p0 + p1)
                        p12 = 0.5*(p1 + p2)
                        p23 = 0.5*(p2 + p3)
                        p30 = 0.5*(p3 + p0)
                        p15 = 0.5*(p1 + p5)
                        p56 = 0.5*(p5 + p6)
                        p62 = 0.5*(p6 + p2)
                        p67 = 0.5*(p6 + p7)
                        p74 = 0.5*(p7 + p4)
                        p45 = 0.5*(p4 + p5)
                        p73 = 0.5*(p7 + p3)
                        p40 = 0.5*(p4 + p0)
                        #
                        w0 = hexahedron_volume(centre, b, p62, f, c, p56, p6, p67) / vol
                        w1 = hexahedron_volume(d, centre, f, p73, p74, c, p67, p7) / vol
                        w2 = hexahedron_volume(p40, e, centre, d, p4, p45, c, p74) / vol
                        w3 = hexahedron_volume(e, p15, b, centre, p45, p5, p56, c) / vol
                        w4 = hexahedron_volume(a, p12, p2, p23, centre, b, p62, f) / vol
                        w5 = hexahedron_volume(p30, a, p23, p3, d, centre, f, p73) / vol
                        w6 = hexahedron_volume(p0, p01, a, p30, p40, e, centre, d) / vol
                        w7 = hexahedron_volume(p01, p1, p12, a, e, p15, b, centre) / vol
                    except ValueError:
                        # One side of our hexahedron has most likely collapsed,
                        # so just take a regular average (not volume-weighted)
                        w0 = 0.125; w1 = 0.125; w2 = 0.125; w3 = 0.125
                        w4 = 0.125; w5 = 0.125; w6 = 0.125; w7 = 0.125
                    #
                    for n in var_list:
                        cells[m][n][i,j,k] = w0*data[m][n][i,j,k] + \
                                             w1*data[m][n][i+1,j,k] + \
                                             w2*data[m][n][i+1,j+1,k] + \
                                             w3*data[m][n][i,j+1,k] + \
                                             w4*data[m][n][i,j,k+1] + \
                                             w5*data[m][n][i+1,j,k+1] + \
                                             w6*data[m][n][i+1,j+1,k+1] + \
                                             w7*data[m][n][i,j+1,k+1]
                    #
                    cells[m]['pos.x'][i,j,k] = centre.x
                    cells[m]['pos.y'][i,j,k] = centre.y
                    cells[m]['pos.z'][i,j,k] = centre.z
                    cells[m]['vol'][i,j,k] = vol
    #
    var_list_cells = ['pos.x', 'pos.y', 'pos.z', 'vol'] + var_list
    return cells, var_list_cells

#------------------------------------------------------------------------
# Functions for writing VTK plot files.

def write_VTK_p_unstructured_file(baseFileName, ni, nj, nk, grids, cells, var_list):
    """
    Write the coordinating file for the partitioned data then
    one partitioned file per block-structured mesh.
    """
    fname = baseFileName + '.pvtu'
    fp = open(fname, 'w')
    fp.write("<VTKFile type=\"PUnstructuredGrid\">\n")
    fp.write("<PUnstructuredGrid GhostLevel=\"0\">")
    fp.write("<PPoints>\n")
    fp.write(" <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    fp.write("</PPoints>\n")
    fp.write("<PCellData>\n")
    for v in var_list:
        fp.write(" <PDataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n" % v)
    fp.write(" <PDataArray Name=\"vol\" type=\"Float32\" NumberOfComponents=\"1\"/>\n")
    if ( 'vel.x' in var_list ) and ( 'vel.y' in var_list ) and ( 'vel.z' in var_list ):
        fp.write(" <PDataArray Name=\"vel\" type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    fp.write("</PCellData>\n")
    for m in range(len(grids)):
        f2name = baseFileName + (".b%02d.vtu" % m)
        fp.write("<Piece Source=\"%s\"/>\n" % f2name)
        write_VTK_unstructured_file(f2name, ni[m], nj[m], nk[m], grids[m], cells[m], var_list)
    fp.write("</PUnstructuredGrid>\n")
    fp.write("</VTKFile>\n")
    fp.close()
    return


def write_VTK_unstructured_file(fname, ni, nj, nk, grid, cell, var_list):
    """
    One of these files per block.
    """
    fp = open(fname, 'w')
    nic = ni - 1 # number of cells in i-direction
    njc = nj - 1
    nkc = nk - 1
    NumberOfPoints = ni*nj*nk
    NumberOfCells = nic*njc*nkc
    fp.write("<VTKFile type=\"UnstructuredGrid\">\n")
    fp.write("<UnstructuredGrid>")
    fp.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" % (NumberOfPoints, NumberOfCells))
    fp.write("<Points>\n")
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\">\n")
    vtx_number = 0
    vtx_id = {}
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                fp.write(" %e %e %e\n" % (grid.x[i,j,k], grid.y[i,j,k], grid.z[i,j,k]))
                vtx_id[(i,j,k)] = vtx_number
                vtx_number += 1
    #
    fp.write(" </DataArray>\n")
    fp.write("</Points>\n")
    fp.write("<Cells>\n")
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\">\n")
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                fp.write(" %d %d %d %d %d %d %d %d\n" % 
                         (vtx_id[(i,j,k)], vtx_id[(i+1,j,k)], 
                          vtx_id[(i+1,j+1,k)], vtx_id[(i,j+1,k)],
                          vtx_id[(i,j,k+1)], vtx_id[(i+1,j,k+1)], 
                          vtx_id[(i+1,j+1,k+1)], vtx_id[(i,j+1,k+1)]))
    fp.write(" </DataArray>\n")
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\">\n")
    # Since all of the point-lists are concatenated, these offsets into the connectivity
    # array specify the end of each cell.
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                fp.write(" %d\n" % (8*(1+i+j*nic+k*(nic*njc))))
    fp.write(" </DataArray>\n")
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\">\n")
    for k in range(nkc):
        for j in range(njc):
            for i in range(nic):
                fp.write(" %d\n" % 12) # VTK_HEXAHEDRON
    fp.write(" </DataArray>\n")
    fp.write("</Cells>\n")
    #
    fp.write("<CellData>\n")
    for name in var_list:
        fp.write(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\">\n" % name)
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    fp.write(" %e\n" % uflowz(cell[name][i,j,k]))
        fp.write(" </DataArray>\n")
    # Assemble velocity vector if available
    if ( 'vel.x' in var_list ) and ( 'vel.y' in var_list ) and ( 'vel.z' in var_list ):
        fp.write(" <DataArray Name=\"vel\" type=\"Float32\" NumberOfComponents=\"3\">\n")
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    fp.write(" %e %e %e\n" % (uflowz(cell['vel.x'][i,j,k]),
                                              uflowz(cell['vel.y'][i,j,k]),
                                              uflowz(cell['vel.z'][i,j,k])))
        fp.write(" </DataArray>\n")
    #
    fp.write("</CellData>\n")
    fp.write("</Piece>\n")
    fp.write("</UnstructuredGrid>\n")
    fp.write("</VTKFile>\n")
    fp.close()
    return

#-----------------------------------------------------------------------
# Functions for writing (native) Eilmer3 grid and flow files.

def write_Eilmer3_flow_file(fp, nni, nnj, nnk, cells, var_list):
    """
    Write a single-block data in Eilmer3 native format.
    """
    fp.write("%20.12e\n" % 0.0)
    fp.write("%s\n" % quoted_string(var_list))
    fp.write("%d %d %d\n" % (nni, nnj, nnk))
    for k in range(nnk):
        for j in range(nnj):
            for i in range(nni):
                for n, name in enumerate(var_list):
                    if n > 0: fp.write(" ")
                    fp.write("%20.12e" % cells[name][i,j,k])
                fp.write("\n")
    return

def write_Eilmer3_files(rootName, nig, njg, nkg, grids, cells, var_list):
    """
    Write the grid and data files as if for Eilmer3.

    These can be picked up by e3post.py
    """
    fp = open(rootName+".times", "w")
    fp.write("# tindx sim_time dt_global\n")
    fp.write("%04d %e %e\n" % (0, 0.0, 1.0e-6))
    fp.close()
    #
    fp = open("block_labels.list", "w")
    fp.write("# indx label\n");
    for jb in range(len(grids)):
        fp.write("%d block-%d\n" % (jb,jb) )
    fp.close()
    #
    # Write one file per block into the grid/t0000/ directory.
    gridPath = os.path.join("grid", "t0000")
    if not os.access(gridPath, os.F_OK):
        os.makedirs(gridPath)
    for jb in range(len(grids)):
        fileName = rootName+(".grid.b%04d.t0000" % jb)
        fileName = os.path.join(gridPath, fileName)
        fp = GzipFile(fileName+".gz", "wb")
        grids[jb].write(fp)
        fp.close()
    #
    # Again, write one flow file per block into the flow/t0000/ directory.
    flowPath = os.path.join("flow", "t0000")
    if not os.access(flowPath, os.F_OK):
        os.makedirs(flowPath)
    for jb in range(len(cells)):
        fileName = rootName+(".flow.b%04d.t0000" % jb)
        fileName = os.path.join(flowPath, fileName)
        fp = GzipFile(fileName+".gz", "wb")
        write_Eilmer3_flow_file(fp, nig[jb]-1, njg[jb]-1, nkg[jb]-1, cells[jb], var_list)
        fp.close()
    return

#-----------------------------------------------------------------------

def main(uoDict):
    """
    Do some real work.
    """
    jobName = uoDict.get("--job", 'laura')
    plot3dGridFile = jobName + '.g'
    plot3dFunctionFile = jobName + '.q'
    plot3dNamesFile = jobName + '.nam'

    var_list_plot3d = create_var_list(plot3dNamesFile)
    # Change from Laura velocity names to Eilmer names.
    var_list_plot3d[var_list_plot3d.index('u')] = 'vel.x'
    var_list_plot3d[var_list_plot3d.index('v')] = 'vel.y'
    var_list_plot3d[var_list_plot3d.index('w')] = 'vel.z'
    print "Reading in Laura grid..."
    nig, njg, nkg, grids = read_plot3d_grid_binary_file(plot3dGridFile)
    print "Done."
    print "Reading in Laura solution..."
    nif, njf, nkf, data = read_laura_function_file(plot3dFunctionFile, var_list_plot3d)
    print "Done."
    assert nig == nif
    assert njg == njf
    assert nkg == nkf
    print "Converting node data to cell data..."
    cells, var_list_cells = node_data_to_cell_data(nig, njg, nkg, grids, data, var_list_plot3d)
    print "Done."
    print "Writing out VTK file..."
    write_VTK_p_unstructured_file(jobName, nig, njg, nkg, grids, cells, var_list_cells)
    print "Done."
    print "Writing Eilmer3 files..."
    write_Eilmer3_files(jobName, nig, njg, nkg, grids, cells, var_list_cells)
    print "Conversion complete."
    return


if __name__ == '__main__':
    print "Begin import_laura_plot3d.py..."
    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or uoDict.has_key("--help"):
        printUsage()
    else:
        main(uoDict)
    print "Done."
    sys.exit(0)

