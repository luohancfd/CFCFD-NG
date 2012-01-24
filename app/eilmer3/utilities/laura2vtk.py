#!/usr/bin/env python
# Author: Rowan J. Gollan
# Date: 26-Apr-2010
# Place: NASA Langley, Hampton, Virginia, USA
#
# This stand-alone script will convert a Laura solution 
# (in Plot3D Fortran unformatted form) to a VTK file for 
# viewing in Paraview.  It can handle solutions
# with multiple blocks and should be able to handle
# solutions with different numbers of of flow variables
# included.
#

from math import sqrt
from numpy import zeros, zeros_like, shape
import re

#
# Global hard-coded variables that may need changing
#

gname = 'laura.g'     # grid file
fname = 'laura.q'     # function file
nfile = 'laura.nam'   # name file
obfn = 'laura'        # base file name for output

# END: hard-coded variables

#-----------------------------------------------------
# Inclusion of a FortranFile class found at:
# http://www.scipy.org/Cookbook/FortranIO/FortranFile
# on 26-Apr-2010
#
# Copyright 2008, 2009 Neil Martinsen-Burrell
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""Defines a file-derived class to read/write Fortran unformatted files.

The assumption is that a Fortran unformatted file is being written by
the Fortran runtime as a sequence of records.  Each record consists of
an integer (of the default size [usually 32 or 64 bits]) giving the
length of the following data in bytes, then the data itself, then the
same integer as before.

Examples
--------

To use the default endian and size settings, one can just do::
    >>> f = FortranFile('filename')
    >>> x = f.readReals()

One can read arrays with varying precisions::
    >>> f = FortranFile('filename')
    >>> x = f.readInts('h')
    >>> y = f.readInts('q')
    >>> z = f.readReals('f')
Where the format codes are those used by Python's struct module.

One can change the default endian-ness and header precision::
    >>> f = FortranFile('filename', endian='>', header_prec='l')
for a file with little-endian data whose record headers are long
integers.
"""

__docformat__ = "restructuredtext en"

import struct
import numpy

class FortranFile(file):

    """File with methods for dealing with fortran unformatted data files"""

    def _get_header_length(self):
        return struct.calcsize(self._header_prec)
    _header_length = property(fget=_get_header_length)

    def _set_endian(self,c):
        """Set endian to big (c='>') or little (c='<') or native (c='@')

        :Parameters:
          `c` : string
            The endian-ness to use when reading from this file.
        """
        if c in '<>@=':
            self._endian = c
        else:
            raise ValueError('Cannot set endian-ness')
    def _get_endian(self):
        return self._endian
    ENDIAN = property(fset=_set_endian,
                      fget=_get_endian,
                      doc="Possible endian values are '<', '>', '@', '='"
                     )

    def _set_header_prec(self, prec):
        if prec in 'hilq':
            self._header_prec = prec
        else:
            raise ValueError('Cannot set header precision')
    def _get_header_prec(self):
        return self._header_prec
    HEADER_PREC = property(fset=_set_header_prec,
                           fget=_get_header_prec,
                           doc="Possible header precisions are 'h', 'i', 'l', 'q'"
                          )

    def __init__(self, fname, endian='@', header_prec='i', *args, **kwargs):
        """Open a Fortran unformatted file for writing.
        
        Parameters
        ----------
        endian : character, optional
            Specify the endian-ness of the file.  Possible values are
            '>', '<', '@' and '='.  See the documentation of Python's
            struct module for their meanings.  The deafult is '>' (native
            byte order)
        header_prec : character, optional
            Specify the precision used for the record headers.  Possible
            values are 'h', 'i', 'l' and 'q' with their meanings from
            Python's struct module.  The default is 'i' (the system's
            default integer).

        """
        file.__init__(self, fname, *args, **kwargs)
        self.ENDIAN = endian
        self.HEADER_PREC = header_prec

    def _read_exactly(self, num_bytes):
        """Read in exactly num_bytes, raising an error if it can't be done."""
        data = ''
        while True:
            l = len(data)
            if l == num_bytes:
                return data
            else:
                read_data = self.read(num_bytes - l)
            if read_data == '':
                raise IOError('Could not read enough data.'
                              '  Wanted %d bytes, got %d.' % (num_bytes, l))
            data += read_data

    def _read_check(self):
        return struct.unpack(self.ENDIAN+self.HEADER_PREC,
                             self._read_exactly(self._header_length)
                            )[0]

    def _write_check(self, number_of_bytes):
        """Write the header for the given number of bytes"""
        self.write(struct.pack(self.ENDIAN+self.HEADER_PREC,
                               number_of_bytes))

    def readRecord(self):
        """Read a single fortran record"""
        l = self._read_check()
        data_str = self._read_exactly(l)
        check_size = self._read_check()
        if check_size != l:
            raise IOError('Error reading record from data file')
        return data_str

    def writeRecord(self,s):
        """Write a record with the given bytes.

        Parameters
        ----------
        s : the string to write

        """
        length_bytes = len(s)
        self._write_check(length_bytes)
        self.write(s)
        self._write_check(length_bytes)

    def readString(self):
        """Read a string."""
        return self.readRecord()

    def writeString(self,s):
        """Write a string

        Parameters
        ----------
        s : the string to write
        
        """
        self.writeRecord(s)

    _real_precisions = 'df'

    def readReals(self, prec='f'):
        """Read in an array of real numbers.
        
        Parameters
        ----------
        prec : character, optional
            Specify the precision of the array using character codes from
            Python's struct module.  Possible values are 'd' and 'f'.
            
        """
        
        _numpy_precisions = {'d': numpy.float64,
                             'f': numpy.float32
                            }

        if prec not in self._real_precisions:
            raise ValueError('Not an appropriate precision')
            
        data_str = self.readRecord()
        num = len(data_str)/struct.calcsize(prec)
        numbers =struct.unpack(self.ENDIAN+str(num)+prec,data_str) 
        return numpy.array(numbers, dtype=_numpy_precisions[prec])

    def writeReals(self, reals, prec='f'):
        """Write an array of floats in given precision

        Parameters
        ----------
        reals : array
            Data to write
        prec` : string
            Character code for the precision to use in writing
        """
        if prec not in self._real_precisions:
            raise ValueError('Not an appropriate precision')
        
        # Don't use writeRecord to avoid having to form a
        # string as large as the array of numbers
        length_bytes = len(reals)*struct.calcsize(prec)
        self._write_check(length_bytes)
        _fmt = self.ENDIAN + prec
        for r in reals:
            self.write(struct.pack(_fmt,r))
        self._write_check(length_bytes)
    
    _int_precisions = 'hilq'

    def readInts(self, prec='i'):
        """Read an array of integers.
        
        Parameters
        ----------
        prec : character, optional
            Specify the precision of the data to be read using 
            character codes from Python's struct module.  Possible
            values are 'h', 'i', 'l' and 'q'
            
        """
        if prec not in self._int_precisions:
            raise ValueError('Not an appropriate precision')
            
        data_str = self.readRecord()
        num = len(data_str)/struct.calcsize(prec)
        return numpy.array(struct.unpack(self.ENDIAN+str(num)+prec,data_str))

    def writeInts(self, ints, prec='i'):
        """Write an array of integers in given precision

        Parameters
        ----------
        reals : array
            Data to write
        prec : string
            Character code for the precision to use in writing
        """
        if prec not in self._int_precisions:
            raise ValueError('Not an appropriate precision')
        
        # Don't use writeRecord to avoid having to form a
        # string as large as the array of numbers
        length_bytes = len(ints)*struct.calcsize(prec)
        self._write_check(length_bytes)
        _fmt = self.ENDIAN + prec
        for item in ints:
            self.write(struct.pack(_fmt,item))
        self._write_check(length_bytes)

#----- end FortranFile class --------------------


#-----------------------------------------
# A bare minimum geometry library to do
# some of the work required below.
#-----------------------------------------

VERY_SMALL_MAGNITUDE = 1.0e-200

class Vector(object):
    __slots__ = 'x', 'y', 'z'
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z
        return

    def __str__(self):
        return "Vector(%g,%g,%g)" % (self.x, self.y, self.z)

    def __repr__(self):
        return self.__str__()

    def __add__(self, other):
        return Vector(self.x+other.x, self.y+other.y, self.z+other.z)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)

    def __pos__(self):
        return Vector(self.x, self.y, self.z)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, Vector):
            return Vector(self.x*other.x, self.y*other.y, self.z*other.z)
        else:
            return Vector(self.x*other, self.y*other, self.z*other)
        
    def __rmul__(self, other):
        return (self * other)

    def __div__(self, other):
        return Vector(self.x/other, self.y/other, self.z/other)

    def sum(self):
        return self.x + self.y + self.z

    def __abs__(self):
        return sqrt((self * self).sum())

    def unit(self):
        mag = abs(self)
        if mag <= 0.0:
            raise ValueError, "Zero magnitude vector has no defined direction."
        return Vector(self.x/mag, self.y/mag, self.z/mag)
    

def dot(a, b):
    return (a * b).sum()

def cross(a, b):
    x = a.y * b.z - a.z * b.y
    y = a.z * b.x - a.x * b.z
    z = a.x * b.y - a.y * b.x
    return Vector(x, y, z)

def quad_properties(p0, p1, p2, p3):
    vector_area = 0.5 * (cross(p0-p3, p2-p3) + cross(p1-p0, p2-p1))
    n = vector_area.unit()
    area = abs(vector_area)
    t1 = (p1 - p0).unit()
    t2 = cross(n, t1)
    centroid = 0.25 * (p0 + p1 + p2 + p3)
    return centroid, n, t1, t2, area

def quad_centroid(p0, p1, p2, p3):
    centroid, n, t1, t2, area = quad_properties(p0, p1, p2, p3)
    return centroid

def tetrahedron_properties(p0, p1, p2, p3):
    volume = dot(p3-p0, cross(p1-p0, p2-p0))/6.0
    centroid = 0.25 * (p0 + p1 + p2 + p3)
    return centroid, volume

def wedge_properties(p0, p1, p2, p3, p4, p5):
    c1, v1 = tetrahedron_properties(p0, p4, p5, p3)
    c2, v2 = tetrahedron_properties(p0, p5, p4, p1)
    c3, v3 = tetrahedron_properties(p0, p1, p2, p5)
    volume = v1 + v2 + v3
    if volume < VERY_SMALL_MAGNITUDE:
        #print "Warning wedge_properties():"
        #print "Very small or negative volume: ", volume
        #print "Setting volume to zero."
        volume = 0.0
        centroid = (c1 + c2 + c3)/3.0
        return centroid, volume

    centroid = (c1*v1 + c2*v2 + c3*v3)/volume
    return centroid, volume

def hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7):
    c1, v1 = wedge_properties(p0, p1, p2, p4, p5, p6)
    c2, v2 = wedge_properties(p0, p2, p3, p4, p6, p7)
    volume = v1 + v2
    if volume < VERY_SMALL_MAGNITUDE:
        #print "Warning hexahedron_properties():"
        #print "Very small or negative volume: ", volume
        #print "Setting volume to zero."
        volume = 0.0
        centroid = 0.5*(c1 + c2)
        return centroid, volume

    centroid = (c1*v1 + c2*v2)/volume
    return centroid, volume

def hexahedron_volume(p0, p1, p2, p3, p4, p5, p6, p7):
    c, v = hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7)
    return v

#------- end geometry library ----------------------------

def create_var_list(nfile):
    nf = open(nfile, "r")
    var_names = []
    while 1:
        line = nf.readline()
        if not line:
            break
        var_names.append(line[1:-1])

    for i in range(len(var_names)):
        # Strip out xml-style structures
        var_names[i] = re.sub('<sub>', '_', var_names[i])
        var_names[i] = re.sub('</sub>', '', var_names[i])
        var_names[i] = re.sub('<greek>', '', var_names[i])
        var_names[i] = re.sub('</greek>', '', var_names[i])

    return var_names

def read_laura_grid_file(gname):
    g = FortranFile(gname)
    nblocks = g.readInts()[0]
    ndim = g.readInts()
    ni = []; nj = []; nk = []
    for m in range(nblocks):
        ni.append(ndim[3*m])
        nj.append(ndim[3*m+1])
        nk.append(ndim[3*m+2])
    
    pts = g.readReals('d')
    grids = []
    pos = 0
    for m in range(nblocks):
        grids.append({})
        grids[m]['x'] = zeros((ni[m], nj[m], nk[m]))
        grids[m]['y'] = zeros_like(grids[m]['x'])
        grids[m]['z'] = zeros_like(grids[m]['x'])
        # Read x values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m]['x'][i,j,k] = pts[pos]
                    pos = pos + 1
        # Read y values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m]['y'][i,j,k] = pts[pos]
                    pos = pos + 1
        # Read z values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m]['z'][i,j,k] = pts[pos]
                    pos = pos + 1
    
    return (ni, nj, nk, grids)

def read_laura_function_file(fname):
    f = FortranFile(fname)
    nblocks = f.readInts()[0]
    ndim = f.readInts()
    ni = []; nj = []; nk = []; nvar = [];
    for m in range(nblocks):
        ni.append(ndim[4*m])
        nj.append(ndim[4*m+1])
        nk.append(ndim[4*m+2])
        nvar.append(ndim[4*m+3])
    
    d = f.readReals('d')
    data = []
    pos = 0
    for m in range(nblocks):
        data.append([])
        for n in range(nvar[m]):
            data[m].append([])
            data[m][n] = zeros((ni[m], nj[m], nk[m]))
            for k in range(nk[m]):
                for j in range(nj[m]):
                    for i in range(ni[m]):
                        data[m][n][i,j,k] = d[pos]
                        pos = pos + 1
    
    return (ni, nj, nk, data)

def node_data_to_cell_data(ni, nj, nk, grids, data):
    # Do on a block-by-block basis
    cells = []
    for m in range(len(grids)):
        cells.append({})
        nic = ni[m] - 1
        njc = nj[m] - 1
        nkc = nk[m] - 1
        for n in range(len(data[m])):
            cells[m][n] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.x'] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.y'] = zeros((nic, njc, nkc), 'd')
        cells[m]['pos.z'] = zeros((nic, njc, nkc), 'd')
        cells[m]['vol'] = zeros((nic, njc, nkc), 'd')

        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    p0 = Vector(grids[m]['x'][i,j,k],       grids[m]['y'][i,j,k],       grids[m]['z'][i,j,k])
                    p1 = Vector(grids[m]['x'][i+1,j,k],     grids[m]['y'][i+1,j,k],     grids[m]['z'][i+1,j,k])
                    p2 = Vector(grids[m]['x'][i+1,j+1,k],   grids[m]['y'][i+1,j+1,k],   grids[m]['z'][i+1,j+1,k])
                    p3 = Vector(grids[m]['x'][i,j+1,k],     grids[m]['y'][i,j+1,k],     grids[m]['z'][i,j+1,k])
                    p4 = Vector(grids[m]['x'][i,j,k+1],     grids[m]['y'][i,j,k+1],     grids[m]['z'][i,j,k+1])
                    p5 = Vector(grids[m]['x'][i+1,j,k+1],   grids[m]['y'][i+1,j,k+1],   grids[m]['z'][i+1,j,k+1])
                    p6 = Vector(grids[m]['x'][i+1,j+1,k+1], grids[m]['y'][i+1,j+1,k+1], grids[m]['z'][i+1,j+1,k+1])
                    p7 = Vector(grids[m]['x'][i,j+1,k+1],   grids[m]['y'][i,j+1,k+1],   grids[m]['z'][i,j+1,k+1])
                    centre, vol = hexahedron_properties(p0, p1, p2, p3, p4, p5, p6, p7)

                    try:
                        a = quad_centroid(p0, p1, p2, p3)
                        b = quad_centroid(p1, p5, p6, p2)
                        c = quad_centroid(p4, p5, p6, p7)
                        d = quad_centroid(p0, p4, p7, p3)
                        e = quad_centroid(p0, p1, p5, p4)
                        f = quad_centroid(p3, p2, p6, p7)

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

                    for n in range(len(data[m])):
                        cells[m][n][i,j,k] = w0*data[m][n][i,j,k] + \
                                             w1*data[m][n][i+1,j,k] + \
                                             w2*data[m][n][i+1,j+1,k] + \
                                             w3*data[m][n][i,j+1,k] + \
                                             w4*data[m][n][i,j,k+1] + \
                                             w5*data[m][n][i+1,j,k+1] + \
                                             w6*data[m][n][i+1,j+1,k+1] + \
                                             w7*data[m][n][i,j+1,k+1]

                            
                    cells[m]['pos.x'][i,j,k] = centre.x
                    cells[m]['pos.y'][i,j,k] = centre.y
                    cells[m]['pos.z'][i,j,k] = centre.z
                    cells[m]['vol'][i,j,k] = vol

    return cells

def uflowz(q):
    """
    Set very small quantities to zero, exactly.

    This is intended primarily to avoid the bad behaviour of VTK
    when reading Float32 values that are *too* small.
    """
    if abs(q) > 1.0e-30:
        return q
    else:
        return 0.0

def write_VTK_p_unstructured_file(bfn, ni, nj, nk, grids, cells, var_list):
    fname = bfn + '.pvtu'
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
    if ( 'u' in var_list ) and ( 'v' in var_list ) and ( 'w' in var_list ):
        fp.write(" <PDataArray Name=\"vel\" type=\"Float32\" NumberOfComponents=\"3\"/>\n")
    fp.write("</PCellData>\n")
    for m in range(len(grids)):
        f2name = bfn + (".b%02d.vtu" % m)
        fp.write("<Piece Source=\"%s\"/>\n" % f2name)
        write_VTK_unstructured_file(f2name, ni[m], nj[m], nk[m], grids[m], cells[m], var_list)
    fp.write("</PUnstructuredGrid>\n")
    fp.write("</VTKFile>\n")
    fp.close()
    return


def write_VTK_unstructured_file(fname, ni, nj, nk, grid, cell, var_list):
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
                fp.write(" %e %e %e\n" % (grid['x'][i,j,k],
                                          grid['y'][i,j,k],
                                          grid['z'][i,j,k]))
                vtx_id[(i,j,k)] = vtx_number
                vtx_number += 1

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

    fp.write("<CellData>\n")
    for n, v in enumerate(var_list):
        fp.write(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\">\n" % v)
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    fp.write(" %e\n" % uflowz(cell[n][i,j,k]))
        fp.write(" </DataArray>\n")
    fp.write(" <DataArray Name=\"vol\" type=\"Float32\" NumberOfComponents=\"1\">\n")
    for k in range(nkc):
        for j in range(njc):
                for i in range(nic):
                    fp.write(" %e\n" % uflowz(cell['vol'][i,j,k]))
    fp.write(" </DataArray>\n")

    # Assemble velocity vector if available
    if ( 'u' in var_list ) and ( 'v' in var_list ) and ( 'w' in var_list ):
        fp.write(" <DataArray Name=\"vel\" type=\"Float32\" NumberOfComponents=\"3\">\n")
        uidx = var_list.index('u')
        vidx = var_list.index('v')
        widx = var_list.index('w')
        for k in range(nkc):
            for j in range(njc):
                for i in range(nic):
                    fp.write(" %e %e %e\n" % (uflowz(cell[uidx][i,j,k]),
                                              uflowz(cell[vidx][i,j,k]),
                                              uflowz(cell[widx][i,j,k])))
        
        fp.write(" </DataArray>\n")
    
    fp.write("</CellData>\n")
    fp.write("</Piece>\n")
    fp.write("</UnstructuredGrid>\n")
    fp.write("</VTKFile>\n")
    fp.close()
    return
    
if __name__ == '__main__':
    var_list = create_var_list(nfile)
    print "Reading in Laura grid..."
    nig, njg, nkg, grids = read_laura_grid_file(gname)
    print "Done."
    print "Reading in Laura solution..."
    nif, njf, nkf, data = read_laura_function_file(fname)
    print "Done."
    assert nig == nif
    assert njg == njf
    assert nkg == nkf
    print "Converting node data to cell data..."
    cells = node_data_to_cell_data(nig, njg, nkg, grids, data)
    print "Done."
    print "Writing out VTK file..."
    write_VTK_p_unstructured_file(obfn, nig, njg, nkg, grids, cells, var_list)
    print "Done."
    print "Conversion complete."

