"""
e3_grid.py: Python grid-generation functions for Elmer3.

Builds on the path and basic geometric functions.

.. Author: P.Jacobs

.. Versions:
   17-Mar-2008 extracted from lib/geometry2/source/blockgrid2d.py and 
               lib/geometry2/source/blockgrid3d.py
"""
#----------------------------------------------------------------------

import sys
import math
from numpy import array, zeros
from libprep3 import *
from e3_defs import *
from cfpylib.util.FortranFile import FortranFile

#----------------------------------------------------------------------

def uflowz(q, tiny=1.0e-30):
    """
    Set very small quantities to zero, exactly.

    This is intended primarily to avoid the bad behaviour of VTK
    when it is reading Float32 values that are *too* small.
    We have also come across unreasonably-small float values 
    in the context of reading GridPro files.
    """
    if abs(q) > tiny:
        return q
    else:
        return 0.0

#----------------------------------------------------------------------
# Helper functions for StructuredGrid class

def locate_VTK_header_line(target, f):
    """
    Locate the line containing target and return the tokens on that line.

    Legacy-format VTK lines use spaces as delimiters.
    """
    found = False; tokens = []
    while not found:
        line = f.readline()
        if len(line) == 0: break # end of file
        line = line.strip()
        if target.lower() in line.lower():
            tokens = line.split()
            found = True; break
    if not found: 
        raise RuntimeError("Did not find %s while reading VTK grid file" % target)
    return tokens


class StructuredGrid(object):
    """
    Storage and service functions for a mesh of points defining
    the cell vertices within a block-structured grid.
    """
    def __init__(self, nijk=None, label=None):
        """
        Create an object containg coordinates representing a grid of points.

        :param nijk: optional tuple of integers specifying the number of points 
           in each direction.
           If the number of vertices are specified in each direction,
           we actually create the storage arrays now.
           Note that the number of grid vertices in each index-direction
           will be one more than the number of finite-volume cells
           in that direction.
        :param label: optional string name for the grid
        """
        if nijk != None:
            self.ni = nijk[0]
            self.nj = nijk[1]
            if len(nijk) > 2: 
                self.nk = nijk[2]
            else:
                self.nk = 1
            # print "New StructuredGrid: ni=", self.ni, "nj=", self.nj, "nk=", self.nk
            self.init_arrays()
        else:
            # print "New StructuredGrid: unknown size"
            # The previous comment worried users of the postprocessing code.
            self.x = self.y = self.z = None
        if label == None:
            self.label = ""
        else:
            self.label = label
        return

    def init_arrays(self):
        """
        Create the storage arrays at the previously-specified sizes.
        """
        self.x = zeros((self.ni, self.nj, self.nk), 'd')
        self.y = zeros((self.ni, self.nj, self.nk), 'd')
        self.z = zeros((self.ni, self.nj, self.nk), 'd')
        self.iblank = zeros((self.ni, self.nj, self.nk), 'i')
        return
    
    def write(self, f):
        """
        Writes the grid to an already open file, f.

        This function essentially defines the Eilmer3 native format.
        """
        f.write("%d %d %d  # ni nj nk\n" % (self.ni, self.nj, self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%20.16e %20.16e %20.16e\n" % 
                            (self.x[i,j,k], self.y[i,j,k], self.z[i,j,k]))
        return

    def read(self, f):
        """
        Reads the grid from an already open file, f.
        """
        # First line contains the number of cells in each direction.
        line = f.readline()
        tokens = line.split()
        self.ni = int(tokens[0]); self.nj = int(tokens[1]); self.nk = int(tokens[2])
        self.init_arrays()
        # Coordinates of vertices follow.
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    line = f.readline()
                    tks = line.split()
                    self.x[i,j,k] = float(tks[0])
                    self.y[i,j,k] = float(tks[1])
                    self.z[i,j,k] = float(tks[2])
        return
    
    def read_from_plot3d_whole_grid(self, f, with_blanking=1):
        """
        Read one block from plot3D whole-grid format (ASCII or text file).

        Note that:

           * the file is already opened.
           *  the size of the grid has already been read.

        This format, without blanking, seems to be used by GridGen (from Pointwise).
        """
        print "Start reading plot3D block in whole-grid format..."
        # print "Read x-coordinates"
        np = self.ni * self.nj * self.nk
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.x[i,j,k] = numbers[i+self.ni*(j+self.nj*k)]
        # print "Read y-coordinates"
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.y[i,j,k] = numbers[i+self.ni*(j+self.nj*k)]
        # print "Read z-coordinates"
        numbers = []
        while len(numbers) < np:
            lineContent = f.readline()
            words = lineContent.split()
            for word in words:
                numbers.append(float(word))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    self.z[i,j,k] = numbers[i+self.ni*(j+self.nj*k)]
        if with_blanking:
            # print "Read blanking"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(int(word))
            for k in range(self.nk):
                for j in range(self.nj):
                    for i in range(self.ni):
                        self.iblank[i,j,k] = numbers[i+self.ni*(j+self.nj*k)]
        print "Finished reading block from plot3d whole-grid format."
        return

    def read_from_plot3d_in_planes(self, f, with_blanking=1, verbosity_level=0):
        """
        Read one block from plot3D in-planes (ASCII or text) format.

        Note that:

           * the file is already opened, f is the file object.
           * the size of the grid has already been read.

        This seems to be the format written by ICEM software
        as used by Bianca and Rowan at EPFL, Lausanne.
        The only place that this format seemed to be documented
        is the original Plot3D manual.
        """
        if verbosity_level > 0: print "Start reading plot3D block in planes format..."
        np = self.ni * self.nj
        for k in range(self.nk):
            # print "Plane:", k, "Read x-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.x[i,j,k] = numbers[i+self.ni*j]
            # print "Plane:", k, "Read y-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.y[i,j,k] = numbers[i+self.ni*j]
            # print "Plane:", k, "Read z-coordinates"
            numbers = []
            while len(numbers) < np:
                lineContent = f.readline()
                words = lineContent.split()
                for word in words:
                    numbers.append(float(word))
            for j in range(self.nj):
                for i in range(self.ni):
                    self.z[i,j,k] = numbers[i+self.ni*j]
            if with_blanking:
                # print "Plane:", k, "Read iblanking"
                numbers = []
                while len(numbers) < np:
                    lineContent = f.readline()
                    words = lineContent.split()
                    for word in words:
                        try:
                            numbers.append(int(word))
                        except ValueError:
                            print "Line:", lineContent
                            raise Exception, "Cannot proceed"
                for j in range(self.nj):
                    for i in range(self.ni):
                        self.iblank[i,j,k] = numbers[i+self.ni*j]
            # print "End plane", k
        if verbosity_level > 0: print "Finished reading block from plot3d in planes format."
        return

    def write_structured_plot3d(self, f):
        """
        Write one block in plot3d format to a text file.

        Note that:

           * the file is already opened, f is the file object
           * the headers for the file have already been written to f.
        """
        # Write out all x-coordinates
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%20.16e\n" % self.x[i,j,k])
        # Write out all y-coordinates
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%20.16e\n" % self.y[i,j,k])
        if self.nk > 1:
            # Write out all z-coordinates
            for k in range(self.nk):
                for j in range(self.nj):
                    for i in range(self.ni):
                        f.write("%20.16e\n" % self.z[i,j,k])
        return

    def read_block_in_VTK_format(self, f, verbosity_level=0):
        """
        Reads the grid from an already open file.

        The lines containing some of the metadata at the start of the file
        are read and then ignored.
        """
        # First line contains a declaration that this is a legacy VTK file.
        locate_VTK_header_line("vtk", f) # expecting "vtk DataFile Version 2.0"
        self.label = f.readline().strip()
        if verbosity_level > 0: print "label=", self.label
        locate_VTK_header_line("ASCII", f)
        locate_VTK_header_line("STRUCTURED_GRID", f)
        tokens = locate_VTK_header_line("DIMENSIONS", f)
        self.ni = int(tokens[1]); self.nj = int(tokens[2]); self.nk = int(tokens[3])
        tokens = locate_VTK_header_line("POINTS", f)
        total_points = int(tokens[1])
        if total_points != self.ni * self.nj * self.nk:
            print "label=", self.label
            print "POINTS=", total_points, "; ni=", ni, "nj=", nj, "nk=", nk
            raise RuntimeError("In VTK STRUCTURED_GRID file: points mismatch.")
        # Proceed with allocating arrays and filling with point data.
        self.init_arrays()
        # Coordinates of vertices follow.
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    line = f.readline()
                    tks = line.split()
                    self.x[i,j,k] = float(tks[0])
                    self.y[i,j,k] = float(tks[1])
                    self.z[i,j,k] = float(tks[2])
        return

    def write_block_in_VTK_format(self, f, verbosity_level=0):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        Note also that this is the legacy VTK format.
        """
        if verbosity_level > 0: print "Begin write block in VTK format:"
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.ni, self.nj, self.nk))
        f.write("POINTS %d float\n" % (self.ni * self.nj * self.nk))
        for k in range(self.nk):
            for j in range(self.nj):
                for i in range(self.ni):
                    f.write("%e %e %e\n" % (self.x[i,j,k], self.y[i,j,k], self.z[i,j,k]))
        if verbosity_level > 0: print "End write block in VTK format."
        return

    def make_grid_from_surface(self, surface, cluster_functions=[None,]*4, verbosity_level=0):
        """
        Given a parametric surface, create a 2D grid via interpolation.

        The order of the cluster functions is N,E,S,W.
        """
        if verbosity_level > 0: print "Begin make grid"
        if self.nk != 1:
            print "Warning: nk=", self.nk, "but should be zero."
            print "The 2D grid will be put in the k=0 slice."
        # Set up distributions of points along each of the nondimensional edges.
        # If a cluster function has not been supplied, default to a linear distribution.
        for icf in range(4):
            if not isinstance(cluster_functions[icf], UnivariateFunction):
                cluster_functions[icf] = LinearFunction()
        rNorth = cluster_functions[0].distribute_parameter_values(self.ni)
        sEast = cluster_functions[1].distribute_parameter_values(self.nj)
        rSouth = cluster_functions[2].distribute_parameter_values(self.ni)
        sWest = cluster_functions[3].distribute_parameter_values(self.nj)
        # Now, work through the mesh, one point at a time,
        # blending the stretched parameter values
        # and creating the actual vertex coordinates in Cartesian space.
        k = 0
        for j in range(self.nj):
            s = float(j) / (self.nj - 1)
            for i in range(self.ni):
                r = float(i) / (self.ni - 1)
                sdash = (1.0-r) * sWest[j] + r * sEast[j] 
                rdash = (1.0-s) * rSouth[i] + s * rNorth[i]
                p = surface.eval(rdash, sdash)
                self.x[i,j,k], self.y[i,j,k], self.z[i,j,k] = p.x, p.y, p.z
            if verbosity_level > 0: print ".",
        if verbosity_level > 0:
            print 
            print "End make 2D grid from surface."
        return

    def make_TFI_grid_from_volume(self, pvolume, cluster_functions=[None,]*12,
                                  verbosity_level=0):
        """
        Given a parametric volume, create the grid via TFI.
        
        The clustering information always comes from the edges.
        A full compliment of 12 should be supplied but, 
        if any are not already specified correctly, 
        they will be set as standard linear functions.
        """
        if verbosity_level > 0: print "Begin make grid."
        # Set up distributions of points along each of the nondimensional edges.
        for i in range(12):
            if not isinstance(cluster_functions[i], UnivariateFunction):
                cluster_functions[i] = LinearFunction()
        #
        r01 = cluster_functions[0].distribute_parameter_values(self.ni)
        r32 = cluster_functions[2].distribute_parameter_values(self.ni)
        s12 = cluster_functions[1].distribute_parameter_values(self.nj)
        s03 = cluster_functions[3].distribute_parameter_values(self.nj)
        #
        r45 = cluster_functions[4].distribute_parameter_values(self.ni)
        r76 = cluster_functions[6].distribute_parameter_values(self.ni)
        s56 = cluster_functions[5].distribute_parameter_values(self.nj)
        s47 = cluster_functions[7].distribute_parameter_values(self.nj)
        #
        t04 = cluster_functions[8].distribute_parameter_values(self.nk)
        t15 = cluster_functions[9].distribute_parameter_values(self.nk)
        t26 = cluster_functions[10].distribute_parameter_values(self.nk)
        t37 = cluster_functions[11].distribute_parameter_values(self.nk)
        #
        # Now, work through the mesh, blending the stretched parameter values
        # and creating the actual vertex coordinates in Cartesian space.
        for k in range(self.nk):
            t = float(k) / (self.nk - 1)
            for j in range(self.nj):
                s = float(j) / (self.nj - 1)
                for i in range(self.ni):
                    r = float(i) / (self.ni - 1)
                    tdash = (1.0-r)*(1.0-s)*t04[k] + r*s*t26[k] + \
                            (1.0-s)*r*t15[k] + s*(1.0-r)*t37[k]
                    sdash = (1.0-t)*(1.0-r)*s03[j] + t*r*s56[j] + \
                            (1.0-t)*r*s12[j] + t*(1-r)*s47[j]
                    rdash = (1.0-s)*(1.0-t)*r01[i] + s*t*r76[i] + \
                            (1.0-s)*t*r45[i] + s*(1.0-t)*r32[i]
                    p = pvolume.eval(rdash, sdash, tdash)
                    self.x[i][j][k], self.y[i][j][k], self.z[i][j][k] = p.x, p.y, p.z
            if verbosity_level > 0: print ".",
            sys.stdout.flush()
        if verbosity_level > 0:
            print 
            print "End make grid."
        return

    def create_subgrid(self, imin, imax, jmin, jmax, kmin=0, kmax=0):
        """
        Returns a new StructuredGrid that is a subgrid of the current grid.

        The limits (imin,imax) are included in the new subgrid.
        """
        ni = imax - imin + 1
        nj = jmax - jmin + 1
        nk = kmax - kmin + 1
        subgrid = StructuredGrid((ni,nj,nk))
        #
        ksub = 0
        for k in range(kmin, kmax+1):
            jsub = 0
            for j in range(jmin, jmax+1):
                isub = 0
                for i in range(imin, imax+1):
                    subgrid.x[isub,jsub,ksub] = self.x[i,j,k]
                    subgrid.y[isub,jsub,ksub] = self.y[i,j,k]
                    subgrid.z[isub,jsub,ksub] = self.z[i,j,k]
                    isub += 1
                jsub += 1
            ksub += 1
        #
        return subgrid
    
    def get_vertex_coords(self, ivtx):
        """
        Returns a tuple of coordinates for a single vertex.

        The indexing for ivtx corresponds to the VTK convention
        for a hexahedral cell (scaled up to the whole block).
        """
        if ivtx == 0:
            i = 0; j = 0; k = 0
        elif ivtx == 1:
            i = self.ni-1; j = 0; k = 0
        elif ivtx == 2:
            i = self.ni-1; j = self.nj-1; k = 0
        elif ivtx == 3:
            i = 0; j = self.nj-1; k = 0
        elif ivtx == 4:
            i = 0; j = 0; k = self.nk-1
        elif ivtx == 5:
            i = self.ni-1; j = 0; k = self.nk-1
        elif ivtx == 6:
            i = self.ni-1; j = self.nj-1; k = self.nk-1
        elif ivtx == 7:
            i = 0; j = self.nj-1; k = self.nk-1
        else:
            raise ValueError, ("vertex index: %d" % ivtx)
        return Vector3(self.x[i,j,k], self.y[i,j,k], self.z[i,j,k])

    def get_vertex_list_for_cell(self, i, j, k=0):
        """
        Returns a list of Vectors defining the vertices of the cell.

        See block.hh for the root definition.
        All 8 points are available in a 3D grid but only p0 through p3
        are available in a 2D grid.
        """
        vtxList = [None,] * 8
        vtxList[0] = Vector(self.x[i,j,k], self.y[i,j,k], self.z[i,j,k])
        vtxList[1] = Vector(self.x[i+1,j,k], self.y[i+1,j,k], self.z[i+1,j,k])
        vtxList[2] = Vector(self.x[i+1,j+1,k], self.y[i+1,j+1,k], self.z[i+1,j+1,k])
        vtxList[3] = Vector(self.x[i,j+1,k], self.y[i,j+1,k], self.z[i,j+1,k])
        if self.nk > 1:
            vtxList[4] = Vector(self.x[i,j,k+1], self.y[i,j,k+1], self.z[i,j,k+1])
            vtxList[5] = Vector(self.x[i+1,j,k+1], self.y[i+1,j,k+1], self.z[i+1,j,k+1])
            vtxList[6] = Vector(self.x[i+1,j+1,k+1], self.y[i+1,j+1,k+1], self.z[i+1,j+1,k+1])
            vtxList[7] = Vector(self.x[i,j+1,k+1], self.y[i,j+1,k+1], self.z[i,j+1,k+1])
        return vtxList
    
    def point_inside_cell(self, i, j, k, x, y, z, dimensions):
        """
        Similar to the point_inside_cell function in cell.cxx
        Returns 1 if the point p is inside or on the cell surface.
        """
        if dimensions == 2:
            # In 2 dimensions,
            # we split the x,y-plane into half-planes and check which side p is on.
            xA = self.x[i+1,j,k]     # A is vtx[1]
            yA = self.y[i+1,j,k] 
            xB = self.x[i+1,j+1,k]   # B is vtx[2] 
            yB = self.y[i+1,j+1,k]
            xC = self.x[i,j+1,k]     # C is vtx[3]
            yC = self.y[i,j+1,k]
            xD = self.x[i,j,k]       # D is vtx[0]
            yD = self.y[i,j,k]
            # Now, check to see if the specified point is on the
            # left of (or on) each boundary line AB, BC, CD and DA.
            if ((x - xB) * (yA - yB) >= (y - yB) * (xA - xB) and
                (x - xC) * (yB - yC) >= (y - yC) * (xB - xC) and
                (x - xD) * (yC - yD) >= (y - yD) * (xC - xD) and
                (x - xA) * (yD - yA) >= (y - yA) * (xD - xA)):
                return 1
            else:
                return 0
        else:
            # 3D
            print "This cell locating function does not function for 3D geometries yet."
            print "You should probably be using Dan Potter's suggest_better_cell()."
            raise RuntimeError("[TODO] Not ready for 3D yet")
        
    def suggest_better_cell(self, i, j, k, x, y, z, dimensions):
        """
        Using the method from cell_finder.cxx (ray-tracing) where a cell step 
        is returned (ie di = change in i index, etc) based on what side of the 
        interfaces the point is on.

        .. Author: Dan (?)
        """
        if dimensions == 2:
            # In 2 dimensions,
            # we split the x,y-plane into half-planes and check which side p is on.
            # 1. define the vertices (see block.hh)
            xA = self.x[i+1,j,k]     # A is vtx[1]
            yA = self.y[i+1,j,k] 
            xB = self.x[i+1,j+1,k]   # B is vtx[2] 
            yB = self.y[i+1,j+1,k]
            xC = self.x[i,j+1,k]     # C is vtx[3]
            yC = self.y[i,j+1,k]
            xD = self.x[i,j,k]       # D is vtx[0]
            yD = self.y[i,j,k]
            # 2. calculate the cell index increments
            di = 0; dj = 0; dk = 0
            if (x - xB) * (yA - yB) >= (y - yB) * (xA - xB): di -= 1
            if (x - xD) * (yC - yD) >= (y - yD) * (xC - xD): di += 1
            if (x - xC) * (yB - yC) >= (y - yC) * (xB - xC): dj -= 1
            if (x - xA) * (yD - yA) >= (y - yA) * (xD - xA): dj += 1
            # 4. return the cell index increments
            return di, dj, dk
        else:
            # In 3 dimensions,
            # we use an algorithm borrowed from CellFinder3D::test_cell()
            # 0. Set geometric book-keeping parameters (NOTE: assuming hex cells)
            nvtx = 8; nfaces = 6; npoints = 3
            hex_vertex_indices = [ [ 2, 6, 7 ],
                                   [ 5, 6, 2 ],
                                   [ 0, 4, 5 ],
                                   [ 0, 3, 7 ],
                                   [ 7, 6, 5 ],
                                   [ 0, 1, 2 ] ]
            # 1. define the vertices (see block.hh)
            vtx = self.get_vertex_list_for_cell(i, j, k)
            # 2. calculate the 'a' parameter that defines what side of the face the point is
            a = [ None ] * nfaces; vp = [ None ] * npoints
            p = Vector3( x, y, z )
            for iface in range(nfaces):
                for i in range(npoints):
                    vp[i] = vtx[ hex_vertex_indices[iface][i] ] - p
                a[iface] = dot( vp[0], cross( vp[1], vp[2] ) )
            # 3. Determine cell index increments from the a dot products
            di = 0; dj = 0; dk = 0
            if a[EAST] < 0.0: di -= 1
            if a[WEST] < 0.0: di += 1
            if a[NORTH] < 0.0: dj -= 1
            if a[SOUTH] < 0.0: dj += 1
            if a[TOP] < 0.0: dk -= 1
            if a[BOTTOM] < 0.0: dk += 1
            # 4. return the cell index increments * scale
            return di, dj, dk
        # end method suggest_better_cell()

    def cell_centre_position(self, i, j, k):
        """
        Given the indices for the cell, Cartesian coordinates of the centre.

        Note that the computation of Block2D.cell_centre_location(i,j,k,gdata)
        over in e3_block.py is a little different (and more involved).
        """
        vtx = self.get_vertex_list_for_cell(i, j, k)
        if self.nk == 1:
            cell_centre = 0.25*(vtx[0]+vtx[1]+vtx[2]+vtx[3])
        else:
            assert self.nk > 1
            cell_centre = 0.125*(vtx[0]+vtx[1]+vtx[2]+vtx[3]+vtx[4]+vtx[5]+vtx[6]+vtx[7])
        return cell_centre.x, cell_centre.y, cell_centre.z
        

    def compute_ghost_cell_positions(self, i, j, k, which_face):
        """
        Given a specific cell on a block boundary, compute the cell-centre
        positions of the adjacent ghost cells, just outside the boundary.
        The ghost cell positions are a continuation of the interior-cell
        positions along the ling through the boundary face. 
        A little like a mirror-image but ignoring the normal direction
        of the face.
        """
        vtx = self.get_vertex_list_for_cell(i, j, k)
        if self.nk == 1:
            cell1_centre = 0.25*(vtx[0]+vtx[1]+vtx[2]+vtx[3])
            if which_face == NORTH:
                face_centre = 0.5*(vtx[3]+vtx[2])
                vtx2 = self.get_vertex_list_for_cell(i, j-1, k)
            elif which_face == EAST:
                face_centre = 0.5*(vtx[1]+vtx[2])
                vtx2 = self.get_vertex_list_for_cell(i-1, j, k)
            elif which_face == SOUTH:
                face_centre = 0.5*(vtx[0]+vtx[1])
                vtx2 = self.get_vertex_list_for_cell(i, j+1, k)
            elif which_face == WEST:
                face_centre = 0.5*(vtx[3]+vtx[0])
                vtx2 = self.get_vertex_list_for_cell(i+1, j, k)
            else:
                raise RuntimeError("compute_ghost_cell_positions(): "
                                   "incorrect face value: %s" % which_face)
            cell2_centre = 0.25*(vtx2[0]+vtx2[1]+vtx2[2]+vtx2[3])
        else:
            assert self.nk > 1
            cell1_centre = 0.125*(vtx[0]+vtx[1]+vtx[2]+vtx[3]+vtx[4]+vtx[5]+vtx[6]+vtx[7])
            if which_face == NORTH:
                face_centre = 0.25*(vtx[3]+vtx[2]+vtx[6]+vtx[7])
                vtx2 = self.get_vertex_list_for_cell(i, j-1, k)
            elif which_face == EAST:
                face_centre = 0.25*(vtx[1]+vtx[2]+vtx[6]+vtx[5])
                vtx2 = self.get_vertex_list_for_cell(i-1, j, k)
            elif which_face == SOUTH:
                face_centre = 0.25*(vtx[0]+vtx[1]+vtx[5]+vtx[4])
                vtx2 = self.get_vertex_list_for_cell(i, j+1, k)
            elif which_face == WEST:
                face_centre = 0.25*(vtx[0]+vtx[4]+vtx[7]+vtx[3])
                vtx2 = self.get_vertex_list_for_cell(i+1, j, k)
            elif which_face == TOP:
                face_centre = 0.25*(vtx[4]+vtx[5]+vtx[6]+vtx[7])
                vtx2 = self.get_vertex_list_for_cell(i, j, k-1)
            elif which_face == BOTTOM:
                face_centre = 0.25*(vtx[0]+vtx[1]+vtx[2]+vtx[3])
                vtx2 = self.get_vertex_list_for_cell(i, j, k+1)
            else:
                raise RuntimeError("compute_ghost_cell_positions(): "
                                   "incorrect face value: %s" % which_face)
            cell2_centre = 0.125*(vtx2[0]+vtx2[1]+vtx2[2]+vtx2[3]+vtx2[4]+vtx2[5]+vtx2[6]+vtx2[7])
        delta1 = face_centre - cell1_centre
        ghost1 = face_centre + delta1
        delta2 = face_centre - cell2_centre
        ghost2 = face_centre + delta2
        return ghost1.x, ghost1.y, ghost1.z, ghost2.x, ghost2.y, ghost2.z

#--------------------------------------------------------------------

def write_plot3d_grid(fname, grid):
    """
    Writes a formatted multiple-block Plot3D grid (ASCII or text) file.

    :param fname: name of the plot3d grid file to create
    :param grid: list of StructuredGrid objects
    :returns: None
    """
    f = open(fname, 'w')
    f.write(" %d\n" % len(grid))
    for i in range(len(grid)):
        if grid[0].nk == 1:
            f.write(" %d %d\n" % (grid[i].ni, grid[i].nj))
        else:
            f.write(" %d %d %d\n" % (grid[i].ni, grid[i].nj, grid[i].nk))
    for g in grid:
        g.write_structured_plot3d(f)
    f.close()
    return

def read_plot3d_grid(fname, dimensions):
    """
    Reads a plot3d grid file, containing one or more block-grids.

    :param fname: name of plot3d grid file
    :param dimensions: can have the value 2 or 3
    :returns: a list StructuredGrid object(s)

    Assumptions:

       * text file, Fortran formatted
       * multiple-block
       * with I-blanking
    """
    f = open(fname, 'r')
    tokens = f.read().split()
    f.close()
    ngrids = int(tokens[0])
    #
    pos = 1
    gridList = []
    for i in range(ngrids):
        ni = int(tokens[pos]); pos += 1
        nj = int(tokens[pos]); pos += 1
        nijk = [ni, nj]
        nk = 1
        if dimensions == 3:
            nk = int(tokens[pos]); pos += 1
            nijk.append(nk)
        # Initialise grid object
        gridList.append(StructuredGrid(nijk))
        # Read x-coordinates
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    gridList[-1].x[i,j,k] = float(tokens[pos]); pos += 1
        # Read y-coordinates
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    gridList[-1].y[i,j,k] = float(tokens[pos]); pos += 1
        # Read z-coordinates (if 3D)
        if dimensions == 3:
            for k in range(nk):
                for j in range(nj):
                    for i in range(ni):
                        gridList[-1].z[i,j,k] = float(tokens[pos]); pos += 1
        # Read (and discard) i-blanking stuff
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    pos += 1
    #
    return gridList

def read_plot3d_grid_from_pointwise(fname, dim, scaleFactor):
    """
    Reads a plot3d grid, returns a list StructuredGrid object(s).
       
    :param fname: name of plot3d grid file
    :param dim: has the value 2 or 3
    :param scaleFactor: allow for correction of dimensions of the grid, if required.
       For example, from millimetres to metres, scaleFactor would be 1000.
    :returns: list of StructuredGrid objects 
    
    Assumptions:

       * formatted
       * multiple-block
       
    .. Modified from 'read_plot3d_grid' for use with POINTWISE plot3d grids.
    .. Authors: LukeD/WilsonC, 02-June-2011
    """
    #
    f = open(fname, 'r')
    tks = f.read().split()
    f.close()
    ngrids = int(tks[0])
    #
    pos = 1
    gridList = []
    nijkList = []
    for i in range(ngrids):
        ni = int(tks[pos]); pos += 1
        nj = int(tks[pos]); pos += 1
        nijkList.append([ni, nj])
        nk = 1
        if dim == 3:
            nk = int(tks[pos]); pos += 1
            nijkList[i].append(nk)
    #
    for indx in range(ngrids):
        #
        # Initialise grid object
        gridList.append(StructuredGrid(nijkList[indx]))
        #
        # Read x-coordinates
        for k in range(nijkList[indx][2]):
            for j in range(nijkList[indx][1]):
                for i in range(nijkList[indx][0]):
                    gridList[indx].x[i,j,k] = float(tks[pos])/scaleFactor; pos += 1
        #
        # Read y-coordinates
        for k in range(nijkList[indx][2]):
            for j in range(nijkList[indx][1]):
                for i in range(nijkList[indx][0]):
                    gridList[indx].y[i,j,k] = float(tks[pos])/scaleFactor; pos += 1
        #
        # Read z-coordinates (if 3D)
        if dim == 3:
            for k in range(nijkList[indx][2]):
                for j in range(nijkList[indx][1]):
                    for i in range(nijkList[indx][0]):
                        gridList[indx].z[i,j,k] = float(tks[pos])/scaleFactor; pos += 1
        #
    #
    return gridList

#--------------------------------------------------------------------

def read_plot3d_grid_binary_file(gname):
    """
    Returns a lists of mesh dimensions and coordinate data read from a single plot3d file.

    This is Rowan's function for reading files from Peter Gnoffo's Laura code.

    :param gname: name of plot3d grid file
    :returns: a tuple containing lists of the i,j and k dimensions of eash block
        and a list of the StructuredGrid objects holding the grid points.

    .. Given the stories about the variations of plot3d formet, it's probably not
       worth the effort of trying to merge this function into the previous
       read_plot3d_grid function.  There seem to be a lot of special/arbitrary cases.
    """
    g = FortranFile(gname)
    nblocks = g.readInts()[0]
    ndim = g.readInts()
    ni = []; nj = []; nk = []
    for m in range(nblocks):
        ni.append(ndim[3*m])
        nj.append(ndim[3*m+1])
        nk.append(ndim[3*m+2])
    #
    pts = g.readReals('d')
    grids = []
    pos = 0
    for m in range(nblocks):
        grids.append(StructuredGrid((ni[m], nj[m], nk[m])))
        # Read x values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m].x[i,j,k] = pts[pos]
                    pos = pos + 1
        # Read y values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m].y[i,j,k] = pts[pos]
                    pos = pos + 1
        # Read z values
        for k in range(nk[m]):
            for j in range(nj[m]):
                for i in range(ni[m]):
                    grids[m].z[i,j,k] = pts[pos]
                    pos = pos + 1
    #
    return (ni, nj, nk, grids)

def read_gridpro_grid(gname, scale=1.0):
    """
    Reads a complete Gridpro grid file, returns a list of StructuredGrids.

    A complete Gridpro grid file contains multiple blocks. This function
    will read through all blocks and store them as StructuredGrid objects.
    These are returned by the function. Gridpro builds grids in the same
    dimensions as the supplied geometry. Care should be taken with 
    Gridpro grids built from CAD geometries which may typically be
    in millimetres. In this case, the required 'scale' would be 0.001
    to convert to metres for use in Eilmer.
    
    :param gname: name of Gridpro grid file
    :param scale: a scale to convert supplied coordinates to metres
    :returns: list of StructuredGrid object(s)
    
    .. Author: Rowan J. Gollan
    .. Date: 16-Aug-2012
    """
    f = open(gname, 'r')
    grids = []
    while True:
        line = f.readline().strip()
        if len(line) == 0: break
        if line[0] == '#': continue
        tks = line.split()
        if not tks: break # Presumably reached end of file
        ni, nj, nk = int(tks[0]), int(tks[1]), int(tks[2])
        grids.append(StructuredGrid([ni, nj, nk]))
        for i in range(ni):
            for j in range(nj):
                for k in range(nk):
                    tks = f.readline().split()
                    grids[-1].x[i,j,k] = uflowz(scale*float(tks[0]))
                    grids[-1].y[i,j,k] = uflowz(scale*float(tks[1]))
                    grids[-1].z[i,j,k] = uflowz(scale*float(tks[2]))
    f.close()
    return grids

#--------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin demo of e3_grid.py."

    print "Generate a grid in a box."
    b0 = StructuredGrid((3,4,5), label="Some box")

    p0 = Node(0.0, 0.0, 0.0); p1 = Node(2.0, 0.0, 0.0)
    p2 = Node(2.0, 1.0, 0.0); p3 = Node(0.0, 1.0, 0.0)
    p4 = Node(0.0, 0.0, 0.5); p5 = Node(2.0, 0.0, 0.5)
    p6 = Node(2.0, 1.0, 0.5); p7 = Node(0.0, 1.0, 0.5)
    boxvol = SimpleBoxVolume(p0, p1, p2, p3, p4, p5, p6, p7, "BOX")
    clusterf = [None,]*12
    clusterf[3] = RobertsClusterFunction(0,1,1.1)
    clusterf[5] = RobertsClusterFunction(0,1,1.1)
    b0.make_TFI_grid_from_volume(boxvol, clusterf)

    print "Write to VTK file."
    fout = open("test_block3d.vtk", "w")
    b0.write_block_in_VTK_format(fout)
    fout.close()

    print "Read from VTK file."
    fin = open("test_block3d.vtk", "r")
    b1 = StructuredGrid()
    b1.read_block_in_VTK_format(fin)
    fin.close()

    print "Extract corners of block."
    for ivtx in range(8):
        print "vertex", ivtx, ":", b1.get_vertex_coords(ivtx)
    
    print "Done."

