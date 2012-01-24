## \file blockgrid2d.py
## \ingroup libgeom2
## \brief Python grid-generation functions for MB_CNS.
##
## \author P.Jacobs
## \version 28-Nov-2005 extracted from gpath.py, scriptit.py as geom_mbcns.py
## \version 15-Jan-2006 now blockgrid2d.py and using libgeom2
"""
Python grid-generation functions for MB_CNS.

Builds on the path and basic geometric functions.
"""
#----------------------------------------------------------------------

import sys
import math
try:
    from numpy import array, zeros
except:
    try:
        from Numeric import array, zeros
    except:
        print "Could import neither numpy nor Numeric."

if __name__ == 'blockgrid2d':
    # then we need to import the geometry routines
    from libgeom2 import *

else:
    #  Assume that this code is part of an execfile call and that
    #  the geometry routines are already loaded.
    print "In file blockgrid2d.py: we expect that geometry2 routines"
    print "are already loaded"
    pass

#----------------------------------------------------------------------

class BlockGrid2D(object):
    """
    Storage and service functions for a mesh of points defining
    the cell vertices within a block-structured grid.
    """
    def __init__(self, ni=None, nj=None, label=None):
        """
        @param ni: number of points in the i-index direction
        @type ni: int
        @param nj: number of points in the j-index direction
        @type nj: int
        @param label: string label for the block
        @type label: string
        
        @note: If the number of vertices are specified in each direction,
            we actually create the storage arrays now.
        @note: The number of grid vertices in each index-direction
            will be one more than the number of finite-volume cells
            in that direction.
        """
        self.ni = ni
        self.nj = nj
        if ni != None and nj != None:
            print "New BlockGrid2D: ni=", ni, "nj=", nj
            self.init_arrays()
        else:
            print "New BlockGrid2D: unknown size"
            self.x = self.y = None
        self.label = label
        return

    def init_arrays(self):
        """
        Create the storage arrays at the previously-specified sizes.
        """
        self.x = zeros((self.ni, self.nj), 'd')
        self.y = zeros((self.ni, self.nj), 'd')
        return
    
    def write_block_in_VTK_format(self, f, header_flag=True):
        """
        Writes the grid to an already open file.
        
        """
        print "Begin write block: label=", self.label
        if header_flag:
            f.write("# vtk DataFile Version 2.0\n")
            f.write("%s\n" % self.label)
            f.write("ASCII\n")
            f.write("\n")
            
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d 1\n" % (self.ni, self.nj))
        f.write("POINTS %d float\n" % (self.ni * self.nj))
        
        for j in range(self.nj):
            for i in range(self.ni):
                f.write("%e %e 0.0\n" % (self.x[i][j], self.y[i][j]))

        print "End write block."
        return

    def write_block_in_tecplot_format(self, f, header_flag=True):
        """
        Writes the grid to an already open file.
        
        """
        print "Begin write block: label=", self.label
        if header_flag:
            f.write("TITLE = \"%s\"" % self.label)
            f.write("VARIABLES = \"X\", \"Y\"",)
        f.write("ZONE I=%i, J=%i, DATAPACKING=BLOCK" % (self.ni, self.nj))

        for i in range(self.ni):
            for j in range(self.nj):
                f.write("%e " % (self.x[i][j]))
            f.write("\n")

        for i in range(self.ni):
            for j in range(self.nj):
                f.write("%e " %(self.y[i][j]))
            f.write("\n")

        print "End write block."
        return

    def write_block_in_classic_mbcns_format(self, f):
        """
        Writes the grid (in old-style format) to an already open file.
        """
        f.write("%d %d\n" % (self.ni-1, self.nj-1)) # number of cells in each dir
        for i in range(self.ni):
            for j in range(self.nj):
                f.write("%20.12e %20.12e\n" % (self.x[i][j], self.y[i][j]))
        return

    def read_block_in_classic_mbcns_format(self, f):
        """
        Reads the grid (in old-style format) from an already open file.
        """
        # First line contains the number of cells in each direction.
        line = f.readline()
        tks = line.split()
        ncelli = int(tks[0])
        ncellj = int(tks[1])
        if self.ni is None or self.nj is None or \
               self.ni != ncelli or self.nj != ncellj:
            self.ni = ncelli + 1
            self.nj = ncellj + 1
            self.init_arrays()
        # Coordinates of vertices follow.
        for i in range(self.ni):
            for j in range(self.nj):
                line = f.readline()
                tks = line.split()
                self.x[i][j] = float(tks[0])
                self.y[i][j] = float(tks[1])
        return

    def read(self, f):
        """
        Alias for read_block_in_classic_mbcns_format.
        """
        self.read_block_in_classic_mbcns_format(f)
        return
    
    def make_grid_from_surface(self, surface, cluster_functions=[None,]*4):
        """
        Given a parametric surface, create the grid via interpolation.
        """
        print "Begin make grid, label=", self.label
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
        for j in range(self.nj):
            s = float(j) / (self.nj - 1)
            for i in range(self.ni):
                r = float(i) / (self.ni - 1)
                sdash = (1.0-r) * sWest[j] + r * sEast[j] 
                rdash = (1.0-s) * rSouth[i] + s * rNorth[i]
                p = surface.eval(rdash, sdash)
                self.x[i][j], self.y[i][j] = p.x, p.y
            print ".",
        print 
        # print "End make grid."
        return

    def create_subgrid(self, imin, imax, jmin, jmax):
        ni = imax - imin + 1
        nj = jmax - jmin + 1
        subgrid = BlockGrid2D(ni, nj)

        isub = 0
        for i in range(imin, imax+1):
            jsub = 0
            for j in range(jmin, jmax+1):
                subgrid.x[isub][jsub] = self.x[i][j]
                subgrid.y[isub][jsub] = self.y[i][j]
                jsub += 1
            isub += 1

        return subgrid


    def create_north_spline(self):
        nsp = 5

        if( self.ni < nsp ):
            nsp = self.ni

        points = [ Vector3( self.x[0][self.nj-1],  self.y[0][self.nj-1] ) ]
        nx = self.ni - 1
        i = 0
        for isp in range(1,nsp-1):
            i += nx / (nsp - isp)
            points.append( Vector3( self.x[i][self.nj-1], self.y[i][self.nj-1]) )
            nx -= i
        points.append( Vector3(self.x[self.ni-1][self.nj-1],  self.y[self.ni-1][self.nj-1]) )

        return Spline(points)

    def create_east_spline(self):
        nsp = 5

        if( self.nj < nsp ):
            nsp = self.nj

        points = [ Vector3( self.x[self.ni-1][0],  self.y[self.ni-1][0] ) ]
        ny = self.nj - 1
        j = 0
        for isp in range(1,nsp-1):
            j += ny / (nsp - isp)
            points.append( Vector3( self.x[self.ni-1][j], self.y[self.ni-1][j]) )
            ny -= j
        points.append( Vector3(self.x[self.ni-1][self.nj-1],  self.y[self.ni-1][self.nj-1]) )

        return Spline(points)

    def create_south_spline(self):
        nsp = 5

        if( self.ni < nsp ):
            nsp = self.ni

        points = [ Vector3( self.x[0][0],  self.y[0][0] ) ]
        nx = self.ni - 1
        i = 0
        for isp in range(1,nsp-1):
            i += nx / (nsp - isp)
            points.append( Vector3( self.x[i][0], self.y[i][0]) )
            nx -= i
        points.append( Vector3(self.x[self.ni-1][0],  self.y[self.ni-1][0]) )

        return Spline(points)

    def create_west_spline(self):
        nsp = 5

        if( self.nj < nsp ):
            nsp = self.nj

        points = [ Vector3( self.x[0][0],  self.y[0][0] ) ]
        ny = self.nj - 1
        j = 0
        for isp in range(1,nsp-1):
            j += ny / (nsp - isp)
            points.append( Vector3( self.x[0][j], self.y[0][j]) )
            ny -= j
        points.append( Vector3(self.x[0][self.nj-1],  self.y[0][self.nj-1]) )

        return Spline(points)

    def points_in_VTK_StructuredGrid_order(self):
        points = []
        for j in range(self.nj):
            for i in range(self.ni):
                points.append([self.x[i][j], self.y[i][j], 0.0])
        return points

    

            

        
#--------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin demo of blockgrid2d.py."
    #
    print "Generate a grid in a box."
    b0 = BlockGrid2D(30, 40, label="Some box")
    p0 = Node(0.0,0.0)
    p1 = Node(2.0,0.0)
    p2 = Node(2.0,3.0)
    p3 = Node(0.0,3.0)
    clusterf = [None, RobertsClusterFunction(0,1,1.1),
                None, RobertsClusterFunction(0,1,1.1)]
    b0.make_grid_from_surface(CoonsPatch(p0, p1, p2, p3, "TEST-GRID"), clusterf)
    #
    print "Write to VTK file."
    fout = open("test_block2d.vtk", "w")
    b0.write_block_in_VTK_format(fout)
    fout.close()
    #
    print "Done."
