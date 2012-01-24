## \file geom_mbcns.py
## \ingroup geom
## \brief Python geometry-specification functions for MB_CNS.
##
## \author P.Jacobs
## \version 28-Nov-2005 extracted from gpath.py, scriptit.py
"""
Face2D definition and grid generation functions for MB_CNS.

Builds on the path and basic geometric functions of libgeom, geom, gpath.
"""
#----------------------------------------------------------------------

# print "Loading gpath.py"
try:
    import libgeom
    # print "Loaded libgeom successfully."
except:
    print "Did not load libgeom C module successfully."

import sys
import math
import Numeric
from geom import *
from gpath import *
from cluster_points import *

#----------------------------------------------------------------------

def interpolate_TFI_2D_mbcns(edge_list, r, s):
    """
    Locate a point p(r,s) using transfinite interpolation.

    @param edge_list: List of 4 paths in mb_cns order [N,E,S,W].
    @type edge_list: list of L{Polyline}-derived objects.
    @param r: interpolation parameter in the ix-index direction
    @type r: float, 0.0<=r<=1.0
    @param s: interpolation parameter in the iy-index direction
    @type s: float, 0.0<=s<=1.0
    """
    assert isinstance(edge_list, list), "Should be a list"
    assert len(edge_list) == 4, "Should have 4 edges."
    for edge in edge_list:
        assert isinstance(edge, Polyline), "edge should be at least a Polyline."
    c2north = edge_list[0].plp;
    c4east  = edge_list[1].plp;
    c1south = edge_list[2].plp;
    c3west  = edge_list[3].plp;
    p = Vector()
    flag = libgeom.coons_patch(c1south, c2north, c3west, c4east,
                               r, s, p.pp)
    return (p.x, p.y, p.z)

#---------------------------------------------------------------------

class BlockGrid2D(object):
    """
    Storage and service functions for a mesh of points defining
    the cell vertices within a block-structured grid.
    """
    def __init__(self, nix=None, niy=None, label=None):
        """
        @param nix: number of points in the ix-index direction
        @type nix: int
        @param niy: number of points in the iy-index direction
        @type niy: int
        @param label: string label for the block
        @type label: string
        
        @note: If the number of vertices are specified in each direction,
            we actually create the storage arrays now.
        @note: The number of grid vertices in each index-direction
            will be one more than the number of finite-volume cells
            in that direction.
        """
        self.nix = nix
        self.niy = niy
        if nix != None and niy != None:
            print "New BlockGrid2D: nix=", nix, "niy=", niy
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
        self.x = Numeric.zeros((self.nix, self.niy), 'd')
        self.y = Numeric.zeros((self.nix, self.niy), 'd')
        return
    
    def write_block_in_VTK_format(self, f):
        """
        Writes the grid to an already open file.

        Note that this function writes all of the VTK header lines so
        that it is implicit that only one block grid goes into the file.
        """
        print "Begin write block: label=", self.label
        f.write("# vtk DataFile Version 2.0\n")
        f.write("%s\n" % self.label)
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d 1\n" % (self.nix, self.niy))
        f.write("POINTS %d float\n" % (self.nix * self.niy))
        for iy in range(self.niy):
            for ix in range(self.nix):
                f.write("%e %e 0.0\n" % (self.x[ix][iy], self.y[ix][iy]))
        print "End write block."
        return

    def write_block_in_classic_mbcns_format(self, f):
        """
        Writes the grid (in old-style format) to an already open file.
        """
        f.write("%d %d\n" % (self.nix-1, self.niy-1)) # number of cells in each dir
        for ix in range(self.nix):
            for iy in range(self.niy):
                f.write("%20.12e %20.12e\n" % (self.x[ix][iy], self.y[ix][iy]))
        return
    
    def make_grid_via_TFI(self, edge_list, cluster_list):
        """
        Given the list of 4 faces (N,E,S,W), create the grid via TFI.
        """
        print "Begin make grid, label=", self.label
        # Set up distributions of points along each of the nondimensional edges.
        rNorth = distribute_parameter_values(self.nix, cluster_list[0],
                                             edge_list[0].length() )
        sEast = distribute_parameter_values(self.niy, cluster_list[1],
                                            edge_list[1].length() )
        print "niy=", self.niy, "len(sEast)=", len(sEast)
        rSouth = distribute_parameter_values(self.nix, cluster_list[2],
                                             edge_list[2].length() )
        sWest = distribute_parameter_values(self.niy, cluster_list[3],
                                            edge_list[3].length() )
        #
        # Now work through the mesh, one point at a time,
        # blending the stretched parameter values
        # and creating the actual vertex coordinates in Cartesian space.
        #
        for iy in range(self.niy):
            s = float(iy) / (self.niy - 1)
            for ix in range(self.nix):
                r = float(ix) / (self.nix - 1)
                sdash = (1.0-r) * sWest[iy] + r * sEast[iy] 
                rdash = (1.0-s) * rSouth[ix] + s * rNorth[ix]
                self.x[ix][iy], self.y[ix][iy], zjunk = \
                                interpolate_TFI_2D_mbcns(edge_list, rdash, sdash)
            print ".",
        print 
        print "End make grid."
        return
        
#--------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin demo of geom_mbcns.py."
    #
    print "Generate a grid in a box."
    b0 = BlockGrid2D(30, 40, label="Some box")
    from gpath import *
    p0 = Node(0.0,0.0)
    p1 = Node(2.0,0.0)
    p2 = Node(2.0,3.0)
    p3 = Node(0.0,3.0)
    edge_list = [Polyline(Line(p3,p2)), Polyline(Line(p1,p2)),
                 Polyline(Line(p0,p1)), Polyline(Line(p0,p3))]
    cluster_list = [(0.01,0.02), (0,1,1.1), None, (0,1,1.1)]
    b0.make_grid_via_TFI(edge_list, cluster_list)
    #
    print "Write to VTK file."
    fout = open("test.vtk", "w")
    b0.write_block_in_VTK_format(fout)
    fout.close()
    #
    print "Done."
