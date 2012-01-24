## \file geom_elmer.py
## \ingroup geom
## \brief Python geometry-specification functions for Elmer.
##
## \author P.Jacobs
## \version 20-Nov-2005 extracted from gpath.py
"""
Edge and surface definition functions for Elmer.

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
from geom import *
from gpath import *

#----------------------------------------------------------------------

class Edge3D(Polyline):
    """
    An Edge3D is a specialized Polyline that contains
    some extra data for mesh generation in the 3D flow simulation code.

    @undocumented: __slots__, write_to_file, cluster_tuple, label,
        path, t0, t1
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """

    # Specifying the allowable slots may help indicate typographical errors
    # in the user's script.
    __slots__ = 'path', 'cluster_tuple', 't0', 't1', 'label'

    def __init__(self, path=None, direction=1, t0=0.0, t1=1.0,
                 cluster_tuple=(0,0,0.0), label=""):
        """
        Initialises an edge consisting of a path, and cluster data
        that may be used for mesh generation

        @param path: may be a single path element or a list of path elements.
            The possible path elements include Line, Arc, Bezier, Polyline and
            Spline objects.
        @param direction: sense in which the path elements are assembled
        @param cluster_tuple: clustering information consisting of
            (to-end-0, to-end-1, beta)
            See roberts.py and roberts.c (distribute_points_1) for
            an explanation of the parameters.
        @param label: optional label (string) for the object 
        """
        # The path may already be a polyline with specified subrange,
        # so test for this case.
        if isinstance(path, Polyline):
            Polyline.__init__(self, [], direction, t0, t1)
            self.append(path, direction)
            # Also, copy the subrange, taking care of specified direction.
            if direction == 1:
                self.t0 = path.t0; self.t1 = path.t1
            else:
                self.t0 = 1.0 - path.t1; self.t1 = 1.0 - path.t0
        else:
            Polyline.__init__(self, path, direction, t0, t1)
        self.cluster_tuple = cluster_tuple
        self.label = label
        return

    def copy(self, direction=1):
        """
        @returns: a separate copy of the Edge3D object, possibly reversed.
        @param direction: Set to -1 to reverse the sense of the path
            for this copy.
        @type direction: int
        """
        newEdge = Edge3D()
        newEdge.append(self, direction)
        # Also, copy the subrange, taking care of specified direction.
        if direction == 1:
            newEdge.t0 = self.t0; newEdge.t1 = self.t1
        else:
            newEdge.t0 = 1.0 - self.t1; newEdge.t1 = 1.0 - self.t0
        newEdge.cluster_tuple = self.cluster_tuple
        newEdge.label = self.label
        return newEdge
    
    def write_to_file(self, filep=sys.stdout):
        """
        Writes the polyline information to the specified file.

        This behaves something like the write_polyline_data() in s_write.inc.
        """
        filep.write("%d %d %e %e %e # (label=%s) end1, end2, beta, t0, t1\n" %
                    (self.cluster_tuple[0], self.cluster_tuple[1],
                     self.cluster_tuple[2], self.t0, self.t1, self.label) )
        Polyline.write_to_file(self, filep)
        return

#----------------------------------------------------------------------

class ParametricSurface(object):
    """
    Abstract class for the surfaces that are parameterised by two parameters u,v.

    At a bare minimum, we expect to be able to call the eval(u,v) method to get
    a point on the surface.
    """
    pass

#----------------------------------------------------------------------

class CoonsPatchSurface(ParametricSurface):
    """
    A CoonsPatchSurface is defined by 4 bounding Polyline paths and
    blended linear interpolation between the edges.
    
    It can be used for the generation of wire-frame 3D blocks by
    sweeping out volumes.

    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, cA, cB, cC, cD):
        """
        Create a ClosedSurfacePatch.

        @param cA: 'South' curve
        @type cA: L{Polyline}-derived object
        @param cB: 'North' curve
        @type cB: L{Polyline}-derived object
        @param cC: 'West' curve
        @type cC: L{Polyline}-derived object
        @param cD: 'East' curve
        @type cD: L{Polyline}-derived object
        
        @note: The logical layout for the bounding curves is::
    
            .            cB
            .    +------->--------+
            .    |                |
            . cC ^                ^ cD
            .    |                |
            .    +------->--------+
            .            cA
    
        These curves must meet at the corners.

        @note: In the 3D simulation code, this CoonsPatchSurface can
            represent any of the 6 faces of a block.
            Curves cA and cB progress in the positive i-index direction
            for Top, Bottom, North and South faces.
            Curves cC and cD progress in the positive j-index direction
            for Top and Bottom faces.
            Curves cA and cB progress in the positive j-index direction
            for West and East faces.
            Curves cC and cD progress in the positive k-index direction
            for North, South, West and East faces. 
        """
        assert isinstance(cA, Polyline), "cA needs to be a Polyline or an Edge3D."
        assert isinstance(cB, Polyline), "cB needs to be a Polyline or an Edge3D."
        assert isinstance(cC, Polyline), "cC needs to be a Polyline or an Edge3D."
        assert isinstance(cD, Polyline), "cD needs to be a Polyline or an Edge3D."
        tolerance = 1.0e-6
        distSW = distance_between_nodes(cA.eval(0.0), cC.eval(0.0))
        assert distSW < tolerance, "CoonsPatchSurface: south-west corner gap=%g" % distSW
        distSE = distance_between_nodes(cA.eval(1.0), cD.eval(0.0))
        assert distSE < tolerance, "CoonsPatchSurface: south-east corner gap=%g" % distSE
        distNW = distance_between_nodes(cB.eval(0.0), cC.eval(1.0))
        assert distNW < tolerance, "CoonsPatchSurface: north-west corner gap=%g" % distNW
        distNE = distance_between_nodes(cB.eval(1.0), cD.eval(1.0))
        assert  distNE < tolerance, "CoonsPatchSurface: north-east corner gap=%g" % distNE
        self.cA = cA.copy()
        self.cB = cB.copy()
        self.cC = cC.copy()
        self.cD = cD.copy()
        # Corners are needed for TFI interpolation.
        self.p00 = self.cA.eval(0.0)
        self.p10 = self.cA.eval(1.0)
        self.p01 = self.cB.eval(0.0)
        self.p11 = self.cB.eval(1.0)
        return

    def eval(self, r, s):
        """
        Locate a point on the CoonsPatchSurface by blended linear interpolation.

        @param r: interpolation parameter for along curves cA and cB
        @type r: float, 0.0<=r<=1.0
        @param s: interpolation parameter for along curves cC and cD
        @type s: float, 0.0<=s<=1.0

        @returns: a L{Vector} value for the point.
        """
        p = Vector()
        if use_libgeom():
            flag = libgeom.coons_patch(self.cA.plp, self.cB.plp,
                                       self.cC.plp, self.cD.plp,
                                       r, s, p.pp)
        else:
            cAr = self.cA.eval(r); cBr = self.cB.eval(r)
            cCs = self.cC.eval(s); cDs = self.cD.eval(s)
            # Although the following is long-winded, am trying to reduce the
            # number of created objects.
            x = (1-s) * cAr.x + s * cBr.x + (1-r) * cCs.x + r * cDs.x - \
                ( (1-r)*(1-s)*self.p00.x + (1-r)*s*self.p01.x +
                  r*(1-s)*self.p10.x + r*s*self.p11.x )
            y = (1-s) * cAr.y + s * cBr.y + (1-r) * cCs.y + r * cDs.y - \
                ( (1-r)*(1-s)*self.p00.y + (1-r)*s*self.p01.y +
                  r*(1-s)*self.p10.y + r*s*self.p11.y )
            z = (1-s) * cAr.z + s * cBr.z + (1-r) * cCs.z + r * cDs.z - \
                ( (1-r)*(1-s)*self.p00.z + (1-r)*s*self.p01.z +
                  r*(1-s)*self.p10.z + r*s*self.p11.z )
            p = Vector(x, y, z)
        return p 

    def extrude(self, cE, direction):
        """
        Extrudes the CoonsPatchSurface to form a wire-frame (closed) volume.

        @param cE: curve along which the extrusion is done.
        @type cE: L{Polyline}-derived object
        @param direction: provides a hint as to which way
           we want to extrude the surface.
        @type direction: string being one of 'i', 'j', or 'k'

        @returns: the list of 12 edges defining a 3D wire-frame block.
        """
        assert isinstance(cE, Polyline), "Extrusion curve needs to be a Polyline or Edge3D."
        if direction == "k":
            d0 = cE.eval(0) - self.cA.eval(0)  # starting offset
            c01 = self.cA.copy()   # start assembling bottom surface
            c01.translate(d0) # move bottom surface to start of extrusion line
            c12 = self.cD.copy()
            c12.translate(d0)
            c32 = self.cB.copy()
            c32.translate(d0)
            c03 = self.cC.copy()
            c03.translate(d0)
            #
            d1 = cE.eval(1) - self.cA.eval(0)  # final offset
            c45 = self.cA.copy()   # start assembling top surface
            c45.translate(d1) # move bottom surface to start of extrusion line
            c56 = self.cD.copy()
            c56.translate(d1)
            c76 = self.cB.copy()
            c76.translate(d1)
            c47 = self.cC.copy()
            c47.translate(d1)
            #
            c04 = cE.copy()
            c15 = cE.copy()
            c15.translate(c01.eval(1) - cE.eval(0))
            c26 = cE.copy()
            c26.translate(c12.eval(1) - cE.eval(0))
            c37 = cE.copy()
            c37.translate(c03.eval(1) - cE.eval(0))
        elif direction == "i":
            d0 = cE.eval(0) - self.cA.eval(0)  # starting offset
            c03 = self.cA.copy()   # start assembling bottom surface
            c03.translate(d0) # move bottom surface to start of extrusion line
            c37 = self.cD.copy()
            c37.translate(d0)
            c47 = self.cB.copy()
            c47.translate(d0)
            c04 = self.cC.copy()
            c04.translate(d0)
            #
            d1 = cE.eval(1) - self.cA.eval(0)  # final offset
            c12 = self.cA.copy()   # start assembling top surface
            c12.translate(d1) # move bottom surface to start of extrusion line
            c26 = self.cD.copy()
            c26.translate(d1)
            c56 = self.cB.copy()
            c56.translate(d1)
            c15 = self.cC.copy()
            c15.translate(d1)
            #
            c01 = cE.copy()
            c32 = cE.copy()
            c32.translate(self.cA.eval(1) - cE.eval(0)) ## FIX ME
            c76 = cE.copy()
            c76.translate(self.cD.eval(1) - cE.eval(0)) ## FIX ME
            c45 = cE.copy()
            c45.translate(self.cC.eval(1) - cE.eval(0)) ## FIX ME
        elif direction == "j":
            d0 = cE.eval(0) - self.cA.eval(0)  # starting offset
            c01 = self.cA.copy()   # start assembling bottom surface
            c01.translate(d0) # move bottom surface to start of extrusion line
            c15 = self.cD.copy()
            c15.translate(d0)
            c45 = self.cB.copy()
            c45.translate(d0)
            c04 = self.cC.copy()
            c04.translate(d0)
            #
            d1 = cE.eval(1) - self.cA.eval(0)  # final offset
            c32 = self.cA.copy()   # start assembling top surface
            c32.translate(d1) # move bottom surface to start of extrusion line
            c26 = self.cD.copy()
            c26.translate(d1)
            c76 = self.cB.copy()
            c76.translate(d1)
            c37 = self.cC.copy()
            c37.translate(d1)
            #
            c03 = cE.copy()
            c12 = cE.copy()
            c12.translate(self.cA.eval(1) - cE.eval(0)) ## FIX ME
            c56 = cE.copy()
            c56.translate(self.cD.eval(1) - cE.eval(0)) ## FIX ME
            c47 = cE.copy()
            c47.translate(self.cC.eval(1) - cE.eval(0)) ## FIX ME
        else:
            print "Invalid direction for extrusion:", direction
        return [c01, c12, c32, c03, c45, c56, c76, c47, c04, c15, c26, c37]

ClosedSurfacePatch = CoonsPatchSurface # for compatibility with old scripts.

#---------------------------------------------------------------------

class ParametricVolume(object):
    """
    A volume is defined by its surfaces and/or its edges.
    """
    def __init__(self, surface_list=[None,]*6, edge_list=[None,]*12):
        """
        @param edge_list: list of 12 paths in Elmer order.
        @type edge_list: list of L{Polyline}-derived objects.
        @param surface_list: list of 6 paths in Elmer order.
        @type surface_list: list of L{ParametricSurface}-derived objects.
        """
        assert isinstance(edge_list, list), "edge_list should be a list"
        assert len(edge_list) == 12, "edge_list should have 12 edges."
        for edge in edge_list:
            assert isinstance(edge, Polyline), "edge should be at least a Polyline."
        assert isinstance(surface_list, list), "surface_list should be a list"
        assert len(surface_list) == 6, "surface_list should have 6 surfaces."
        for f in surface_list:
            if f != None:
                assert isinstance(f, ParametricSurface), \
                       "surface should be at least a ParametricSurface."
        # The pure-Python module makes use of the surfaces in the first instance.
        fSouth = surface_list[0]; fBottom = surface_list[1]
        fWest = surface_list[2]; fEast = surface_list[3]
        fNorth = surface_list[4]; fTop = surface_list[5]
        # If any particular surface is missing, assemble it as a CoonsPatchSurface.
        # See page 42 of 2005 workbook for notation.
        c = edge_list
        if fSouth == None:
            fSouth = CoonsPatchSurface(c[0], c[4], c[8], c[9])
        if fBottom == None:
            fBottom = CoonsPatchSurface(c[0], c[2], c[3], c[1])
        if fWest == None:
            fWest = CoonsPatchSurface(c[3], c[7], c[8], c[11])
        if fEast == None:
            fEast = CoonsPatchSurface(c[1], c[5], c[9], c[10])
        if fNorth == None:
            fNorth = CoonsPatchSurface(c[2], c[6], c[11], c[10])
        if fTop == None:
            fTop = CoonsPatchSurface(c[4], c[6], c[7], c[5])
        # Should check for consistency between edges and surfaces.
        self.edge_list = edge_list
        self.surface_list = [fSouth, fBottom, fWest, fEast, fNorth, fTop]
        # The corners are needed for TFI interpolation.
        # Maybe we should do the following with surface evaluations.
        self.p000 = c[0].eval(0.0); self.p100 = c[0].eval(1.0)
        self.p010 = c[2].eval(0.0); self.p110 = c[2].eval(1.0)
        self.p001 = c[4].eval(0.0); self.p101 = c[4].eval(1.0)
        self.p011 = c[6].eval(0.0); self.p111 = c[6].eval(1.0)
        return
    
    def eval(self, r, s, t):
        """
        Locate a point p(r,s,t) using transfinite interpolation.

        @param r: interpolation parameter in the i-index direction
        @type r: float, 0.0<=r<=1.0
        @param s: interpolation parameter in the j-index direction
        @type s: float, 0.0<=s<=1.0
        @param t: interpolation parameter in the k-index direction
        @type t: float, 0.0<=t<=1.0
        """
        if use_libgeom():
            # The C-module libgeom knows only the wire-frame type of volume.
            p = Vector()
            c = self.edge_list
            flag = libgeom.TFI_3D_no_array(c[0].plp, c[1].plp, c[2].plp, c[3].plp,
                                           c[4].plp, c[5].plp, c[6].plp, c[7].plp,
                                           c[8].plp, c[9].plp, c[10].plp, c[11].plp,
                                           r, s, t, p.pp)
        else:
            # The pure-Python module makes use of the surfaces in the first instance.
            fSouth = self.surface_list[0]; fBottom = self.surface_list[1]
            fWest = self.surface_list[2]; fEast = self.surface_list[3]
            fNorth = self.surface_list[4]; fTop = self.surface_list[5]
            # Although the following is long-winded, we are trying to reduce the
            # number of temporary objects within the Vector functions.
            fW = fWest.eval(s,t); fE = fEast.eval(s,t)
            fS = fSouth.eval(r,t); fN = fNorth.eval(r,t)
            fB = fBottom.eval(r,s); fT = fTop.eval(r,s)
            omr = 1.0 - r; oms = 1.0 - s; omt = 1.0 - t
            BigCx = omr*oms*omt*self.p000.x + omr*oms*t*self.p001.x + \
                    omr*s*omt*self.p010.x   + omr*s*t*self.p011.x + \
                    r*oms*omt*self.p100.x   + r*oms*t*self.p101.x + \
                    r*s*omt*self.p110.x     + r*s*t*self.p111.x
            x = 0.5 * ( omr * fW.x + r * fE.x +
                        oms * fS.x + s * fN.x +
                        omt * fB.x + t * fT.x ) - 0.5 * BigCx
            BigCy = omr*oms*omt*self.p000.y + omr*oms*t*self.p001.y + \
                    omr*s*omt*self.p010.y   + omr*s*t*self.p011.y + \
                    r*oms*omt*self.p100.y   + r*oms*t*self.p101.y + \
                    r*s*omt*self.p110.y     + r*s*t*self.p111.y
            y = 0.5 * ( omr * fW.y + r * fE.y +
                        oms * fS.y + s * fN.y +
                        omt * fB.y + t * fT.y ) - 0.5 * BigCy
            BigCz = omr*oms*omt*self.p000.z + omr*oms*t*self.p001.z + \
                    omr*s*omt*self.p010.z   + omr*s*t*self.p011.z + \
                    r*oms*omt*self.p100.z   + r*oms*t*self.p101.z + \
                    r*s*omt*self.p110.z     + r*s*t*self.p111.z
            z = 0.5 * ( omr * fW.z + r * fE.z +
                        oms * fS.z + s * fN.z +
                        omt * fB.z + t * fT.z ) - 0.5 * BigCz
            p = Vector(x, y, z)
        return p

#----------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin test of geom_elmer.py"
    set_libgeom_flag(0)
    print "use_libgeom()=", use_libgeom()

    print "Try transfinite interpolation inside a box."
    from simple_boxes import *
    v1 = ParametricVolume(edge_list=simpleBox())
    print v1.eval(0.5, 0.25, 0.3)

    print "----------------------------"
    print "Try a CoonsPatchSurface and Extrusion."
    n0 = Node(0.0, 0.0)
    n1 = Node(1.0, 0.0)
    n2 = Node(1.0, 1.0)
    n3 = Node(0.0, 1.0)
    cs1 = CoonsPatchSurface(Edge3D([Line(n0,n1),]), # south
                            Edge3D([Line(n3,n2),]), # north
                            Edge3D([Line(n0,n3),]), # west
                            Edge3D([Line(n1,n2),])) # east
    print "p(0.5,0.25)=", cs1.eval(0.5, 0.25)

    n4 = Node(0.0, 0.0, 1.0)
    r = 0.5; s = 0.25; t = 0.3
    for direction in ["k", "i", "j"]:
        print "Extrude in ", direction, "-direction"
        v2 = ParametricVolume(edge_list=cs1.extrude(Edge3D([Line(n0,n4),]),direction))
        print "p(", r, s, t, ")=", v2.eval(r, s, t)
    
    print "Done."
