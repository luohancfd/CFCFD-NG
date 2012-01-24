## \file gpath.py
## \ingroup geom
## \brief Python geometry-specification functions for mb_cns and Elmer.
##
## \author P.Jacobs
## \version 1.0 August-November 2004
## \version 1.1 31-January-2005 -- independent of Elmer and mb_cns
"""
The user-specified geometry data is organised via the following
classes: L{Node}, L{Line}, L{Arc}, L{Bezier}, L{Spline}, and L{Polyline}.
This module builds on the L{Vector} and L{Node} classes provided by
the module L{geom} and provides curvilinear path-building classes.
 
The path elements L{Line}, L{Arc}, L{Bezier} are a mix of
Python top-level classes and lower-level C functions on arrays of points.
Although it would have been much neater and more maintainable to use a pure
Python implementation, we wanted to use the same basic code for
both the C and the Python programs.

The compound L{Polyline} and L{Spline} objects are also a mix of
Python classes and C functions on GPathPolyline structures.

@undocumented: cartesian_displacement
"""
#----------------------------------------------------------------------

import sys
from copy import copy
import math
from zero_solvers import secant
from geom import *
try:
    # print "Loading gpath.py"
    import libgeom
    # print "Loaded libgeom successfully."
except:
    print "Did not load libgeom C module successfully."
from bezier_spline import bezier_3_spline

#----------------------------------------------------------------------

def cartesian_displacement(dx, dy, dz):
    """
    Returns the Cartesian displacement coordinates even if
    the first argument happens to be a Node or a Vector already.

    This is a helper function used by the translate methods below.
    These functions expect the displacement as Carestian components
    but it may be convenient for the user-scripts to supply the
    displacement as a vector (or equivalently) as a Node position
    relative to the origin.
    """
    if isinstance(dx, Vector):
        return dx.x, dx.y, dx.z
    else:
        return dx, dy, dz

#----------------------------------------------------------------------

class Line(object):
    """
    Defines a straight-line segment.

    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, a, b):
        """
        Define the directed line from L{Node} a to L{Node} b.

        @param a: Starting point on line.
        @type a: L{Node} object
        @param b: Finishing point on line.
        @type b: L{Node} object
        """
        assert isinstance(a, Node)
        assert isinstance(b, Node)
        self.np = 2
        if use_libgeom():
            self.pap = libgeom.create_point_3D_array(self.np)
            ppa = libgeom.get_point_3D_ptr(self.pap, 0)
            libgeom.point_3D_set_label(ppa, a.label)
            ppa.x = a.x; ppa.y = a.y; ppa.z = a.z
            ppb = libgeom.get_point_3D_ptr(self.pap, 1)
            libgeom.point_3D_set_label(ppb, b.label)
            ppb.x = b.x; ppb.y = b.y; ppb.z = b.z
            self.elementType = libgeom.GPATH_LINE
        else:
            self.pa = a.copy()
            self.pb = b.copy()
        return

    def copy(self, direction=1):
        """
        Make a copy of the Line, possibly reversed from the original.

        @param direction: set to 0 to get the direction of the copy reversed.
        @type direction: int
        """
        if use_libgeom():
            raise NotImplementedError, "Line.copy() not implemented in libgeom."
        if direction:
            a, b = self.pa, self.pb
        else:
            b, a = self.pa, self.pb
        return Line(a, b)
    
    def translate(self, dx, dy=0.0, dz=0.0):
        """
        Displace the Line.

        @param dx: displacement Vector representing (dx, dy, dz)
            or x-component of displacement.
        @type dx: L{Vector} or float
        @param dy: y-component of displacement (if dx was a scalar)
        @type dy: float
        @param dz: z-component of displacement (if dx was a scalar)
        @type dz: float
        """
        dx, dy, dz = cartesian_displacement(dx, dy, dz)
        if use_libgeom():
            libgeom.line_translate(self.pap, dx, dy, dz)
        else:
            self.pa.translate(dx, dy, dz)
            self.pb.translate(dx, dy, dz)
        return
        
    def eval(self, t):
        """
        Locate a point on the line.

        @param t: interpolation parameter.
        @type t: float,  0<=t<=1.0
        @returns: a L{Vector} for the point location.
        """
        loc = Vector()
        if use_libgeom():
            libgeom.line_eval(self.pap, t, loc.pp)
        else:
            loc.x = (1 - t) * self.pa.x + t * self.pb.x
            loc.y = (1 - t) * self.pa.y + t * self.pb.y
            loc.z = (1 - t) * self.pa.z + t * self.pb.z
        return loc

    def length(self):
        """
        @returns: the length of the line.
        """
        if use_libgeom():
            return libgeom.line_length(self.pap)
        else:
            return abs(self.pa - self.pb)

    def __str__(self):
        """
        String representation of the line.
        """
        if use_libgeom():
            ppa = libgeom.get_point_3D_ptr(self.pap, 0)
            ppb = libgeom.get_point_3D_ptr(self.pap, 1)
            xa = ppa.x; ya = ppa.y; za = ppa.z; labela = ppa.label
            xb = ppb.x; yb = ppb.y; zb = ppb.z; labelb = ppb.label
        else:
            xa = self.pa.x; ya = self.pa.y; za = self.pa.z; labela = self.pa.label
            xb = self.pb.x; yb = self.pb.y; zb = self.pb.z; labelb = self.pb.label
        text = "Line(" + \
               ("Node(%g,%g,%g,\"%s\")," % (xa, ya, za, labela)) + \
               ("Node(%g,%g,%g,\"%s\"))" % (xb, yb, zb, labelb))
        return text

    def __repr__(self):
        return self.__str__()


class Arc(object):
    """
    Defines a circular-arc from L{Node} a to L{Node} b with centre at c.

    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, a, b, c):
        """
        @param a: Starting point for arc.
        @type a: L{Node} object
        @param b: Finish point for arc. 
        @type b: L{Node} object
        @param c: Centre of curvature.
        @type c: L{Node} object

        @note: The radii c-->a and c-->b must match closely.
        """
        assert isinstance(a, Node), "First point should be a Node object."
        assert isinstance(b, Node), "Final point should be a Node object."
        assert isinstance(c, Node), "Centre point should be a Node object."
        self.np = 3
        dx1 = a.x - c.x
        dy1 = a.y - c.y
        dz1 = a.z - c.z
        r1 = math.sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
        dx2 = b.x - c.x
        dy2 = b.y - c.y
        dz2 = b.z - c.z
        r2 = math.sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2)
        assert abs(r1 - r2)/(r1 + 1.0) < 1.0e-5
        if use_libgeom():
            self.pap = libgeom.create_point_3D_array(self.np)
            ppa = libgeom.get_point_3D_ptr(self.pap, 0)
            libgeom.point_3D_set_label(ppa, a.label)
            ppa.x = a.x; ppa.y = a.y; ppa.z = a.z
            ppb = libgeom.get_point_3D_ptr(self.pap, 1)
            libgeom.point_3D_set_label(ppb, b.label)
            ppb.x = b.x; ppb.y = b.y; ppb.z = b.z
            ppc = libgeom.get_point_3D_ptr(self.pap, 2)
            libgeom.point_3D_set_label(ppc, c.label)
            ppc.x = c.x; ppc.y = c.y; ppc.z = c.z
            self.elementType = libgeom.GPATH_ARC
        else:
            self.pa = a.copy()
            self.pb = b.copy()
            self.pc = c.copy()
        return

    def copy(self, direction=1):
        """
        Make a copy of the Arc, possibly reversed from the original.

        @param direction: set to 0 to get the direction of the copy reversed.
        @type direction: int
        """
        if use_libgeom():
            raise NotImplementedError, "Arc.copy() not implemented in libgeom."
        if direction:
            a, b, c = self.pa, self.pb, self.pc
        else:
            b, a, c = self.pa, self.pb, self.pc
        return Arc(a, b, c)

    def translate(self, dx, dy=0.0, dz=0.0):
        """
        Displace the Arc in Cartesian space.

        @param dx: displacement Vector representing (dx, dy, dz)
            or x-component of displacement.
        @type dx: L{Vector} or float
        @param dy: y-component of displacement (if dx was a scalar)
        @type dy: float
        @param dz: z-component of displacement (if dx was a scalar)
        @type dz: float
        """
        dx, dy, dz = cartesian_displacement(dx, dy, dz)
        if use_libgeom():
            libgeom.arc_translate(self.pap, dx, dy, dz)
        else:
            self.pa.translate(dx, dy, dz)
            self.pb.translate(dx, dy, dz)
            self.pc.translate(dx, dy, dz)
        return

    def evaluate_position_and_length(self, t):
        """
        Both the position of the point and the length of the arc are evaluated
        using mostly the same process of transforming to the plane local to the arc.
        """
        loc = Vector()
        L = 0.0
        # Vector objects are not recorded in the Node.nodeList
        a = Vector(self.pa.x, self.pa.y, self.pa.z)
        b = Vector(self.pb.x, self.pb.y, self.pb.z)
        c = Vector(self.pc.x, self.pc.y, self.pc.z)
        ca = a - c; ca_mag = abs(ca)
        cb = b - c; cb_mag = abs(cb)
        # print "ca_mag=", ca_mag, "cb_mag=", cb_mag
        if abs(ca_mag - cb_mag) > 1.0e-5:
            print "Arc.eval(): radii do not match ca=", ca, " cb=", cb
            return loc, L
        # First vector in plane.
        t1 = ca.copy(); t1 = t1.unit(); # print "t1=", t1
        # Compute unit normal to plane of all three points.
        n = cross(ca, cb); # print "n=", n
        if abs(n) > 0.0:
            n = n.unit()
        else:
            print "Arc.eval(): cannot find plane of three points."
            return loc, L
        # Third (orthogonal) vector is in the original plane.
        t2 = cross(n, t1); # print "t2=", t2
        # Now transform to local coordinates so that we can do 
        # the calculation of the point along the arc in 
        # the local xy-plane, with ca along the x-axis. */
        cb_local = Vector(dot(cb,t1), dot(cb,t2), dot(cb,n))
        # print "cb_local=", cb_local
        if abs(cb_local.z) > 1.0e-6:
            print "Arc.eval(): problems with transformation cb_local=", cb_local
            return loc, L
        # Angle of the final point on the arc is in the range -pi < th <= +pi.
        theta = math.atan2(cb_local.y, cb_local.x)
        # The length of the circular arc.
        L = theta * cb_mag
        # print "whole arc: theta=", theta, "L=", L
        # Move the second point around the arc in the local xy-plane.
        theta *= t
        cb_local.x = math.cos(theta) * cb_mag
        cb_local.y = math.sin(theta) * cb_mag
        # print "in local plane: cb_local=", cb_local
        # Transform back to global xyz coordinates
        # and remember to add the centre coordinates.
        loc.x = cb_local.x * t1.x + cb_local.y * t2.x + cb_local.z * n.x + c.x
        loc.y = cb_local.x * t1.y + cb_local.y * t2.y + cb_local.z * n.y + c.y
        loc.z = cb_local.x * t1.z + cb_local.y * t2.z + cb_local.z * n.z + c.z
        # print "in global coords: loc=", loc
        return loc, L
    
    def eval(self, t):
        """
        Locate a point on the arc.

        @param t: interpolation parameter.
        @type t: float,  0<=t<=1.0
        @returns: a L{Vector} for the point location.
        """
        if use_libgeom():
            loc = Vector()
            libgeom.arc_eval(self.pap, t, loc.pp)
        else:
            loc, L = self.evaluate_position_and_length(t)
        return loc

    def length(self):
        """
        @returns: the length of the arc.
        """
        if use_libgeom():
            L = libgeom.arc_length(self.pap)
        else:
            loc, L = self.evaluate_position_and_length(1.0)
        return L

    def __str__(self):
        """
        String representation of arc.
        """
        if use_libgeom():
            ppa = libgeom.get_point_3D_ptr(self.pap, 0)
            ppb = libgeom.get_point_3D_ptr(self.pap, 1)
            ppc = libgeom.get_point_3D_ptr(self.pap, 2)
            xa = ppa.x; ya = ppa.y; za = ppa.z; labela = ppa.label
            xb = ppb.x; yb = ppb.y; zb = ppb.z; labelb = ppb.label
            xc = ppc.x; yc = ppc.y; zc = ppc.z; labelc = ppb.label
        else:
            xa = self.pa.x; ya = self.pa.y; za = self.pa.z; labela = self.pa.label
            xb = self.pb.x; yb = self.pb.y; zb = self.pb.z; labelb = self.pb.label
            xc = self.pc.x; yc = self.pc.y; zc = self.pc.z; labelc = self.pc.label
        text = "Arc(" + \
               ("Node(%g,%g,%g,\"%s\")," % (xa, ya, za, labela)) + \
               ("Node(%g,%g,%g,\"%s\")," % (xa, ya, za, labelb)) + \
               ("Node(%g,%g,%g,\"%s\"))" % (xc, yc, zc, labelc))
        return text

    def __repr__(self):
        return self.__str__()


class Arc3(Arc):
    """
    Defines a circular-arc from L{Node} a through L{Node} b  ending at L{Node} c.

    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, a, b, c):
        """
        @param a: Starting point for arc.
        @type a: L{Node} object
        @param b: Intermediate point on arc.
        @type b: L{Node} object
        @param c: Finish point for arc.
        @type c: L{Node} object

        @note: The points must not be colinear.
        """
        # 2D version of this function by Adriaan Window,
        # generalised to 3D by PJ.
        assert isinstance(a, Node), "First point should be a Node object."
        assert isinstance(b, Node), "Intermediate point should be a Node object."
        assert isinstance(c, Node), "Final point should be a Node object."
        # Compute normal vector to the plane of the circle.
        n = cross(b-a, b-c)
        assert abs(n) > 1.0e-11, "Points appear colinear."
        # The centre of the circle lies along the bisector of ab (and bc).
        ab_mid = 0.5 * (a+b)
        nab = cross(b-a, n)
        def locate_centre(s):
            return ab_mid + (s * nab)
        def error_in_radius(s):
            centre = locate_centre(s)
            return abs(a - centre) - abs(c - centre)
        s = secant(error_in_radius, 1.0, 1.01)
        centre = locate_centre(s)
        assert isinstance(centre, Vector), \
               "Couldn't locate centre,secant result= %s" % centre
        super(Arc3,self).__init__(a, c, Node(centre.x, centre.y, centre.z))
        return

    
class Bezier(object):
    """
    Defines a Bezier polynomial curve.

    @ivar N: The order of the curve is N = len(B)-1.
    @type N: int
    
    @note: The curve goes through the end-points but that
        the intermediate points generally do not lie on the curve. 
    
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, B):
        """
        Defines a Bezier polynomial of order N=len(B)-1.

        @param B: The control nodes of the curve,  B[0] --> B[-1].
        @type B: list of L{Node} objects
        """
        assert type(B) == type(list())
        self.np = len(B)
        assert self.np >= 2
        self.N = self.np - 1
        for i in range(self.np):
            assert isinstance(B[i], Node)
        if use_libgeom():
            self.pap = libgeom.create_point_3D_array(self.np)
            for i in range(self.np):
                pp = libgeom.get_point_3D_ptr(self.pap, i)
                libgeom.point_3D_set_label(pp, B[i].label)
                pp.x = B[i].x; pp.y = B[i].y; pp.z = B[i].z
            self.elementType = libgeom.GPATH_BEZIER
        else:
            self.B = []
            for node in B:
                # Note that we duplicate the Node objects completely.
                self.B.append(node.copy())
        return

    def copy(self, direction=1):
        """
        Make a copy of the Bezier, possibly reversed from the original.

        @param direction: set to 0 to get the direction of the copy reversed.
        @type direction: int
        """
        if use_libgeom():
            raise NotImplementedError, "Bezier.copy() not implemented in libgeom."
        nodeList = copy(self.B)  # copy the list of references to the Node objects
        if not direction: nodeList.reverse()
        return Bezier(nodeList)

    def translate(self, dx, dy=0.0, dz=0.0):
        """
        Displacee the Bezier curve.

        @param dx: displacement Vector representing (dx, dy, dz)
            or x-component of displacement.
        @type dx: L{Vector} or float
        @param dy: y-component of displacement (if dx was a scalar)
        @type dy: float
        @param dz: z-component of displacement (if dx was a scalar)
        @type dz: float
        """
        dx, dy, dz = cartesian_displacement(dx, dy, dz)
        if use_libgeom():
            libgeom.bezier_translate(self.N, self.pap, dx, dy, dz)
        else:
            for node in self.B:
                node.translate(dx, dy, dz)
        return

    def eval(self, t):
        """
        Locate a point on the Bezier curve.

        @param t: interpolation parameter.
        @type t: float,  0<=t<=1.0
        @returns: a L{Vector} for the point location.
        """
        loc = Vector()
        if use_libgeom():
            libgeom.bezier_eval(self.N, self.pap, t, loc.pp)
        else:
            n = len(self.B) - 1   # order of the polynomial
            if n == 3:
                # Third-order Bezier curves are expected to be very common.
                omt = 1 - t
                loc.x = omt*omt*omt*self.B[0].x + 3.0*omt*omt*t*self.B[1].x + \
                        3.0*omt*t*t*self.B[2].x + t*t*t*self.B[3].x
                loc.y = omt*omt*omt*self.B[0].y + 3.0*omt*omt*t*self.B[1].y + \
                        3.0*omt*t*t*self.B[2].y + t*t*t*self.B[3].y
                loc.z = omt*omt*omt*self.B[0].z + 3.0*omt*omt*t*self.B[1].z + \
                        3.0*omt*t*t*self.B[2].z + t*t*t*self.B[3].z
            else:
                # The general case should be less frequent.
                Qx = []; Qy = []; Qz = []
                # The following work arrays copy the original Nodes.
                # Note that we deliberately avoid creating Node objects.
                for node in self.B:
                    Qx.append(node.x)
                    Qy.append(node.y)
                    Qz.append(node.z)
                # Now, generate one new level at a time,
                # over-writing the work array at each level.
                for k in range(n):
                    for i in range(n-k):
                        Qx[i] = (1.0 - t) * Qx[i] + t * Qx[i+1]
                        Qy[i] = (1.0 - t) * Qy[i] + t * Qy[i+1]
                        Qz[i] = (1.0 - t) * Qz[i] + t * Qz[i+1]
                loc.x, loc.y, loc.z = Qx[0], Qy[0], Qz[0]
                del Qx, Qy, Qz
        return loc

    def length(self):
        """
        @returns: the length of the Bezier curve.
        @note: This is obtained approximately by sampling the curve.
        """
        if use_libgeom():
            L = libgeom.bezier_length(self.N, self.pap)
        else:
            L = 0.0
            n_seg = 10
            dt = 1.0 / n_seg
            p0 = self.eval(0.0)
            for i in range(n_seg):
                t = (i+1) * dt
                p1 = self.eval(t)
                L += abs(p1 - p0)
                p0 = p1
        return L

    def __str__(self):
        """
        String representation of the Bezier curve.
        """
        string_rep = "Bezier(["
        for i in range(self.np):
            if use_libgeom():
                pp = libgeom.get_point_3D_ptr(self.pap, i)
            else:
                pp = self.B[i]
            string_rep += ("Node(%g,%g,%g,\"%s\")," % (pp.x, pp.y, pp.z, pp.label))
        string_rep += "])"
        return string_rep

    def __repr__(self):
        return self.__str__()


class Polyline(object):
    """
    Polylines are composed of a number of gpath elements.

    This is also the data-structure used in the C-functions
    that define the edges of the grid.

    @undocumented: write_to_file, set_t0, get_t0, set_t1, get_t1,
        __slots__, plp
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """

    __slots__ = 'plp', 't0', 't1', '_elements', '_lengths', '_t0', '_t1'
    
    def __init__(self, pathElements=[], direction=1, t0=0.0, t1=1.0):
        """
        Create, possibly initialising from a list of gpath elements.

        @param pathElements: either a list of gpath objects or a single such object
        @param direction: sense in which to assemble the pathElements
        @param t0: starting value for the evaluation parameter (defines a subpath)
        @param t1: end value for the evaluation parameter (defines a subpath)
        """
        if use_libgeom():
            # The C-module has its own internal storage.
            self.plp = libgeom.gpath_init(None)
        else:
            # Places to put things for the pure-Python implementation.
            self._elements = []
            self._lengths = []
        if pathElements is None: pathElements = []
        if not isinstance(pathElements, list): pathElements = [pathElements,]
        if direction == -1: pathElements.reverse()
        for element in pathElements:
            self.append(element, direction)
        self.t0 = t0
        self.t1 = t1
        return

    def append(self, other, direction=1):
        """
        We can append simple gpath elements or another Polyline.

        By the time we have finished appending the pieces, we have
        a Polyline that consists of only simple elements.
        Appended Polylines are expanded into their simpler pieces.
        
        @param other: The item to append.
        @type other: a L{Polyline} object or other gpath object
        @param direction: A value of -1 will reverse the sense of the
            appended object.
        @type direction: int
        """
        if isinstance(other, Line) or isinstance(other, Arc) \
               or isinstance(other, Bezier):
            # A copy of a simple gpath element is appended.
            if use_libgeom():
                libgeom.gpath_add_element(self.plp, other.elementType,
                                          other.np, other.pap, direction)
            else:
                e = other.copy(direction)
                self._elements.append(e)
                self._lengths.append(e.length())
        elif isinstance(other, Polyline):
            # We expand each appended Polyline and append the simpler elements.
            if use_libgeom():
                libgeom.gpath_append_polyline(self.plp, other.plp, direction)
            else:
                eList = copy(other._elements)
                if not direction: eList.reverse()
                for e in eList:
                    self.append(e, direction)
        else:
            raise Exception, "Do not know how to append object to add to polyline."
        return

    def get_t0(self):
        if use_libgeom():
            return self.plp.t0
        else:
            return self._t0
    def set_t0(self, t0):
        if use_libgeom():
            self.plp.t0 = t0
        else:
            self._t0 = t0
        return
    t0 = property(fget=get_t0, fset=set_t0, doc="Lower bound for subrange.")

    def get_t1(self):
        if use_libgeom():
            return self.plp.t1
        else:
            return self._t1
    def set_t1(self, t1):
        if use_libgeom():
            self.plp.t1 = t1
        else:
            self._t1 = t1
        return
    t1 = property(fget=get_t1, fset=set_t1, doc="Upper bound for subrange.")

    def copy(self, direction=1):
        """
        @returns: a separate copy of the Polyline, possibly reversed.
        @param direction: Set to -1 to reverse the sense of the copy.
        @type direction: int
        """
        newPolyline = Polyline()
        newPolyline.append(self, direction)
        # Also, copy the subrange, taking care of specified direction.
        if direction == 1:
            newPolyline.t0 = self.t0; newPolyline.t1 = self.t1
        else:
            newPolyline.t0 = 1.0 - self.t1; newPolyline.t1 = 1.0 - self.t0
        return newPolyline

    def nelements(self):
        """
        @returns: the number of elements in Polyline.
        """
        if use_libgeom():
            return self.plp.ne
        else:
            return len(self._elements)

    def translate(self, dx, dy=0.0, dz=0.0):
        """
        Displace the Polyline.

        @param dx: displacement Vector representing (dx, dy, dz)
            or x-component of displacement.
        @type dx: L{Vector} or float
        @param dy: y-component of displacement (if dx was a scalar)
        @type dy: float
        @param dz: z-component of displacement (if dx was a scalar)
        @type dz: float
        """
        dx, dy, dz = cartesian_displacement(dx, dy, dz)
        if use_libgeom():
            libgeom.gpath_polyline_translate(self.plp, dx, dy, dz)
        else:
            for e in self._elements:
                e.translate(dx, dy, dz)
        return

    def length(self):
        """
        @returns: the Polyline length
        @note: The length will be updated with the addition of each new element.
        """
        if use_libgeom():
            L = self.plp.length
        else:
            L = sum(self._lengths)
        return L
    
    def eval(self, t):
        """
        Locate a point on the Polyline path.

        @param t: interpolation parameter.
        @type t: float,  0<=t<=1.0
        @returns: a L{Vector} for the point location.

        @note: A subset of the full Polyline may be selected with
               attributes t0, t1.
        """
        if use_libgeom():
            loc = Vector()
            libgeom.gpath_polyline_eval(self.plp, t, loc.pp)
        else:
            t = self.t0 + t * (self.t1 - self.t0) # map 0.0-1.0 over the subrange
            totalL = sum(self._lengths)
            L0 = 0.0
            for i in range(self.nelements()):
                e = self._elements[i]
                L = self._lengths[i]
                ta = L0 / totalL
                tb = (L0 + self._lengths[i]) / totalL
                t_local = (t - ta) / (tb - ta)
                if t_local <= 1.0:
                    # This will handle slightly negative values of t.
                    loc = e.eval(t_local)
                    break
                if i == self.nelements() - 1:
                    # We have reached the last element, force evaluation.
                    loc = e.eval(t_local)
                L0 += L
        return loc

    def __str__(self):
        text = "Polyline([\n"
        if use_libgeom():
            text += str(self.plp) + "\n"
        else:
            for e in self._elements:
                text += "    " + str(e) + ",\n"
        text += ", " + str(self.t0) + ", " + str(self.t1) + "])\n"
        return text
    
    def write_to_file(self, fp=sys.stdout):
        """
        Write all polyline elements to the specified file and
        return the number of characters written or -1 on error.
        """
        if use_libgeom():
            nc = libgeom.gpath_write_all_elements_to_file(fp, self.plp)
            if nc == -1:
                raise Exception, "Problems writing Polyline elements to file."
        else:
            fp.write(self)
            nc += len(str(text))
        return nc

    
class Spline(Polyline):
    """
    Defines a cubic-spline path.

    @undocumented: __slots__, plp
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    def __init__(self, B):
        """
        Creates the spline as a set of Bezier segments, then casts it as a Polyline.

        @param B: interpolation points, B[0] --> B[-1].
        @type B: list of L{Node} objects.

        @note: The internal representation is a set of N cubic Bezier segments
            that have the B[j] nodes as end points.

        """
        assert type(B) == type(list())
        np = len(B)
        assert np >= 2
        for i in range(np):
            assert isinstance(B[i], Node)
        bcp = bezier_3_spline(np-1, B)
        bezier_segments = []
        for i in range(np-1):
            p0 = bcp[i*4]
            p1 = bcp[i*4 + 1]
            p2 = bcp[i*4 + 2]
            p3 = bcp[i*4 + 3]
            segment = [Node(p0.x, p0.y, p0.z),
                       Node(p1.x, p1.y, p1.z),
                       Node(p2.x, p2.y, p2.z),
                       Node(p3.x, p3.y, p3.z)]
            bezier_segments.append(Bezier(segment))
        Polyline.__init__(self, bezier_segments)
        return
    
#----------------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin test of gpath.py"
    set_libgeom_flag(0)
    print "use_libgeom()=", use_libgeom()

    a = Node(-1.0, 0.0, label="a")
    b = Node(0.0, 1.0, label="b")
    c = Node(0.0, 0.0, 0.0, "c")
    d = c.copy()
    d.translate(Vector(1.0))
    d.label="something particularly long so that it gets truncated by libgeom"
    print "Example Nodes:", a, b, c, d

    ab = Line(a, b)
    print "Example Line: ab=", ab
    print "   length=", ab.length(), " midpoint=", ab.eval(0.5)
    abc = Arc(a, b, c)
    print "Example Arc: abc=", abc
    print "   length=", abc.length(), " midpoint=", abc.eval(0.5)
    abcd = Bezier( [a, b, c, d] )
    print "Example Bezier: abcd=", abcd
    print "   length=", abcd.length(), " midpoint=", abcd.eval(0.5)

    poly = Polyline([ab, abc, abcd])
    print "Example Polyline: poly=", poly

    print "------------------------------------"

    print "Set up a spline approximation to cos(t)"
    pointList = []
    N = 20
    for i in range(N+1):
        tpi2 = 2.0 * math.pi * float(i) / N
        pointList.append(Node(tpi2, math.cos(tpi2)))
    spline1 = Spline(pointList)
    # spline1.write_to_file()
    print "spline1=", spline1
    N = 10
    for i in range(N+1):
        t = float(i)/N
        tpi2 = t * math.pi * 2.0
        p = spline1.eval(t)
        print "spline1(", t, ")=", p, "cos(", p.x, ")=", math.cos(p.x)

    spline1.t0 = 0.8
    spline1.t1 = 1.0
    t = 0.5
    print "in subrange 0.8--1.0, spline1(", t, ")=", spline1.eval(t)

    print "----------------------------"
    
    print "Try translations in cartesian space."
    ab.translate(0.1, 0.2, 0.3)
    print "ab translated by increment", ab
    abc.translate(d)
    print "abc translated by d", abc
    abcd.translate(d)
    print "abcd translated by d", abcd
    dx, dy, dz = 0.0, 1.0, 1.5
    spline1.translate(dx, dy, dz)
    print "spline1 translated by (", dx, dy, dz, ")", spline1
    N = 5
    for i in range(N+1):
        t = float(i)/N
        tpi2 = t * math.pi * 2.0
        p = spline1.eval(t)
        print "spline1(", t, ")=", p, "cos(", p.x, ")+",dy,"=", math.cos(p.x)+dy

    print "----------------------------"

    print "Done."
