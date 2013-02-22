/** \file libgeom2.i
 *  \ingroup libgeom2 
 *  \brief SWIG interface header file for C++ geometry library.
 *  \author PJ
 *  \version 29-Dec-2005 initial coding
 *  \version 31-Dec-2005 eliminated label from Vector3 but kept it in Node3.
 *  \version 13-Jan-2006 up to volume classes.
 *  \version 17-Jan-2006 Node class moved from C++ into Python code.
 *  \version 21-Feb-2013 Exception handler added.
 */
%define GEOM2_DOCSTRING
"Python interface to the geometry and path-defining classes, now in implemented C++"
%enddef
%module(docstring=GEOM2_DOCSTRING) libgeom2

%include "exception.i"
%exception {
    try {
        $action
    }
    // Could put custom exception catches here.
    // catch ( const std::exception & e ) {
    //    std::cout << "hello from my exception handler" << endl;
    //    SWIG_exception(SWIG_RuntimeError, e.what());
    // }
    // but, for the moment, we're happy with the standard set.
    SWIG_CATCH_STDEXCEPT
}
%include "std_string.i"
%include "std_vector.i"

%{
#include <string>
#include <vector>
#include <stdio.h>
#include "geom.hh"
#include "../../nm/source/fobject.hh"
#include "gpath.hh"
#include "pypath.hh"
#include "surface.hh"
#include "volume.hh"
#include "pysurface.hh"
#include "pyvolume.hh"
%}

#ifdef SWIGPYTHON
%typemap(in) FILE * {
    if (!PyFile_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "expected PyFile");
        return NULL;
    }
    $1=PyFile_AsFile($input);
}
#endif

%rename(ostream_print) operator<<( ostream &os, const Vector3 &v );
%rename(Vector3_unary_plus) operator+( const Vector3 &v );
%rename(Vector3_unary_minus) operator-( const Vector3 &v );
%rename(Vector3_add_vv) operator+( const Vector3 &v1, const Vector3 &v2 );
%rename(Vector3_sub_vv) operator-( const Vector3 &v1, const Vector3 &v2 );
%rename(Vector3_mul_vd) operator*( const Vector3 &v1, double v2 );
%rename(Vector3_mul_dv) operator*( double v1, const Vector3 &v2 );
%rename(Vector3_div_vr) operator/( const Vector3 &v1, double v2 );

%ignore Vector3::operator=;
%include "geom.hh"

%extend Vector3 {
    char *__str__() {
        static char tmp[1024];
	strncpy(tmp, self->str().c_str(), 1023);
	tmp[1023] = '\0';
	return tmp;
    }
};

%pythoncode %{
# Connect these renamed operators to the appropriate special methods.
# We seem to have to define Python functions that call the wrapped
# C++ functions and then set the special attributes to these Python functions.

def Vector3_iadd_v(self, other):
    # Yes, this seems wasteful but I cannot get access the C++ += operator.
    new = self + Vector3(other)
    return new
Vector3.__iadd__ = Vector3_iadd_v

def Vector3_isub_v(self, other):
    new = self - Vector3(other)
    return new
Vector3.__isub__ = Vector3_isub_v

def Vector3_pos(self):
    return Vector3_unary_plus(self)
Vector3.__pos__ = Vector3_pos

def Vector3_neg(self):
    return Vector3_unary_minus(self)
Vector3.__neg__ = Vector3_neg

def Vector3_add(self, other):
    return Vector3_add_vv(self, Vector3(other))
Vector3.__add__ = Vector3_add

def Vector3_sub(self, other):
    return Vector3_sub_vv(self, Vector3(other))
Vector3.__sub__ = Vector3_sub

def Vector3_mul(self, other):
    return Vector3_mul_vd(self, other)
Vector3.__mul__ = Vector3_mul

def Vector3_rmul(self,other):
    return Vector3_mul_dv(other, self)
Vector3.__rmul__ = Vector3_rmul

def Vector3_div(self, other):
    return Vector3_div_vr(self, other)
Vector3.__div__ = Vector3_div

# Want the Vector3 class to be also known as Vector in Python.
Vector = Vector3
%}


%pythoncode %{
# Some Red-Green-Blue colour tuples.
vrml_pure_red   = (1.0, 0.0, 0.0)
vrml_pure_green = (0.0, 1.0, 0.0)
vrml_pure_blue  = (0.0, 0.0, 1.0)
vrml_white      = (1.0, 1.0, 1.0)
vrml_black      = (0.0, 0.0, 0.0)
vrml_yellow     = (1.0, 1.0, 0.0)
vrml_cyan       = (0.0, 1.0, 1.0)
vrml_magenta    = (1.0, 0.0, 1.0)
vrml_light_gray = (0.75, 0.75, 0.75)
vrml_medium_gray= (0.5, 0.5, 0.5)
vrml_dark_gray  = (0.25, 0.25, 0.25)
vrml_dark_red   = (0.5, 0.0, 0.0)
vrml_dark_green = (0.0, 0.5, 0.0)
vrml_dark_blue  = (0.0, 0.0, 0.5)

class Node(Vector3):
    """
    Defines a (labelled) nodal-point in 3D space

    Often these Nodes will be subsequently used in the definition
    of Path segments.

    @undocumented: label, nodeList
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    nodeList = []  # a collection so that we can later render the nodes
                   # to MetaPost or whatever.
    
    def __init__(self, x=0.0, y=0.0, z=0.0, label=""):
        # Usually, we expect to be given Cartesian coordinates.
        if isinstance(x, Vector3):
            # We may have been handed a Vector3 object.
            x, y, z = x.x, x.y, x.z
        Vector3.__init__(self, x, y, z)
        self.label = label
        Node.nodeList.append(self)
        return

    def __str__(self):
        return "Node(%g,%g,%g, \"%s\")" % (self.x, self.y, self.z, self.label)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        """
        Create and return a new instance which is a copy of the original.
        """
        return Node(self.x, self.y, self.z, self.label)

    def clone(self):
        """
        Create and return a new instance which is a copy of the original.

	This function matches the name more commonly used in the C++ code.
        """
        return Node(self.x, self.y, self.z, self.label)
    
    def translate(self, dx, dy=0.0, dz=0.0):
        """
        Translate node position by displacement Vector dx
        or by displacement (dx, dy, dz) in Cartesian coordinates.

        @param dx: displacement
        @type dx: either L{Vector} or float
        @param dy: (optional) Cartesian displacement in the y-direction
        @type dy: float
        @param dz: (optional) Cartesian displacement in the z-direction
        @type dz: float
        """
        if isinstance(dx, Vector3):
            dx, dy, dz = dx.x, dx.y, dx.z
        self.x += dx; self.y += dy; self.z += dz
        return self

    def vrml_str(self, rgb=vrml_black, radius=0.001, offset=0.001, font_size=0.005 ):
        ost = ["\nTransform {",]
        ost.append("  # Node3 label= \"%s\"" % self.label)
        ost.append("  translation %f %f %f" % (self.x, self.y, self.z))
        ost.append("  children [")
        ost.append("    Shape {")
        ost.append("      appearance Appearance {")
        ost.append("        material Material { diffuseColor %g %g %g }" % rgb)
        ost.append("      }")
        ost.append("      geometry Sphere { radius %f }" % radius)
        ost.append("    }")
        ost.append("    Transform {")
        ost.append("      # label")
        ost.append("      translation %f %f 0.0" % (offset, offset))
        ost.append("      children [")
        ost.append("        Shape {")
        ost.append("          appearance Appearance {")
        ost.append("            material Material { diffuseColor %g %g %g }" % rgb)
        ost.append("          }")
        ost.append("          geometry Text { ")
        ost.append("            string \"%s\"" % self.label)
        ost.append("            fontStyle FontStyle {")
        ost.append("              size %g" % font_size )
        ost.append("              family \"SANS\"")
        ost.append("              justify \"BEGIN\"")
        ost.append("            }")
        ost.append("          }")
        ost.append("        }")
        ost.append("      ]")
        ost.append("    }")
        ost.append("  ]")
        ost.append("}\n")
        return "\n".join(ost)
%}

// The following magic allows us to pass a list of Vector3 objects
// to be collected as a C++ vector in the Bezier constructor.
%template(vectorVector3) std::vector<Vector3>;
%template(vectorVector3Ptr) std::vector<Vector3*>;
%template(vectorPathPtr) std::vector<Path*>;
%template(vectorXPathPtr) std::vector<XPath*>;
#ifndef STD_VECTOR_TEMPLATES_ALREADY_DEFINED
%template(vectori) std::vector<int>;
%template(vectord) std::vector<double>;
#define STD_VECTOR_TEMPLATES_ALREADY_DEFINED
#endif


%rename(ostream_print_UnivariateFunction) operator<<( ostream &os, const UnivariateFunction &f );
%rename(ostream_print_BivariateFunction) operator<<( ostream &os, const BivariateFunction &f );
%include"../../nm/source/fobject.hh"


%rename(ostream_print_Path) operator<<( ostream &os, const Path &v );
%ignore XPoly::operator=;
%include "gpath.hh"
%include "pypath.hh"

%extend Path {
    char *__str__() {
        static char tmp[4096];
	strncpy(tmp, self->str().c_str(), 4095);
	tmp[4095] = '\0';
	return tmp;
    }
    char *get_label() {
        static char tmp[128];
	strncpy(tmp, self->label.c_str(), 127);
	tmp[127] = '\0';
	return tmp;
    }
    void set_label( const string new_label ) {
	self->label = new_label;
	return;
    }
};

%pythoncode %{
def Path_get_name(self):
    return self.get_name()
def Path_set_name(self, name):
    self.set_name(name)
    return
Path.name = property(fget=Path_get_name, fset=Path_set_name)

def Path_vrml_str(self, rgb=vrml_black, n=10):
    """
    Returns the string representation in VRML format.
    """
    ost = "\nShape {"
    ost += "\n    # Path label= \"" + self.label + "\""
    ost += "\n    appearance Appearance {"
    ost += "\n        material Material { diffuseColor %g %g %g }" % rgb
    ost += "\n    }"
    dt = 1.0/(n-1)
    ost += "\n    geometry IndexedLineSet {"
    ost += "\n        coord Coordinate {"
    ost += "\n            point ["
    for i in range(n):
        p = self.eval(i*dt)
       	ost += "\n                %g %g %g," % (p.x, p.y, p.z)
    p = self.eval(1.0) 
    ost += "\n                %g %g %g" % (p.x, p.y, p.z)
    ost += "\n            ]"
    ost += "\n        }"
    ost += "\n        coordIndex ["
    ost += "\n            "
    for i in range(n):
	ost += "%d, " % i
    ost += "-1"
    ost += "\n        ]"
    ost += "\n    } # end indexedLineSet"
    ost += "\n} # end Shape";
    return ost
Path.vrml_str = Path_vrml_str


def Polyline2( *args ):
    """
    Constructs a Polyline from a arbitrary number of arguments
    that may be Vectors (points) or Path elements.

    Example usage:
    >>> from libgeom2 import *
    >>> a = Polyline2(Vector(1.0,2.0),
    ...               Line(Vector(1.0,3.0),Vector(1.0,4.0)),
    ...               Vector(1.0,5.0))
    >>> print a.eval(0.5)
    Vector3(1, 3.5, 0)
    """
    from cfpylib.util.flatten import flatten
    seq = flatten(args)
    # print "received args=", seq
    segment_list = []
    prev_item = None
    for item in seq:
        if prev_item is None:
            # print "First item", item
            if isinstance(item, Path): segment_list.append(item)
        elif isinstance(prev_item, Vector3) and isinstance(item, Vector3):
            # print "Point to point", prev_item, item
            segment_list.append(Line(prev_item, item))
        elif isinstance(prev_item, Path) and isinstance(item, Vector3):
            # print "Path to point", prev_item, item
            segment_list.append(Line(prev_item.eval(1.0), item))
        elif isinstance(prev_item, Vector3) and isinstance(item, Path):
            # print "Point to path", prev_item, item
            segment_list.append(Line(prev_item, item.eval(0.0)))
            segment_list.append(item)
        elif isinstance(prev_item, Path) and isinstance(item, Path):
            # print "Path to path", prev_item, item
            old_end = prev_item.eval(1.0)
            new_start = item.eval(0.0)
            if vabs(old_end - new_start) > 1.0e-7:
                # There is a gap that needs to be filled.
                segment_list.append(Line(old_end, new_start))
            segment_list.append(item)
        else:
            print "Polyline2(): unknown item(s):", prev_item, item
        # Remember the current item for the next iteration.
        prev_item = item
    # print "Final segment list:"
    # for item in segment_list:
    #     print item, "start", item.eval(0.0), "end", item.eval(1.0)
    return Polyline(segment_list)

def Spline2( fname ):
    """
    Contructs a spline from a file containing x(,y(,z)) coordinates.

    This function takes a filename and processes it assuming that each
    line contatins (x,y,z) triples (space-delimited).  If any values are
    missing on a given line, they are assumed to be 0.0.  The x,y,z-triples
    are gathered and used to create a Spline. This Spline is returned
    to the caller.

    Inputs:
    -------
    fname = file name, text file with x,y,z triples (one per line)

    Returns:
    --------
    Spline object
    """

    fp = open(fname, "r")
    points = []
    for line in fp.readlines():
        tks = line.split()
        if len(tks) == 0:
            continue
        x = float(tks[0])
        try:
            y = float(tks[1])
        except:
            y = 0.0
        try:
            z = float(tks[2])
        except:
            z = 0.0
        points.append(Vector3(x, y, z))

    return Spline(points)
%}


%template(vectorParametricSurfacePtr) std::vector<ParametricSurface*>;
%rename(ostream_print_ParametricSurface) operator<<( ostream &os, const ParametricSurface &surf );
%include "surface.hh"
%include "pysurface.hh"

%extend ParametricSurface {
    char *__str__() {
        // Just in case the string gets very long, 
        // this function should safely pass the first part of it 
        // through the SWIG wrapper.
        static char tmp[4096];
	strncpy(tmp, self->str().c_str(), 4095);
	tmp[4095] = '\0';
	return tmp;
    }
    char *get_label() {
        static char tmp[128];
	strncpy(tmp, self->label.c_str(), 127);
	tmp[127] = '\0';
	return tmp;
    }
    void set_label( const string new_label ) {
	self->label = new_label;
	return;
    }
};

%pythoncode %{
def ParametricSurface_get_name(self):
    return self.get_name()
def ParametricSurface_set_name(self, name):
    self.set_name(name)
    return
ParametricSurface.name = property(fget=ParametricSurface_get_name, 
                                  fset=ParametricSurface_set_name)

def ParametricSurface_vrml_str(self, rgb=vrml_black, ni=10, nj=10, draw_as_mesh=1):
    """
    Returns a string representing the surface in VRML format.
    """
    ost = "# ParametricSurface: " + self.label
    # The following process is pretty inefficient and causes the
    # surface to be cloned 6 times, however, we don't expect to call
    # upon it frequently.
    f_r = LinearFunction()
    f_s = LinearFunction()
    if draw_as_mesh:
        dr = 1.0/(ni-1)
	for i in range(ni):
            r = i * dr
            f = LinearFunction(0.0, r)
            pth = PathOnSurface(self, f, f_s)
            ost += "\n# Path at fixed r=%g for meshline" % r
            ost += pth.vrml_str(rgb, nj)
        ds = 1.0/(nj-1)
	for j in range(nj):
            s = j * ds
            f = LinearFunction(0.0, s)
            pth = PathOnSurface(self, f_r, f)
            ost += "\n# Path at fixed s=%g for meshline" % s
            ost += pth.vrml_str(rgb, ni)
    else:
        # Draw an boundary paths with a cross in the middle.
        f_zero = LinearFunction(0.0, 0.0)
        f_half = LinearFunction(0.0, 0.5)
        f_one = LinearFunction(0.0, 1.0)
        ost += "\n# cA boundary"
        pA = PathOnSurface(self, f_r, f_zero)
        ost += pA.vrml_str(rgb, ni)
        ost += "\n# cB boundary"
        pB = PathOnSurface(self, f_r, f_one)
        ost += pB.vrml_str(rgb, ni)
        ost += "\n# cC boundary"
        pC = PathOnSurface(self, f_zero, f_s)
        ost += pC.vrml_str(rgb, nj)
        ost += "\n# cD boundary"
        pD = PathOnSurface(self, f_one, f_s)
        ost += pD.vrml_str(rgb,nj)
        ost += "\n# Now, put draw lines through the midpoint."
        pmid1 = PathOnSurface(self, f_half, f_s, "mid1", 0.05, 0.95)
        ost +=pmid1.vrml_str(rgb, nj)
        pmid2 = PathOnSurface(self, f_r, f_half, "mid2", 0.05, 0.95)
        ost += pmid2.vrml_str(rgb, ni)
    ost += "\n# end ParametricSurface"
    return ost
ParametricSurface.vrml_str = ParametricSurface_vrml_str
%}


%rename(ostream_print_ParametricVolume) operator<<( ostream &os, const ParametricVolume &v );
%include "volume.hh"
%include "pyvolume.hh"

%extend ParametricVolume {
    char *__str__() {
        static char tmp[4096];
	strncpy(tmp, self->str().c_str(), 4095);
	tmp[4095] = '\0';
	return tmp;
    }
    char *get_label() {
        static char tmp[128];
	strncpy(tmp, self->label.c_str(), 127);
	tmp[127] = '\0';
	return tmp;
    }
    void set_label( const string new_label ) {
	self->label = new_label;
	return;
    }
};

%pythoncode %{
def ParametricVolume_get_name(self):
    return self.get_name()
def ParametricVolume_set_name(self, name):
    self.set_name(name)
    return
ParametricVolume.name = property(fget=ParametricVolume_get_name, 
                                 fset=ParametricVolume_set_name)

def ParametricVolume_vrml_str(self, rgb=vrml_black, ni=10, nj=10, nk=10, draw_as_mesh=1):
    """
    Returns the VRML string representation of the volume as the
    concatenation of the strings for the 6 bounding surfaces.
    """
    # FIX-ME, I think...
    # Need to render the subsection of the volume as VRML surfaces
    # together with the original boundary surfaces rendered as wire-frame boxes. 
    ost = "# ParametricVolume: " + self.label
    ost += "\n# (Note that this is the full volume, even if a subsection is specified.)"
    ost += "\n# south boundary"
    ost += self.south.vrml_str(rgb, ni, nk, draw_as_mesh)
    ost += "\n# bottom boundary"
    ost += self.bottom.vrml_str(rgb, ni, nj, draw_as_mesh)
    ost += "\n# west boundary"
    ost += self.west.vrml_str(rgb, nj, nk, draw_as_mesh)
    ost += "\n# east boundary"
    ost += self.east.vrml_str(rgb, nj, nk, draw_as_mesh)
    ost += "\n# north boundary"
    ost += self.north.vrml_str(rgb, ni, nk, draw_as_mesh)
    ost += "\n# top boundary"
    ost += self.top.vrml_str(rgb, ni, nj, draw_as_mesh)
    ost += "\n# end ParametricVolume"
    return ost
ParametricVolume.vrml_str = ParametricVolume_vrml_str


def makeSimpleBox(xPos=0.0, yPos=0.0, zPos=0.0, xSize=1.0, ySize=1.0, zSize=1.0):
    """
    Creates a simple box volume from a starting point and
    a size in each direction.
    """
    p0 = Node(xPos,       yPos,       zPos)
    p1 = Node(xPos+xSize, yPos,       zPos)
    p2 = Node(xPos+xSize, yPos+ySize, zPos)
    p3 = Node(xPos,       yPos+ySize, zPos)
    p4 = Node(xPos,       yPos,       zPos+zSize)
    p5 = Node(xPos+xSize, yPos,       zPos+zSize)
    p6 = Node(xPos+xSize, yPos+ySize, zPos+zSize)
    p7 = Node(xPos,       yPos+ySize, zPos+zSize)
    return SimpleBoxVolume(p0, p1, p2, p3, p4, p5, p6, p7)
%}

%pythoncode %{
def RCF(a,b,beta):
    return RobertsClusterFunction(a, b, beta)

def HCF(dL0,dL1):
    return HypertanClusterFunction(dL0,dL1)

def BRCF(a,b,c,d,beta0,beta1,gamma):
    RCF_a = RCF(a,b,beta0)
    RCF_b = RCF(c,d,beta1)
    return DiscontinuousUnivariateFunction(gamma,RCF_a,RCF_b)

def TRCF( a, b, c, d, e, f, beta0, beta1, beta2, gamma0, gamma1 ):
    rcf0 = RCF(a,b,beta0)
    rcf1 = RCF(c,d,beta1)
    rcf2 = RCF(e,f,beta2)
    brcf01 = DiscontinuousUnivariateFunction(gamma0,rcf0,rcf1)
    return DiscontinuousUnivariateFunction(gamma1,brcf01,rcf2)

def BHCF( dL0, dL1, dL2, gamma ):
    hcf0 = HCF(dL0,dL1)
    hcf1 = HCF(dL1,dL2)
    return DiscontinuousUnivariateFunction(gamma,hcf0,hcf1)

def BHRCF( beta, dL0, dL1, gamma ):
    rcf = RCF(0,1,beta)
    hcf = HCF(dL0,dL1)
    return DiscontinuousUnivariateFunction(gamma,rcf,hcf)

def VCF(dL0, dL1, L, n):
    return ValliammaiFunction(dL0,dL1,L,n)

def BVCF(beta0,beta1,beta2,gamma,R,n):
    VCF_a = VCF(beta0*R,beta1*R,R*gamma,int(n*gamma))
    VCF_b = VCF(beta1*R,beta2*R,R*(1-gamma),int(n*(1-gamma)))
    return DiscontinuousUnivariateFunction(gamma,VCF_a,VCF_b)
%}
