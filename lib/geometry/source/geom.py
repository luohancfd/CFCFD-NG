## \file geom.py
## \ingroup geom
## \brief Python geometry-specification functions for mb_cns and Elmer.
##
## \author P.Jacobs and R. Gollan
## \version 1.0 August-November 2004
## \version 1.1 31-January-2005  -- independent of mb_cns and Elmer
## \version 1.2 1-Feb-2005 -- include Rowan's functions to overload operators
## \version 1.3 2-Dec-2005 -- pure-Python option so that we no longer need to
##                            be tied to the libgeom C library.
"""
Provides basic 3D/2D L{Vector} and L{Node} classes for constructing
geometric descriptions.

For 2D modelling, the z-coordinate can be omitted so that it takes
its default value of 0.0.
"""

#---------------------------------------------------------------------

def set_libgeom_flag(value):
    """
    Adjusts the flag that determines if we want calls to to go to
    the C-version of libgeom or to the pure Python code instead.

    @param value: 1=use C libgeom calls; 1=use Python only
    @type value: int
    """
    global use_libgeom_C_flag
    use_libgeom_C_flag = value
    return

def use_libgeom():
    global use_libgeom_C_flag
    return use_libgeom_C_flag

try:
    # print "Loading geom.py"
    import libgeom
    # print "Loaded libgeom successfully."
    set_libgeom_flag(1)  # default is to use C lib if available
except:
    print "Did not load libgeom C module successfully."
    set_libgeom_flag(0)

from math import sqrt

#----------------------------------------------------------------------

class Vector(object):
    """
    Defines a vector 3D space.

    The vector is created in the C module data space and
    new-style object properties to access the C-module values.

    @undocumented: __del__, getX, setX, getY, setY, getZ, setZ,
        getLabel, setLabel
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """

    __slots__ = 'x', 'y', 'z', 'pp', 'label', 'label_string'
    
    def __init__(self, x=0.0, y=0.0, z=0.0, label=""):
        """
        Create a Vector from its Cartesian components.

        @param x: x-component
        @type x: float       
        @param y: y-component
        @type y: float       
        @param z: z-component
        @type z: float
        @param label: optional text label that will appear in the MetaPost file.
        @type label: string
        """
        if use_libgeom():
            self.pp = libgeom.create_point_3D(x, y, z)
            # Note that we do not release memory for Vector in C module just in case
            # the C functions try to use the points after the python interpreter
            # has thrown away the Vector object.
            # Yes, this is a memory leak but we shall tolerate it.
            if self.pp == None:
                print "Creation of point at (%e, %e, %e) failed." % (x, y, z)
                return
        else:
            self.pp = [x, y, z]
        self.label_string = label
        return

    def setX(self, x):
        if use_libgeom():
            self.pp.x = float(x)
        else:
            self.pp[0] = float(x)
        return
    def getX(self):
        if use_libgeom():
            return self.pp.x
        else:
            return self.pp[0]
    x = property(getX, setX, doc="x-coordinate of the vector")

    def setY(self, y):
        if use_libgeom():
            self.pp.y = float(y)
        else:
            self.pp[1] = float(y)
        return
    def getY(self):
        if use_libgeom():
            return self.pp.y
        else:
            return self.pp[1]
    y = property(getY, setY, doc="y-coordinate of the vector")

    def setZ(self, z):
        if use_libgeom():
            self.pp.z = float(z)
        else:
            self.pp[2] = float(z)
        return
    def getZ(self):
        if use_libgeom():
            return self.pp.z
        else:
            return self.pp[2]
    z = property(getZ, setZ, doc="z-coordinate of the vector")

    def setLabel(self, label):
        if use_libgeom():
            libgeom.point_3D_set_label(self.pp, str(label))
        else:
            self.label_string = str(label)
        return
    def getLabel(self):
        if use_libgeom():
            return libgeom.point_3D_get_label(self.pp)
        else:
            return self.label_string
    label = property(getLabel, setLabel, doc="text label for the vector")
    
    def __str__(self):
        return "Vector(%g,%g,%g,\"%s\")" % (self.x, self.y, self.z, self.label)

    def __repr__(self):
        return self.__str__()

    def copy(self):
        """
        Create and return a new instance which is a copy of the original.
        """
        return Vector(self.x, self.y, self.z, self.label)

    def __add__(self, other):
        """Element by element addition."""
        if not isinstance(other,Vector):
            raise TypeError, "Should have supplied Vectors to +."
        return Vector(self.x+other.x, self.y+other.y, self.z+other.z)

    def __neg__(self):
        """Negation of all elements for -x."""
        return Vector(-self.x, -self.y, -self.z)

    def __pos__(self):
        """+x"""
        return Vector(self.x, self.y, self.z)

    def __sub__(self, other):
        """Element-by-element subtraction."""
        return self + (-other)

    def __mul__(self, other):
        """Element-by-element multiplication."""
        if isinstance(other,Vector):
            return Vector(self.x * other.x, self.y * other.y, self.z * other.z)
        else:
            return Vector(self.x * other, self.y * other, self.z * other)

    def __rmul__(self, other):
        return (self * other)

    def __div__(self, other):
        raise NotInplementedError, "Division not implemented for Vectors"

    def __rdiv__(self, other):
        raise NotInplementedError, "Division not implemented for Vectors"

    def sum(self):
        """Add elements together."""
        return (self.x + self.y + self.z)

    def __abs__(self):
        """Absolute value is the magnitude of the Vector."""
        return sqrt((self * self).sum())

    def unit(self):
        """@return: unit vector with same direction as this vector."""
        mag = abs(self)
        if mag <= 0.0:
            raise ValueError, "Zero magnitude vector has no defined direction."
        return Vector(self.x/mag, self.y/mag, self.z/mag)

def dot(a, b):
    """
    Vector dot product.
    
    @returns: scalar product.
    """
    return (a * b).sum()

def cross(a, b):
    """
    Vector cross product.

    @returns: vector product.
    """
    x = a.y * b.z - a.z * b.y
    y = a.z * b.x - a.x * b.z
    z = a.x * b.y - a.y * b.x
    return Vector(x, y, z)

#---------------------------------------------------------------------

class Node(Vector):
    """
    Defines a nodal-point in 3D space to be subsequently used
    in the definition of line segments.

    @undocumented: label, nodeList
    @undocumented: __delattr__, __getattribute__, __hash__, __new__,
        __reduce__, __reduce_ex__, __repr__, __setattr__, __str__
    """
    nodeList = []  # a collection so that we can later render the nodes
                   # to MetaPost or whatever.
    
    def __init__(self, x=0.0, y=0.0, z=0.0, label=""):
        Vector.__init__(self, x, y, z, label)
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
        if isinstance(dx, Vector):
            dx, dy, dz = dx.x, dx.y, dx.z
        if use_libgeom():
            libgeom.point_3D_translate(self.pp, dx, dy, dz)
        else:
            self.x += dx; self.y += dy; self.z += dz
        return
    
#----------------------------------------------------------------------

def distance_between_nodes(a, b):
    """Return the distance between Nodes a and b."""
    return abs(a - b)

#----------------------------------------------------------------------

if __name__ == '__main__':
    print "Begin test of geom.py"
    set_libgeom_flag(0)
    print "use_libgeom()=", use_libgeom()

    print "Part 1: Vectors"
    for i in range(100):
        a = Vector(1,2,0)
    # Hopefully the garbage collector would be exercised.
    print "Examples: a=", a, "2*(+a)=", 2*(+a)
    b = a + a + a
    print "After sum: b=", b
    c = b - Vector(1.0, 2.2, 3.14159)
    print "Subtraction: c=", c
    print "abs(a)=", abs(a)
    print "sqrt(dot(a,a))=", sqrt(dot(a,a))
    d = cross(a,b)
    print "a cross b = ", d
    print "unit-v of c=", c.unit()
    print "abs(c.unit())=", abs(c.unit())
    
    print "Part 2: Nodes"
    na = Node(10, 20, 30, "point-a")
    nb = Node(-2, -4)
    print "Examples: na=", na, "nb=", nb
    print "Vector add: na+a=", na + a
    nc = na.copy()
    nc.label = "new-point-c"
    print "Copy and relabel: nc=", nc, "na=", na
    na.translate(2*a)
    print "after Vector translate: na=", na
    na.translate(0.1, 0.2, 0.3)
    print "after component translate: na=", na
    
    print "End test of geom.py"
