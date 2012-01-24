# minimal_geometry.py
#-----------------------------------------
# A bare minimum geometry library to do
# some of the work required by Rowan's laura2vtk.py
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
