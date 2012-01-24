## \file simple_cylinders.py
## \ingroup geom
## \brief Some handy functions for geometry specification with cylinders.
##
## It is intended that a user script could import this module
## and make use of the functions to generate lists of block edges.

from geom import *
from gpath import *
from geom_elmer import *

def simpleCylinder(xPos=0.0, yPos=0.0, zPos=0.0, Rmin=1.0, Rmax=2.0, ySize=1.0):
    """
    Creates an edge list for a quarter cylinder.
    
    See workbook page 45, 28-Nov-04 for orientation.
    The end-face is in the (x,z)-plane and the cylinder axis is along the y-axis.
    """
    p0 = Node(xPos-Rmin, yPos, zPos + 0.0)
    p1 = Node(xPos, yPos, zPos+Rmin)
    p2 = Node(xPos, yPos+ySize, zPos+Rmin)
    p3 = Node(xPos-Rmin, yPos+ySize, zPos)
    pc01 = Node(xPos, yPos, zPos)
    
    p4 = Node(xPos-Rmax, yPos, zPos + 0.0)
    p5 = Node(xPos, yPos, zPos+Rmax)
    p6 = Node(xPos, yPos+ySize, zPos+Rmax)
    p7 = Node(xPos-Rmax, yPos+ySize, zPos)
    pc32 = Node(xPos, yPos+ySize, zPos)

    pl0 = Edge3D(Arc(p0, p1, pc01))
    pl1 = Edge3D(Line(p1, p2))
    pl2 = Edge3D(Arc(p3, p2, pc32))
    pl3 = Edge3D(Line(p0, p3))
    pl4 = Edge3D(Arc(p4, p5, pc01))
    pl5 = Edge3D(Line(p5, p6))
    pl6 = Edge3D(Arc(p7, p6, pc32))
    pl7 = Edge3D(Line(p4, p7))
    pl8 = Edge3D(Line(p0, p4))
    pl9 = Edge3D(Line(p1, p5))
    pl10 = Edge3D(Line(p2, p6))
    pl11 = Edge3D(Line(p3, p7))
    return [pl0, pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, pl10, pl11]


def simpleCylinderExtruded(xPos=0.0, yPos=0.0, zPos=0.0, Rmin=1.0, Rmax=2.0, ySize=1.0):
    """
    Creates an edge list for a quarter cylinder.
    
    This function first creates the end-surface and then
    extrudes it along the cylinder axis in the j-direction
    (which happens to coincide with the y-axis).
    """
    p0 = Node(xPos-Rmin, yPos, zPos + 0.0)
    p1 = Node(xPos, yPos, zPos+Rmin)
    p4 = Node(xPos-Rmax, yPos, zPos + 0.0)
    p5 = Node(xPos, yPos, zPos+Rmax)
    pc01 = Node(xPos, yPos, zPos)  # centre of circular paths

    pl_01 = Edge3D(Arc(p0, p1, pc01))
    pl_45 = Edge3D(Arc(p4, p5, pc01))
    pl_04 = Edge3D(Line(p0, p4))
    pl_15 = Edge3D(Line(p1, p5))
    
    p3 = Node(xPos-Rmin, yPos+ySize, zPos)
    pl_03 = Edge3D(Line(p0, p3))

    end_surface = ClosedSurfacePatch(pl_01, pl_45, pl_04, pl_15)
    edge_list = end_surface.extrude(pl_03, "j")
    return edge_list


