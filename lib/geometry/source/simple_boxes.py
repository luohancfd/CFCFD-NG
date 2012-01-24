## \file simple_boxes.py
## \ingroup geom
## \brief Some handy functions for geometry specification with boxes.
##
## It is intended that a user script could import this module
## and make use of the functions to generate lists of block edges.

from geom import *
from gpath import *
from geom_elmer import *

def simpleBoxCorners(xPos=0.0, yPos=0.0, zPos=0.0, xSize=1.0, ySize=1.0, zSize=1.0):
    """\brief Creates a corner coordinate list for a simple box."""
    p0 = Node(xPos,       yPos,       zPos)
    p1 = Node(xPos+xSize, yPos,       zPos)
    p2 = Node(xPos+xSize, yPos+ySize, zPos)
    p3 = Node(xPos,       yPos+ySize, zPos)
    p4 = Node(xPos,       yPos,       zPos+zSize)
    p5 = Node(xPos+xSize, yPos,       zPos+zSize)
    p6 = Node(xPos+xSize, yPos+ySize, zPos+zSize)
    p7 = Node(xPos,       yPos+ySize, zPos+zSize)
    return [p0, p1, p2, p3, p4, p5, p6, p7]

def simpleBoxEdges(p):
    """\brief Creates an edge list from a corner coordinate list for a simple box."""
    pl0 = Edge3D(Line(p[0], p[1]))
    pl1 = Edge3D(Line(p[1], p[2]))
    pl2 = Edge3D(Line(p[3], p[2]))
    pl3 = Edge3D(Line(p[0], p[3]))
    pl4 = Edge3D(Line(p[4], p[5]))
    pl5 = Edge3D(Line(p[5], p[6]))
    pl6 = Edge3D(Line(p[7], p[6]))
    pl7 = Edge3D(Line(p[4], p[7]))
    pl8 = Edge3D(Line(p[0], p[4]))
    pl9 = Edge3D(Line(p[1], p[5]))
    pl10 = Edge3D(Line(p[2], p[6]))
    pl11 = Edge3D(Line(p[3], p[7]))
    return [pl0, pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, pl10, pl11]

def simpleBox(xPos=0.0, yPos=0.0, zPos=0.0, xSize=1.0, ySize=1.0, zSize=1.0):
    """\brief Compute box edges from sizes."""
    return simpleBoxEdges(simpleBoxCorners(xPos, yPos, zPos, xSize, ySize, zSize))

