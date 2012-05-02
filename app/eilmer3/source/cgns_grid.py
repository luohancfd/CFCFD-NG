# cgns_grid.py
"""
Read structured grids from CGNS files.

Authors: Paul Petrie-Repar, QGECE, 15 Feb 2011
         Peter J, 23 Feb 2011
         Sancho J, 8 Feb 2012

This module makes use of Oliver Borm's CGNS package for Python 2.x
http://sourceforge.net/projects/python-cgns/

We have a copy of that package in cfcfd3/extern/python-cgns/.

Sancho J just introduced some modifications to read 2D meshes
Sancho J moved bc translation part to main eilmer3 script
"""

from CGNS import CGNS
import sys, os
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid


def check_CGNS_error_flag(errorFlag, messageText=""):
    """
    Give some information on CGNS failure.
    """
    if (errorFlag != 0):
        print "CGNS Error %d %s" % (errorFlag, messageText)
        CGNS.cg_error_exit()
    return

def getFaceName(face, bcPnts):
    
    if (bcPnts[0] == bcPnts[3]):
        if (bcPnts[0] == 1):
            return "WEST"
        else: 
            return "EAST"
    elif (bcPnts[1] == bcPnts[4]):
        if (bcPnts[1] == 1):
            return "SOUTH"
        else: 
            return "NORTH"
    elif (bcPnts[2] == bcPnts[5]):
        if (bcPnts[2] == 1):
            return "BOTTOM"
        else: 
            return "TOP"
    else:
        print "Cannot determine Face"

# modification by jorge sancho to get 2D meshes---------------
def getFaceName2D(face, bcPnts):
    
    if (bcPnts[0] == bcPnts[2]):
        if (bcPnts[0] == 1):
            return "WEST"
        else: 
            return "EAST"
    elif (bcPnts[1] == bcPnts[3]):
        if (bcPnts[1] == 1):
            return "SOUTH"
        else: 
            return "NORTH"
    else:
        print "Cannot determine Face"
#--------------------------------------------------------------

def read_ICEM_CGNS_grids(cgnsFileName, labelStem='test', gridScale=1.0):
    """
    Dip into the CGNS file and extract the structured grids from within.

    Returns a dictionary containing the interesting data.
    """
    cgnsData = dict() # a convenient place to collect everything

    print "open CGNS file: ", cgnsFileName, " scale: ", gridScale
    filePtr = CGNS.intp()
    ier = CGNS.cg_open(cgnsFileName, CGNS.CG_MODE_READ, filePtr)
    fileValue = filePtr.value()
    check_CGNS_error_flag(ier, "on opening file")

    basePtr = CGNS.intp()
    basePtr.assign(1)
    numberZonesPtr = CGNS.intp()
    ier = CGNS.cg_nzones(filePtr.value(), basePtr.value(), numberZonesPtr)
    check_CGNS_error_flag(ier, "on getting number of zones")
    numberZones = numberZonesPtr.value()
    cgnsData['nblock'] = numberZones
    cgnsData['grids'] = []
    cgnsData['bcs'] = []
    # print "cgnsData=", cgnsData

    zonePtr = CGNS.intp()
    zoneSizeArray = CGNS.intArray(9)
    range_min = CGNS.intArray(3)
    range_max = CGNS.intArray(3)

    for zone in range(numberZones):
        zonePtr.assign(zone+1)
        ier, zoneName = CGNS.cg_zone_read(filePtr.value(), basePtr.value(), zonePtr.value(),
                                          zoneSizeArray)
        check_CGNS_error_flag(ier, "on getting zone pointer")
        #
        ni = zoneSizeArray[0]; nj = zoneSizeArray[1]; nk = zoneSizeArray[2]
        if (zoneSizeArray[4] == 0): #modification by JorgeSancho to include 2D meshes
            nk = 1
#            print "This block is 2D" 
        # print "ni=", ni, "nj=", nj, "nk=", nk
        range_min[0] = 1; range_max[0] = ni
        range_min[1] = 1; range_max[1] = nj
        range_min[2] = 1; range_max[2] = nk
        # Read the coordinate arrays.
        numberNodes = ni * nj * nk
        coordX = CGNS.doubleArray(numberNodes)
        coordY = CGNS.doubleArray(numberNodes)
        coordZ = CGNS.doubleArray(numberNodes)
        ier, coordName = CGNS.cg_coord_read(filePtr.value(), basePtr.value(), zonePtr.value(),
                                            "CoordinateX", CGNS.RealDouble, 
                                            range_min, range_max, coordX)
        check_CGNS_error_flag(ier, "on getting x-coordinates")
        ier, coordName = CGNS.cg_coord_read(filePtr.value(), basePtr.value(), zonePtr.value(),
                                            "CoordinateY", CGNS.RealDouble, 
                                            range_min, range_max, coordY)
        check_CGNS_error_flag(ier, "on getting y-coordinates")
        ier, coordName = CGNS.cg_coord_read(filePtr.value(), basePtr.value(), zonePtr.value(),
                                            "CoordinateZ", CGNS.RealDouble, 
                                            range_min, range_max, coordZ)
        check_CGNS_error_flag(ier, "on getting z-coordinates")
        # Now that we've read the coordinates, repack them into Eilmer's data structure.
        g = StructuredGrid((ni,nj,nk), label='%s%04d'%(labelStem, zone))
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    indx = k*nj*ni + j*ni + i
                    g.x[i,j,k] = coordX[indx] * gridScale
                    g.y[i,j,k] = coordY[indx] * gridScale
                    g.z[i,j,k] = coordZ[indx] * gridScale
        cgnsData['grids'].append(g)

###     read boundary conditions

        if (nk == 1): #modification by JorgeSancho to include 2D meshes
            count = 4 
        else: 
            count = 6

        for face in range(count):
            ier = CGNS.cg_goto(fileValue, basePtr.value(), "Zone_t", (zone + 1), "ZoneBC_t", 1, "BC_t", (face+1), "end")
            ier, familyName = CGNS.cg_famname_read()

            bcPntsPointer = CGNS.intArray(count)
            normalListPointer = CGNS.intArray(3)
            CGNS.cg_boco_read(fileValue, basePtr.value(), (zone+1), (face+1), bcPntsPointer,  normalListPointer);

            print zone, face, familyName    
            bc = {}
            bc['block'] = zone

            if (count == 4): #modification by JorgeSancho to include 2D meshes
                bc['face'] = getFaceName2D(face, bcPntsPointer)
            else: 
                bc['face'] = getFaceName(face, bcPntsPointer)
#
#          do not add internal and wall boundaries
#
#          J.Sancho modification: Now bocos are defined in main file
           
            familyName.capitalize()
            bc['type'] =familyName

        # end for zone...
    ier = CGNS.cg_close(filePtr.value())
    check_CGNS_error_flag(ier, "on closing file")

    if (nk == 1): #modification by JorgeSancho to include 2D meshes
        print 'The mesh is 2d'
    print cgnsData['bcs']

    return cgnsData


if __name__ == '__main__':
    print "cgns_import.py Demonstration..."
    if len(sys.argv) < 2:
        print "Usage: python import_cgns_grid.py CGNS_filename"
        sys.exit()
    fileName = sys.argv[1]
    print "CGNS file name: ", fileName
    dataDict = read_ICEM_CGNS_grids(fileName)
    nb = dataDict['nblock']
    print "number of blocks found=", nb
    print "writing VTK grid files:"
    for jb in range(nb):
        g = dataDict['grids'][jb]
        print '   ', g.label, 'ni=', g.ni, 'nj=', g.nj, 'nk=', g.nk
        f = open(g.label + ".vtk", 'w')
        g.write_block_in_VTK_format(f)
        f.close()
    print "Done."

