#!/usr/bin/env python
# convert Eilmer3 solution to CGNS
# started by Paul Petrie-Repar 13 Jul 10 for UQ

from CGNS import CGNS

import ConfigParser
import sys, os, gzip
from getopt import getopt

#from bc_defs import BoundaryCondition
sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow

    

def getNumberVerticesAtFace(grid, face):
    num = 0
    if (face == "east" or face == "west"):
        num = grid.nj * grid.nk
    elif (face == "north" or face == "south"):
        num = grid.ni * grid.nk
    elif (face == "top" or face == "bottom"):
        num = grid.ni * grid.nk
    return num

def getConnVerts(grid, face, index):
#     index = i2 * n1 + i1    
    n1 = 0
    if (face == "east" or face == "west"):
        n1 = grid.nj
    elif (face == "north" or face == "south"):
        n1 = grid.nk
    elif (face == "top" or face == "bottom"):
        n1 = grid.ni
        
    i2 = int(index/n1)
    i1 = index - i2 * n1

    i=0
    j=0
    k = 0

    if (face == "west"):
        i=0
        j=i1
        k=i2
    elif (face == "east"):
        i=grid.ni-1
        j=i1
        k=i2
    elif (face == "south"):
        i=i2
        j=0
        k=i1
    elif (face == "north"):
        i=i2
        j=grid.nj-1
        k=i1
    elif (face == "bottom"):
        i=i1
        j=i2
        k=0
    elif (face == "top"):
        i=i1
        j=i2
        k=grid.nk-1

#    print index, n1, i1, i2, i,j,k
         
    return i,j,k


def getFaceRanges(grid, face):
    if (face == 'east'):
        ilo = grid.ni - 1
        ihi = grid.ni - 1
        jlo = 0
        jhi = grid.nj - 1
        klo = 0
        khi = grid.nk - 1
    elif (face == 'west'):
        ilo = 0
        ihi = 0
        jlo = 0
        jhi = grid.nj - 1
        klo = 0
        khi = grid.nk - 1
    elif (face == 'south'):
        ilo = 0
        ihi = grid.ni -1
        jlo = 0
        jhi = 0
        klo = 0
        khi = grid.nk - 1
    elif (face == 'north'):
        ilo = 0
        ihi = grid.ni -1
        jlo = grid.nj - 1
        jhi = grid.nj - 1
        klo = 0
        khi = grid.nk - 1
    elif (face == 'bottom'):
        ilo = 0
        ihi = grid.ni - 1
        jlo = 0
        jhi = grid.nj - 1
        klo = 0
        khi = 0
    elif (face == 'top'):
        ilo = 0
        ihi = grid.ni - 1
        jlo = 0
        jhi = grid.nj - 1
        klo = grid.nk - 1 
        khi = grid.nk - 1 
    
    return ilo, ihi, jlo, jhi, klo, khi 


def addBoundary(bcs, block, face, name):
    bcs['numberBoundaries'] = bcs['numberBoundaries'] + 1
    bcs['block'].append(block) 
    bcs['name'].append(name)
    bcs['face'].append(face) 

def addBlockConnection(conns, block, face, connectedBlock, connFace):

    conns['block_left'].append(block)
    conns['block_right'].append(connectedBlock)

    conns['face_left'].append(face)
    conns['face_right'].append(connFace)
    print "Join", conns['numberConnections'], "block" , block, "face", face, "connected to block", connectedBlock, "face", connFace

    conns['numberConnections'] = conns['numberConnections'] + 1 


def readConfigFile(fileName, blocks):
    bcs = dict()
    blockConnections = dict()
    
    bcs['numberBoundaries'] = 0
    
    bcs['block'] = []
    bcs['face'] = []
    bcs['name'] = []
    bcs['label'] = []

    blockConnections['numberConnections'] = 0
     
    blockConnections['block_left'] = []
    blockConnections['block_right'] = []

    blockConnections['face_left'] = []
    blockConnections['face_right'] = []
     

    
    Config = ConfigParser.ConfigParser()
    Config.read(fileName)

    blockFaces = ['east', 'west', 'south', 'north', 'bottom', 'top']
    wallCount = 0
    luaCount = 0
    for block in blocks:
        for face in blockFaces:
            section = "block/" + str(block) + "/face/" + face
            bnd_type = Config.get(section, "bc")
            isWall = Config.get(section, "is_wall")
            connectedBlock = int(Config.get(section, "other_block"))
            bclabel = Config.get(section, "label")
            if (isWall == '1'):
                if (bclabel == ""):
                    bclabel = "Wall" + str(wallCount)
                    wallCount = wallCount+1
                addBoundary(bcs, block, face, bclabel)
            elif (bnd_type == '16'):
                if (bclabel == ""):
                    bclabel = "Lua" + str(luaCount)
                    luaCount = luaCount+1
                addBoundary(bcs, block, face, bclabel)
            elif (connectedBlock >= 0):
                if (connectedBlock > block):  # we add each connection only once
                    connFace = Config.get(section, "other_face")
                    neighbour_orientation = int(Config.get(section, "neighbour_orientation"))
                    if (neighbour_orientation != 0):
                        print "Neighbour Orientation: ", neighbour_orientation
                        print "Add code here"
                        sys.exit()
                    else:
                        addBlockConnection(blockConnections, block, face, connectedBlock, connFace)
            else:
                print "Unknown face connection or boundary condition for", section
                sys.exit()
                
    print blockConnections
  
    return bcs, blockConnections 

def determineSign(low, high):
    if (high < low):
        return -1
    
    return 1    
    
def checkCGNSerror(ier):
    if (ier != 0):
        CGNS.cg_error_print()

def getCGNSboundaryType(name):
    if(name == "inlet"):
        return CGNS.BCInflow
    elif (name == "outlet"):
        return CGNS.BCOutflow
    elif (name.find("blade") != -1 or name.find("hub") != -1 or name.find("shroud") != -1):
        return CGNS.BCWall
    else:
        return CGNS.BCTypeNull

def getBlockLabels(fileName): 
    Config = ConfigParser.ConfigParser()
    Config.read(fileName)
    numberBlocks = int(Config.get("global_data", "nblock"))
    print "Number of Blocks:", numberBlocks
    blockLabels = []

    blocks = range(0, numberBlocks)
    for block in blocks:
        section = "block/" + str(block) 
        blockLabels.append(Config.get(section, "label"))

    print blockLabels
    return blockLabels

def printUsage():
    print ""
    print "Usage: elm2cgns.py" + \
          " jobName" + \
          " cgnsFile" + \
          " -s startBlock -e endBlock]" 
    print ""
    return

def writeZone(fileValue, zoneName, blocks, boundaryConditions, blockConnections):

    tindx = '9999'
    
    flow = []
    grid = []

    for block_index in range(len(blocks)):
    
        block = blocks[block_index]
        bindx = '000' + str(block) 
    
    #  read blocks
    
        fileName = 'flow/t'+tindx+'/'+ jobName + '.flow.b'+bindx+'.t'+tindx+'.gz'
        fp = gzip.open(fileName, "r")
        flow.append(StructuredGridFlow())
        flow[block_index].read(fp)
        fp.close()
        print "flow data: ni=", flow[block_index].ni, "nj=", flow[block_index].nj, "nk=", flow[block_index].nk
        # The grid is always at tindx 0.
        fileName = 'grid/t0000/'+ jobName + '.grid.b'+bindx+'.t0000.gz'
        fp = gzip.open(fileName, "r")
        grid.append(StructuredGrid())
        grid[block_index].read(fp)
        fp.close()
        print "Block ", block, "grid data: ni=", grid[block_index].ni, "nj=", grid[block_index].nj, "nk=", grid[block_index].nk
    
    numberCells = 0
    numberStructuredPoints = 0
    startBlockVertex =[]
    startBlockCell = []
    for block_index in range(len(blocks)):
    
        startBlockVertex.append(numberStructuredPoints)
        startBlockCell.append(numberCells)  
        numberCells = numberCells + flow[block_index].ni * flow[block_index].nj * flow[block_index].nk
        numberStructuredPoints = numberStructuredPoints + grid[block_index].ni * grid[block_index].nj * grid[block_index].nk
    
    print "Number Cells: ", numberCells
    print "Number Structured Vertices: ", numberStructuredPoints
        
    # copy co-ordinates to lists
    
    structured_x = []
    structured_y = []
    structured_z = []
    
    for block_index in range(len(blocks)):
        for k in range(grid[block_index].nk):
            for j in range(grid[block_index].nj):
                for i in range(grid[block_index].ni):
                    structured_x.append(grid[block_index].x[i,j,k]) 
                    structured_y.append(grid[block_index].y[i,j,k]) 
                    structured_z.append(grid[block_index].z[i,j,k])
    
    # Hard wire block connections should be ready in"
    
    
    
    
    #create Unstructured Vertex Map
    
    #print blockConnections
    
    vertexMap = range(numberStructuredPoints)
    
    for join in range(blockConnections['numberConnections']):
        
    
        block_left = blockConnections['block_left'][join]
        block_left_index = blocks.index(block_left)         
        block_right = blockConnections['block_right'][join]
        block_right_index = blocks.index(block_right)         
    
        face_left = blockConnections['face_left'][join]
        face_right = blockConnections['face_right'][join]
    
        numberJoinVertices_left = getNumberVerticesAtFace(grid[block_left_index], face_left)
        numberJoinVertices_right = getNumberVerticesAtFace(grid[block_right_index], face_right)
       
        try:
            (numberJoinVertices_left == numberJoinVertices_right)
        except:
            print "Error"
        else:
            print "Join : ", join, " ", numberJoinVertices_left, numberJoinVertices_right
    
        n_i_l = grid[block_left_index].ni  
        n_ij_l = n_i_l * grid[block_left_index].nj 
    
        n_i_r = grid[block_right_index].ni  
        n_ij_r = n_i_r * grid[block_right_index].nj 
    
    
        for i in range(numberJoinVertices_left):   
            i_l, j_l, k_l = getConnVerts(grid[block_left_index], face_left, i)
            i_r, j_r, k_r = getConnVerts(grid[block_right_index], face_right, i)
    
                    
            vert_l = startBlockVertex[block_left_index] + k_l * n_ij_l + j_l * n_i_l + i_l 
            vert_r = startBlockVertex[block_right_index] + k_r * n_ij_r + j_r * n_i_r + i_r
                    
    #                print "Map vert ", vert_l, " to ", vert_r     
            vertexMap[vert_l] = vert_r
    
    #create list of unstructuredPoints and remove multiple references in vertMap  
    
    #vertexMap = createVertexMap(structured_x, structured_y, structured_z)
    
    
    
    unstructuredPoints = []
    
    for i in range(len(vertexMap)):
        if ( vertexMap[i] == i):
            unstructuredPoints.append(i)
        else:
            x1 = 0
            refPoint = vertexMap[i] 
    #        print "Ref", i, refPoint, vertexMap[refPoint] 
            while (refPoint != vertexMap[refPoint]):
                refPoint = vertexMap[refPoint]
    #            print "double ref", x1, refPoint
                try:
                    (x1 < 8)
                except:
                    print "Error nest too deep"
            vertexMap[i] = refPoint
           
                    
            
    numberUnstructuredPoints = len(unstructuredPoints) 
    
    print "Number Unstructured VerticesCells: ", numberUnstructuredPoints
    
    # create inverse Map
    
    inverseMap = range(numberStructuredPoints)
    for i in range(numberUnstructuredPoints):
        inverseMap[unstructuredPoints[i]] = i
    
    for i in range(numberStructuredPoints):
        if (i != vertexMap[i]):
            inverseMap[i] = inverseMap[vertexMap[i]]
        
        
        
    
    #create zone
    
    zoneArrayp= CGNS.intArray(3)
    #vertex size
    zoneArrayp[0] = numberUnstructuredPoints
    zoneArrayp[1] = numberCells
    zoneArrayp[2] = 0
    
    zonePointer = CGNS.intp()
    ier = CGNS.cg_zone_write ( fileValue, bValue, zoneName, zoneArrayp, CGNS.Unstructured, zonePointer )
    checkCGNSerror(ier)
    zoneValue = zonePointer.value()
    
    #write co-ordinates
    
    CoordX= CGNS.doubleArray(numberUnstructuredPoints)
    CoordY= CGNS.doubleArray(numberUnstructuredPoints)
    CoordZ= CGNS.doubleArray(numberUnstructuredPoints)
    
    for i in range(numberUnstructuredPoints):
        CoordX[i] = structured_x[unstructuredPoints[i]]
        CoordY[i] = structured_y[unstructuredPoints[i]]
        CoordZ[i] = structured_z[unstructuredPoints[i]]
    
    
    coordPointer = CGNS.intp()
    ier = CGNS.cg_coord_write ( fileValue, bValue, zoneValue, CGNS.RealDouble, "CoordinateX", CoordX, coordPointer )
    checkCGNSerror(ier)
    ier = CGNS.cg_coord_write ( fileValue, bValue, zoneValue, CGNS.RealDouble, "CoordinateY", CoordY, coordPointer )
    checkCGNSerror(ier)
    ier = CGNS.cg_coord_write ( fileValue, bValue, zoneValue, CGNS.RealDouble, "CoordinateZ", CoordZ, coordPointer )
    checkCGNSerror(ier)
    
    #write cell connections
    
    Elements= CGNS.intArray(numberCells * 8)
    x1 = 0
    for block_index in range(len(blocks)):
        n_i = grid[block_index].ni  
        n_ij = n_i * grid[block_index].nj 
    
        for k in range(flow[block_index].nk):
            for j in range(flow[block_index].nj):
                for i in range(flow[block_index].ni):
    
                    structPoint = startBlockVertex[block_index] + k * n_ij + j * n_i + i
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + k * n_ij + j * n_i + i + 1
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + k * n_ij + (j + 1) * n_i + i + 1
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + k * n_ij + (j + 1) * n_i + i 
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + (k + 1) * n_ij + j * n_i + i
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + (k + 1) * n_ij + j * n_i + i + 1
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + (k + 1) * n_ij + (j + 1) * n_i + i + 1
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
                    structPoint = startBlockVertex[block_index] + (k + 1) * n_ij + (j + 1) * n_i + i 
                    unstrPoint = inverseMap[structPoint]
                    Elements[x1] = unstrPoint + 1
                    x1 = x1 + 1
    
    
    sectionPointer = CGNS.intp()
    ier = CGNS.cg_section_write(fileValue, bValue, zoneValue, "Cells", CGNS.HEXA_8, 1, numberCells, 0, Elements, sectionPointer)
    checkCGNSerror(ier)

    
    # write solution
    
    solnPointer = CGNS.intp()
    ier = CGNS.cg_sol_write ( fileValue, bValue, zoneValue, "FlowSolution", CGNS.CellCenter, solnPointer ) 
    checkCGNSerror(ier)

    solnValue = solnPointer.value()
    
    
    density = CGNS.doubleArray(numberCells)
    velx = CGNS.doubleArray(numberCells)
    vely = CGNS.doubleArray(numberCells)
    velz = CGNS.doubleArray(numberCells)
    pressure = CGNS.doubleArray(numberCells)
    
    x1 = 0
    for block_index in range(len(blocks)):
        n_i = grid[block_index].ni  
        n_ij = n_i * grid[block_index].nj 
    
        for k in range(flow[block_index].nk):
            for j in range(flow[block_index].nj):
                for i in range(flow[block_index].ni):
    
                    density[x1] = flow[block_index].data['rho'][i,j,k] 
                    velx[x1] = flow[block_index].data['vel.x'][i,j,k] 
                    vely[x1] = flow[block_index].data['vel.y'][i,j,k] 
                    velz[x1] = flow[block_index].data['vel.z'][i,j,k] 
                    pressure[x1] = flow[block_index].data['p'][i,j,k] 
                    x1 = x1 + 1
    
    fieldPointer = CGNS.intp()
    ier = CGNS.cg_field_write ( fileValue, bValue, zoneValue, solnValue, CGNS.RealDouble, "Density", density, fieldPointer )
    ier = CGNS.cg_field_write ( fileValue, bValue, zoneValue, solnValue, CGNS.RealDouble, "VelocityX", velx, fieldPointer )
    ier = CGNS.cg_field_write ( fileValue, bValue, zoneValue, solnValue, CGNS.RealDouble, "VelocityY", vely, fieldPointer )
    ier = CGNS.cg_field_write ( fileValue, bValue, zoneValue, solnValue, CGNS.RealDouble, "VelocityZ", velz, fieldPointer )
    ier = CGNS.cg_field_write ( fileValue, bValue, zoneValue, solnValue, CGNS.RealDouble, "Pressure", pressure, fieldPointer )
    
    
    #write boundary conditions to CGNS file
    
    element_index = numberCells
    
    for bnd in range(boundaryConditions['numberBoundaries']):
    
        block = boundaryConditions['block'][bnd]
        block_index = blocks.index(block)
        face = boundaryConditions['face'][bnd]
        bndName = boundaryConditions['name'][bnd]
    
        ilo, ihi, jlo, jhi, klo, khi = getFaceRanges(grid[block_index], face)
        
    
        i_size = abs(ihi - ilo) #  cells not vertices
        j_size = abs(jhi - jlo) #  cells not vertices 
        k_size = abs(khi - klo)  #  cells not vertices
        
        i_sign = determineSign(ilo, ihi) 
        j_sign = determineSign(jlo, jhi)
        k_sign = determineSign(klo, khi)
    
        n_i = grid[block_index].ni  
        n_ij = n_i * grid[block_index].nj 
    
       
        if (face == "west" or face == "east"):
            numberCellsAtBoundary =  j_size * k_size 
            numberPointsAtBoundary =  (j_size + 1) * (k_size + 1) 
        elif (face == "bottom" or face == "top"):
            numberCellsAtBoundary =  i_size * j_size 
            numberPointsAtBoundary =  (i_size + 1) * (j_size + 1) 
        elif (face == "north" or face == "south"):
            numberCellsAtBoundary =  i_size * k_size
            numberPointsAtBoundary =  (i_size + 1) * (k_size + 1) 
    
        
        parentData = CGNS.intArray(numberCellsAtBoundary * 4)
    
        array_index = 0
        parentIndex = 0
        if (face == "west" or face == "east"):
    
            i_b = ilo 
            for k in range(k_size):
                for j in range(j_size):
                    j_b = jlo + j * j_sign 
                    k_b = klo + k * k_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + (k_b + 1) * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + (k_b + 1) * n_ij + (j_b + 1) * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + (j_b + 1) * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    parentData[parentIndex] = startBlockCell[block_index] + k_b * (n_ij -1) + j_b * (n_i-1) + i_b + 1
                    parentData[numberCellsAtBoundary + parentIndex ] = 0
                        
                    if (face == "west"):
                        parentData[2 * numberCellsAtBoundary + parentIndex] = 1
                    else: 
                        parentData[2 * numberCellsAtBoundary + parentIndex] = 2
    
                    parentData[3 * numberCellsAtBoundary + parentIndex] = 0
                    parentIndex = parentIndex + 1
    
        elif (face == "bottom" or face == "top"):
        
            k_b = klo 
            for j in range(j_size):
                for i in range(i_size):
                    i_b = ilo + i * i_sign 
                    j_b = jlo + j * j_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + (j_b + 1) * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + (j_b + 1) * n_i + (i_b + 1) 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + (i_b + 1) 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
    
        elif (face == "north" or face == "south"):
        
            j_b = jlo 
            for k in range(k_size):
                for i in range(i_size):
                    i_b = ilo + i * i_sign 
                    k_b = klo + k * k_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b + 1 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + (k_b + 1) * n_ij + j_b * n_i + (i_b + 1) 
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
                    structPoint = startBlockVertex[block_index] + (k_b + 1) * n_ij + j_b * n_i + i_b  
                    unstrPoint = inverseMap[structPoint]
                    Elements[array_index] = unstrPoint + 1
                    array_index = array_index + 1
    
    
    
        start = element_index + 1 
        end = start + numberCellsAtBoundary - 1
    
        element_index = end
        
        ier = CGNS.cg_section_write(fileValue, bValue, zoneValue, bndName, CGNS.QUAD_4, start, end, 0, Elements, sectionPointer)
        checkCGNSerror(ier)
#        sectionValue = sectionPointer.value() 
        
        for i in range(numberCellsAtBoundary):
            Elements[4*i] = start + i
        
    #    CGNS.cg_parent_data_write(fileValue, bValue, zoneValue, sectionValue, parentData)
    
        parentData.__del__
        
    #    write boundary condition
        boundaryPoints = CGNS.intArray(numberPointsAtBoundary)
    
        if (face == "west" or face == "east"):
    
            array_index = 0
            i_b = ilo 
            for k in range(k_size+1):
                for j in range(j_size+1):
                    j_b = jlo + j * j_sign 
                    k_b = klo + k * k_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    boundaryPoints[array_index] = unstrPoint + 1
    #                print arrayIndex, j_b, k_b, structPoint, unstrPoint
                    array_index = array_index + 1
    
        elif (face == "bottom" or face == "top"):
    
            array_index = 0
            k_b = klo 
            for j in range(j_size+1):
                for i in range(i_size+1):
                    j_b = jlo + j * j_sign 
                    i_b = ilo + i * i_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    boundaryPoints[array_index] = unstrPoint + 1
    #                print arrayIndex, j_b, k_b, structPoint, unstrPoint
                    array_index = array_index + 1
    
        elif (face == "north" or face == "south"):
    
            array_index = 0
            j_b = jlo 
            for k in range(k_size+1):
                for i in range(i_size+1):
                    k_b = klo + k * k_sign 
                    i_b = ilo + i * i_sign 
    
                    structPoint = startBlockVertex[block_index] + k_b * n_ij + j_b * n_i + i_b 
                    unstrPoint = inverseMap[structPoint]
                    boundaryPoints[array_index] = unstrPoint + 1
    #                print arrayIndex, j_b, k_b, structPoint, unstrPoint
                    array_index = array_index + 1
      
        
        boundaryPointer = CGNS.intp()
        bocotype = getCGNSboundaryType(bndName)
        ier = CGNS.cg_boco_write(fileValue, bValue, zoneValue, bndName, bocotype, CGNS.PointList, numberPointsAtBoundary, boundaryPoints, boundaryPointer)
    
        checkCGNSerror(ier)
        
        print "Boundary", bnd, bndName, "block", block, "cells", numberCellsAtBoundary
    
########################################################################################################    


# Start Main Program

# input data

shortOptions = "s:e:"
longOptions = []

if __name__ == '__main__':
    print "Begin Eilmer3 to cgns converter ..."
    userOptions, args = getopt(sys.argv[1:], shortOptions, longOptions)
    print args, userOptions
    if len(args) != 2:
        printUsage()
        sys.exit()

jobName = args[0]
print "Read grid from job: ", jobName
cgnsFileName = args[1]
print "Write grid to file: ", cgnsFileName

configFileName = jobName + ".config"

blockLabels = getBlockLabels(configFileName) 

domains = set(blockLabels)
numberDomains = len(domains)

print "Number Domains: ", numberDomains, domains


#if (len(userOptions) == 2):
#    
#    for opt, arg in userOptions:
#        if opt in ('-s'):
#            block_start = int(arg)
#        elif opt in ('-e'):
#            block_end = int(arg)

#else:
#    block_start = 0
#    block_end = numberBlocks
 
    
#jobName = sys.argv[1]
#cgnsFileName = sys.argv[2]
#block_start = int(sys.argv[3])
#block_end = int(sys.argv[4])




# open cgns file

filePointer = CGNS.intp()
ier = CGNS.cg_open(cgnsFileName, CGNS.CG_MODE_WRITE, filePointer)
fileValue = filePointer.value()

#Definition of globals

cell_dim = 3
phys_dim = 3
numberZones = 1
sizeZones = cell_dim*phys_dim

# write base data class + units


b = CGNS.intp()
ier = CGNS.cg_base_write ( fileValue, "Base", cell_dim, phys_dim, b )
bValue = b.value()

ier = CGNS.cg_goto ( fileValue, bValue, "end" )
ier = CGNS.cg_state_write ( "Dimensional" )
ier = CGNS.cg_dataclass_write ( CGNS.Dimensional )
ier = CGNS.cg_units_write ( CGNS.Kilogram, CGNS.Meter, CGNS.Second, CGNS.Kelvin, CGNS.Radian )


for domain in domains:
    print domain

    blocks = []
    for block_index in range(len(blockLabels)):
        if(blockLabels[block_index] == domain):
            blocks.append(block_index)


    boundaryConditions, blockConnections = readConfigFile(configFileName, blocks)

    writeZone(fileValue, domain, blocks, boundaryConditions, blockConnections)

ier = CGNS.cg_close ( filePointer.value() )
