#!/usr/bin/env python

import ConfigParser
from getopt import getopt
import sys, gzip, os
import math

sys.path.append(os.path.expandvars("$HOME/e3bin"))
from e3_grid import StructuredGrid
from e3_flow import StructuredGridFlow
from libprep3 import *


def getTurboData(fileName):

    Config = ConfigParser.ConfigParser()
    Config.read(fileName)

    numberBlocks = int(Config.get("global_data", "nblock"))
    print "Number of Blocks:", numberBlocks
    dimensions = int(Config.get("global_data", "dimensions"))

    data = dict()

    data['numberBlocks'] = numberBlocks
    data['dimensions'] = dimensions
    data['inlet'] = []
    data['outlet'] = []

    if (dimensions == 2): 
        blockFaces = ['east', 'west', 'south', 'north']
    else:
        blockFaces = ['east', 'west', 'south', 'north', 'bottom', 'top']

    for block in range(numberBlocks):
        for face in blockFaces:
            section = "block/" + str(block) + "/face/" + face
            bnd_type = Config.get(section, "bc")
            isWall = Config.get(section, "is_wall")
            connectedBlock = int(Config.get(section, "other_block"))
            bclabel = Config.get(section, "label")
            if (bclabel == "INLET"):
                data['inlet'].append({'block':block, 'face':face})
            elif (bclabel == "OUTLET"):
                data['outlet'].append({'block':block, 'face':face})
  
    return data


def readFlowData(jobName, numberBlocks):

    tindx = '9999'
    zeros = "000"
    
    flow = []
    grid = []

    for block in range(numberBlocks):
    
        nZeros = 4 - len(str(block))

        bindx = zeros[0:nZeros] + str(block) 
    
    #  read blocks
    
        fileName = 'flow/t'+tindx+'/'+ jobName + '.flow.b'+bindx+'.t'+tindx+'.gz'
        fp = gzip.open(fileName, "r")
        flow.append(StructuredGridFlow())
        flow[block].read(fp)
        fp.close()
#        print "flow data: ni=", flow[block].ni, "nj=", flow[block].nj, "nk=", flow[block].nk
        # The grid is always at tindx 0.
        fileName = 'grid/t0000/'+ jobName + '.grid.b'+bindx+'.t0000.gz'
        fp = gzip.open(fileName, "r")
        grid.append(StructuredGrid())
        grid[block].read(fp)
        fp.close()
#        print "Block ", block, "grid data: ni=", grid[block].ni, "nj=", grid[block].nj, "nk=", grid[block].nk

    return flow, grid        

def getFaceArea(grid, face, i, j, k):

    p0,p1,p2,p3,p4,p5,p6,p7 = grid.get_vertex_list_for_cell(i,j,k)
    if (grid.nk == 1):  ##  2D grid
        if (face == "west"):
            c0 = p0
            c1 = p3
        elif (face == "east"):
            c0 = p1
            c1 = p2
        elif (face == "north"):
            c0 = p2
            c1 = p3
        elif (face == "south"):
            c0 = p0
            c1 = p1
        else:
            print "Add code here"
            sys.exit()
        surface_area = vabs(c0-c1)        
    else:   ## 3D grid
        if (face == "west"):
            c0 = p0
            c1 = p3
            c2 = p7
            c3 = p4
        elif (face == "east"):
            c0 = p1
            c1 = p2
            c2 = p6
            c3 = p5
        elif (face == "south"):
            c0 = p0
            c1 = p1
            c2 = p5
            c3 = p4
        elif (face == "north"):
            c0 = p2
            c1 = p3
            c2 = p7
            c3 = p6
        elif (face == "bottom"):
            c0 = p0
            c1 = p1
            c2 = p2
            c3 = p3
        elif (face == "top"):
            c0 = p4
            c1 = p5
            c2 = p6
            c3 = p7
        else:
            print "Add code here"
            sys.exit()
        surface_area = quad_area(c0, c1, c2, c3)

#    surface_centroid = quad_centroid(p0, p1, p2, p3)
#    surface_normal = quad_normal(p0, p1, p2, p3)

    return surface_area

def printAveragedFlowAtBoundary(boundary, flow, grid):

    averPres = 0.0
    averRho = 0.0
    averVel = 0.0
    totalArea = 0.0
    averSpeedSound = 0.0

    rho = 0.0
    pres = 0.0
    vel = 0.0

    for bc in boundary:
        block = bc['block']
        face = bc['face']

        if (bc['face'] == "west"):
            iRange = [0]
            jRange = range(flow[block].nj) 
            kRange = range(flow[block].nk) 
        elif (bc['face'] == "east"):
            iRange = [flow[block].ni-1]
            jRange = range(flow[block].nj) 
            kRange = range(flow[block].nk) 
        elif (bc['face'] == "south"):
            iRange = range(flow[block].ni)
            jRange = [0]
            kRange = range(flow[block].nk) 
        elif (bc['face'] == "north"):
            iRange = range(flow[block].ni)
            jRange = [flow[block].nj-1] 
            kRange = range(flow[block].nk) 
        elif (bc['face'] == "bottom"):
            iRange = range(flow[block].ni)
            jRange = range(flow[block].nj)
            kRange = [0]
        elif (bc['face'] == "top"):
            iRange = range(flow[block].ni)
            jRange = range(flow[block].nj)
            kRange = [flow[block].nk-1] 
        else:
            print "Add code here"
            sys.exit(-1)
#        print iRange, jRange, kRange
        for k in kRange:
            for j in jRange:
                for i in iRange:
                    rho = flow[block].data['rho'][i,j,k] 
                    velx = flow[block].data['vel.x'][i,j,k] 
                    vely = flow[block].data['vel.y'][i,j,k] 
                    velz = flow[block].data['vel.z'][i,j,k] 
                    pres = flow[block].data['p'][i,j,k] 
                    speedSound = flow[block].data['a'][i,j,k] 
                    faceArea = getFaceArea(grid[block], face, i, j, k)
 #                   print block, face, i, j, k, faceArea, pres
                    vel = math.sqrt(velx * velx + vely * vely + velz * velz)
                    averRho = averRho + rho * faceArea
                    averPres = averPres + rho * faceArea * pres
                    averVel = averVel + rho * faceArea * vel
                    averSpeedSound = averSpeedSound + rho * faceArea * speedSound
                    totalArea = totalArea + faceArea
        

    averPres = averPres / averRho
    averVel = averVel / averRho
    averSpeedSound = averSpeedSound / averRho
    averRho = averRho / totalArea

    averMach = averVel / averSpeedSound

    print "Area: ", totalArea
    print "Average pressure: ", averPres
    print "Average velocity: ", averVel
    print "Average Speed of Sound: ", averSpeedSound
    print "Average Mach number: ", averMach



#######################################################

# Start Main Program

# input data

jobName = sys.argv[1]
print "Read grid from job: ", jobName


configFileName = jobName + ".config"

turboData = getTurboData(configFileName) 
#print turboData
flow, grid = readFlowData(jobName, turboData['numberBlocks'])
print "INLET"
printAveragedFlowAtBoundary(turboData['inlet'], flow, grid)
print "OUTLET"
printAveragedFlowAtBoundary(turboData['outlet'], flow, grid)


#print turboData
