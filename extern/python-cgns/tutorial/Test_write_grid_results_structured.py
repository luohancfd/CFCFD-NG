#!/usr/bin/env python
# -*- coding: utf-8 -*-

#***************************************************************************
#*   Copyright (C) 2011 by Oliver Borm                                     *
#*  oli.borm@web.de                                                        *
#*                                                                         *
#*   This program is free software; you can redistribute it and/or modify  *
#*   it under the terms of the GNU General Public License as published by  *
#*   the Free Software Foundation; either version 3 of the License, or     *
#*   (at your option) any later version.                                   *
#*                                                                         *
#*   This program is distributed in the hope that it will be useful,       *
#*   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
#*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
#*   GNU General Public License for more details.                          *
#*                                                                         *
#*   You should have received a copy of the GNU General Public License     *
#*   along with this program; if not, write to the                         *
#*   Free Software Foundation, Inc.,                                       *
#*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
#***************************************************************************

# Author: Oliver Borm
# Date: January 2011

import CGNS

# Definition of globals
#-----------------------
cell_dim = 3
phys_dim = 3
numberZones = 1

NVertexR = 5
NVertexZ = 15
NVertexdepth = 3
numberNodes = NVertexZ*NVertexR*NVertexdepth

length = 2.0
width = 1.0
depth = 0.5

debug = 2

dateiName = "abc0.cgns"
blockName = "block_"

########################################
# Testing of CGNS file creation+writing
########################################

# Opening and creation of CGNSBase_t
#------------------------------------
sizeZones   = cell_dim*phys_dim

filePointer = CGNS.intp()
ier = CGNS.cg_open(dateiName, CGNS.CG_MODE_WRITE, filePointer)
fileValue = filePointer.value()
b = CGNS.intp()
ier = CGNS.cg_base_write ( fileValue, "Base", cell_dim, phys_dim, b )
ier = CGNS.cg_close ( filePointer.value() )

# Opening and placement at CGNSBase_t end
#-----------------------------------------
bPointer = CGNS.intp()
bPointer.assign(1)
bValue = bPointer.value()

ier = CGNS.cg_open ( dateiName, CGNS.CG_MODE_MODIFY, filePointer)
fileValue = filePointer.value()
ier = CGNS.cg_goto ( fileValue, bValue, "end" )

# Creates Subknot of Type 'Zone_t' (Chapter 6.2)
#------------------------------------------------
zoneArrayp= CGNS.intArray(sizeZones)
# Vertex size
zoneArrayp[0] = NVertexdepth #X-Coordinate
zoneArrayp[1] = NVertexR     #Y-Coordinate
zoneArrayp[2] = NVertexZ     #Z-Coordinate
# Cell size
zoneArrayp[3] = zoneArrayp[0]-1
zoneArrayp[4] = zoneArrayp[1]-1
zoneArrayp[5] = zoneArrayp[2]-1
# Boundary vertex size (always zero for structured grids)
zoneArrayp[6] = 0
zoneArrayp[7] = 0
zoneArrayp[8] = 0

if debug >=2:
  for i in range(sizeZones):
    print zoneArrayp[i]

# Initialization of CFD mesh and fields
#---------------------------------------
CoordX= CGNS.doubleArray(numberNodes)
CoordY= CGNS.doubleArray(numberNodes)
CoordZ= CGNS.doubleArray(numberNodes)

sol_rho= CGNS.doubleArray(numberNodes)
sol_pre= CGNS.doubleArray(numberNodes)
sol_tem= CGNS.doubleArray(numberNodes)
sol_ent= CGNS.doubleArray(numberNodes)
sol_vlx= CGNS.doubleArray(numberNodes)
sol_vly= CGNS.doubleArray(numberNodes)
sol_vlz= CGNS.doubleArray(numberNodes)
sol_csd= CGNS.doubleArray(numberNodes)

n = 0
for i in range(NVertexZ):
  for j in range(NVertexR):
    for k in range(NVertexdepth):
      if NVertexdepth > 1:
        CoordX[n] = depth/(NVertexdepth-1)*k
      else:
        CoordX[n] = 0.0
              
      CoordY[n] = width/(NVertexR-1)*j
      CoordZ[n] = length/(NVertexZ-1)*i
      
      sol_rho[n] = 1.0      +0.0125*n
      sol_pre[n] = 101300.0 +0.111*n
      sol_tem[n] = 293.0    +0.123*n
      sol_ent[n] = 10.0
      sol_vlx[n] = 40.0 +1.234*n
      sol_vly[n] = 10.0 +4.567*n
      sol_vlz[n] = 5.0  +8.901*n
      sol_csd[n] = 340.0
      
      n += 1

if debug >=3:
  for n in range(numberNodes):
    print CoordX[n], CoordY[n], CoordZ[n],sol_rho[n]

# Writes Coordinates (Chapter 11.1)
#-----------------------------------
zoneSizep = CGNS.intp()
ier = CGNS.cg_zone_write ( fileValue, bValue, blockName+str(bValue) , zoneArrayp, CGNS.Structured, zoneSizep )

zoneSizepValue = zoneSizep.value()

if debug >=2: print "Write Coordinates"

cIndex = CGNS.intp()

ier = CGNS.cg_coord_write ( fileValue, bValue, zoneSizepValue, CGNS.RealDouble, "CoordinateX", CoordX, cIndex )
if debug >=2: print ier, cIndex.value()

ier = CGNS.cg_coord_write ( fileValue, bValue, zoneSizepValue, CGNS.RealDouble, "CoordinateY", CoordY, cIndex )
if debug >=2: print ier, cIndex.value()

ier = CGNS.cg_coord_write ( fileValue, bValue, zoneSizepValue, CGNS.RealDouble, "CoordinateZ", CoordZ, cIndex )
if debug >=2: print ier, cIndex.value()


# Writes FlowSolution
#---------------------
if debug >=2: print "Write Solution"
sIndex = CGNS.intp()
ier = CGNS.cg_sol_write ( fileValue, bValue, zoneSizepValue, "FlowSolution", CGNS.Vertex, sIndex ) #// (Chapter 12.1)
if debug >=2 and ier!=0: print ier
sValue = sIndex.value()

#@@@ Does not work yet
#ier = CGNS.cg_gridlocation_write ( CGNS.Vertex ) # // (Chapter 9.1)
#if debug >=2 and ier!=0: print ier

ier = CGNS.cg_goto ( fileValue, bValue, "Zone_t", zoneSizepValue, "FlowSolution_t", sValue, "end" )
if debug >=2 and ier!=0: print ier


# SIDS-compliant FlowSolution // (Chapter 12.1)
if debug >=2: print "Write Solutionfields"
fIndex = CGNS.intp()
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "Density", sol_rho, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "Pressure", sol_pre, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "Temperature", sol_tem, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "Entropy", sol_ent, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "VelocityX", sol_vlx, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "VelocityY", sol_vly, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "VelocityZ", sol_vlz, fIndex )
if debug >=2: print ier
ier = CGNS.cg_field_write ( fileValue, bValue, zoneSizepValue, sValue, CGNS.RealDouble, "VelocitySound", sol_csd, fIndex )
if debug >=2: print ier

# Reference
ier = CGNS.cg_goto ( fileValue, bValue, "end" )
ier = CGNS.cg_state_write ( "Dimensional" ) # (Chapter 10.1)
ier = CGNS.cg_close ( fileValue )


#####################################
# Testing of CGNS file modification
#####################################
filePointer = CGNS.intp()
ier = CGNS.cg_open ( dateiName, CGNS.CG_MODE_MODIFY, filePointer)
fileValue = filePointer.value()

bPointer = CGNS.intp()
bPointer.assign(1)
bValue = bPointer.value()

ier = CGNS.cg_goto ( fileValue, bValue, "end" )
ier = CGNS.cg_dataclass_write ( CGNS.Dimensional )
ier = CGNS.cg_units_write ( CGNS.Kilogram, CGNS.Meter, CGNS.Second, CGNS.Kelvin, CGNS.Radian )

nzones = CGNS.intp()
nsolutions = CGNS.intp()
nfields = CGNS.intp()

exp = CGNS.doubleArray(5)
exp[0] = 1.0
exp[1] = -3.0
exp[2] = 0.0
exp[3] = 0.0
exp[4] = 0.0
if debug >=3:
  for i in range(5):
    print "exp[",i,"]=",exp[i]
                        
dataType = CGNS.DataType_tp()
dataType.assign(CGNS.RealDouble)

ier = CGNS.cg_nzones ( fileValue, bValue, nzones )

print "Number of zones :",nzones.value()
for Z in range(1,nzones.value()+1):
  
  ier = CGNS.cg_nsols ( fileValue, bValue, Z, nsolutions )
  print "  Number of solutions in zone Z=",Z,':',nsolutions.value()
  for S in range(1,nsolutions.value()+1):
            
    ier = CGNS.cg_nfields ( fileValue, bValue, Z, S, nfields )
    print "    Number of fields of solution S=",S,':',nfields.value()
    for F in range(1,nfields.value()+1):

      ier,dataFieldName = CGNS.cg_field_info ( fileValue, bValue, Z, S, F, dataType.cast() )
      if debug >=2: print "      Field name :",dataFieldName

      print "      fileValue, bValue, Z, S, F :",fileValue, bValue, Z, S, F
      
      ier = CGNS.cg_goto ( fileValue, bValue, "Zone_t", Z, "FlowSolution_t", S, "DataArray_t", F, "end" )
      ier = CGNS.cg_exponents_write ( CGNS.RealDouble, exp )
      #ier = CGNS.cg_expfull_write ( CGNS.RealDouble, exp )
      #if debug >=2 and ier!=0: print "        ",ier



# Closes CGNS file (Chapter 3.1)
ier = CGNS.cg_close ( fileValue )
