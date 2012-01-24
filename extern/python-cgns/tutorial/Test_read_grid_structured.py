#!/usr/bin/env python
# -*- coding: utf-8 -*-

#***************************************************************************
#*   Copyright (C) 2011 by Lionel Gamet                                    *
#*   Lionel.Gamet@fluorem.com                                              *
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

# Author: Lionel Gamet
# Date: January 2011

import sys
import numpy
import CGNS

# Initializations
#-----------------
# Change Flow'Design file name here
# The geometry "boite" is that of a cube with (idim,jdim,kdim) = (5,5,5).
# "boite_CGNSADF.pgeo"  has been created with Turb'Flow(R) compiled against CGNS-ADF library.
# "boite_CGNSHDF5.pgeo" has been created with Turb'Flow(R) compiled against CGNS-HDF5 library.
# Choose appropriate name depending on Low Level Library choice (HDF5 or ADF) :
FileName = "boite_CGNSADF.pgeo"
#FileName = "boite_CGNSHDF5.pgeo"

# Debug level
DebugLevel=1

# Index transformation function to handle 2D meshes
#---------------------------------------------------
def indextransform(transform,transpose,i1,j1,k1):    
  matrix=[[0,0,0],[0,0,0],[0,0,0]]
  if transpose == 0:
    matrix[0][transform[0]-1] = 1
    matrix[1][transform[1]-1] = 1
    matrix[2][transform[2]-1] = 1
  elif transpose == 1:
    matrix[transform[0]-1][0] = 1
    matrix[transform[1]-1][1] = 1
    matrix[transform[2]-1][2] = 1
  else :
    print "@@@indextransform: Incorrect second argument: Expecting 1(TRUE) or 0(FALSE)"
    return (0,0,0)
  i2 = matrix[0][0]*i1 + matrix[0][1]*j1 + matrix[0][2]*k1
  j2 = matrix[1][0]*i1 + matrix[1][1]*j1 + matrix[1][2]*k1
  k2 = matrix[2][0]*i1 + matrix[2][1]*j1 + matrix[2][2]*k1
  return (i2,j2,k2)


class ERREUR: pass
try:

  # Open Turb'Flow CGNS file for read-only
  #----------------------------------------
  ptr_F_index = CGNS.intp()
  ier = CGNS.cg_open(FileName,CGNS.CG_MODE_READ,ptr_F_index)
  if ier != 0:
    print "@@@ERROR: Cannot open file: %s" %(FileName)
    raise ERREUR()
  F_index = ptr_F_index.value()
  print "Opening file",FileName

  # Access to node CGNSBase_t
  index_base = 1
  index_zone = 1
  ier = CGNS.cg_goto(F_index,index_base,"end")
  if ier != 0: raise ERREUR()

  # Read identifier character string
  ier,Identname,Ident = CGNS.cg_descriptor_read(1)
  if ier != 0: raise ERREUR()
  if Ident != "CGNM3S1":
    print "@@@ERROR: Incorrect identifier for a TurbFlow CGNS file: %s" %(Ident)
    raise ERREUR()

  # Read dimensions
  #-----------------
  # Read mesh physical dimensions
  ier = CGNS.cg_goto(F_index,index_base,"end")
  if ier != 0: raise ERREUR()
  ptr_celldim  = CGNS.intp()
  ptr_physdim  = CGNS.intp()
  ier,basename = CGNS.cg_base_read(F_index,index_base,ptr_celldim,ptr_physdim)
  if ier != 0: raise ERREUR()
  celldim = ptr_celldim.value()
  physdim = ptr_physdim.value()

  # Access to node Zone_t
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,"end")
  if ier != 0: raise ERREUR()

  # Read mesh dimensions (number of nodes, number of cells, number of boundary nodes)
  size = CGNS.intArray(9)
  ier,zonename = CGNS.cg_zone_read(F_index,index_base,index_zone,size)
  if ier != 0: raise ERREUR()
  # Set number of nodes
  iDimL=1
  jDimL=1
  kDimL=1
  for n in range(physdim):
    if n==0: iDimL=size[n]
    if n==1: jDimL=size[n]
    if n==2: kDimL=size[n]
  (iDim,jDim,kDim) = (0,0,0)

  # Identifies the number of user data
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,"end")
  if ier != 0: raise ERREUR()
  ptr_nuserdata = CGNS.intp()
  ier = CGNS.cg_nuser_data(ptr_nuserdata)
  if ier != 0: raise ERREUR()
  nuserdata = ptr_nuserdata.value()
  index_FluoremData = 1
  ier,nodename = CGNS.cg_user_data_read(index_FluoremData)
  if nodename != "FluoremData":
    print "@@@ERROR: Could not find node FluoremData"
    raise ERREUR()

  # Iblank nodes detection does not work
  iblankedL=0
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,"ZoneGridConnectivity_t",1,"OversetHoles_t",1,"end")
  if ier == 0: iblankedL=1

  print "Done reading dimensions"

  # Reading user defined arrays information
  #-----------------------------------------
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,
                     "UserDefinedData_t",index_FluoremData,"end")
  if ier != 0: raise ERREUR()
  ptr_narrays = CGNS.intp()
  ier = CGNS.cg_narrays(ptr_narrays)
  if ier != 0: raise ERREUR()
  narrays = ptr_narrays.value()

  ptr_dtype = CGNS.DataType_tp()
  ptr_DataDimension = CGNS.intp()
  DimensionVector = CGNS.intArray(3)
  for n in range(1,narrays+1):
    ier,nodename = CGNS.cg_array_info(n,ptr_dtype.cast(),ptr_DataDimension,DimensionVector)
    if ier != 0: raise ERREUR()
    # Read index permutation array
    if nodename == "IndexTransform":
      MIndexTransform = CGNS.intArray(3)
      ier = CGNS.cg_array_read(n,MIndexTransform)
      if ier != 0: raise ERREUR()
      (iDim,jDim,kDim) = indextransform(MIndexTransform,1,iDimL,jDimL,kDimL)
      print "Structured mesh has dimensions",iDim,"x",jDim,"x",kDim
    # Read boundary periodicity vector
    if nodename == "BCPeriodicityVector":
      PasNCGN = CGNS.doubleArray(9)
      ier = CGNS.cg_array_read_as(n,CGNS.RealDouble,PasNCGN)
      if ier != 0: raise ERREUR()
    # Read boundary normal vector
    if nodename == "BCNormalVector":
      NorNCGN = CGNS.doubleArray(18)
      ier = CGNS.cg_array_read_as(n,CGNS.RealDouble,NorNCGN)
      if ier != 0: raise ERREUR()

  print "Done reading user defined arrays"

  # Reading temporal information
  #------------------------------
  ier = CGNS.cg_goto(F_index,index_base,"BaseIterativeData_t",1,"end")
  if ier != 0: raise ERREUR()
  ier = CGNS.cg_narrays(ptr_narrays)
  if ier != 0: raise ERREUR()
  narrays = ptr_narrays.value()

  for n in range(1,narrays+1):
    ier,nodename = CGNS.cg_array_info(n,ptr_dtype.cast(),ptr_DataDimension,DimensionVector)
    if ier != 0: raise ERREUR()
    # Read time value
    if nodename == "TimeValues":
      tmparray = CGNS.doubleArray(1)
      ier = CGNS.cg_array_read_as(n,CGNS.RealDouble,tmparray)
      if ier != 0: raise ERREUR()
      time = tmparray[0]
    # Read number of iterations per cycle
    if nodename == "IterationPerCycleValues":
      tmparray = CGNS.intArray(1)
      ier = CGNS.cg_array_read_as(n,CGNS.Integer,tmparray)
      if ier != 0: raise ERREUR()
      nbritcyc = tmparray[0]
  
  # Read number of iterations
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,"end")
  if ier != 0: raise ERREUR()
  ptr_nbrite = CGNS.intp()
  ier,tmpstring = CGNS.cg_convergence_read(ptr_nbrite)
  if ier != 0: raise ERREUR()
  nbrite = ptr_nbrite.value()
  print "Done reading temporal information"

  # Reading title and comment
  #---------------------------
  ier = CGNS.cg_goto(F_index,index_base,"Zone_t",index_zone,
                     "UserDefinedData_t",index_FluoremData,"end")
  if ier != 0: raise ERREUR()
  ptr_ndescrs = CGNS.intp()
  ier = CGNS.cg_ndescriptors(ptr_ndescrs)
  if ier != 0: raise ERREUR()
  ndescrs = ptr_ndescrs.value()

  for n in range(1,ndescrs+1):
    ier,nodename,tmpstring = CGNS.cg_descriptor_read(n)
    if ier != 0: raise ERREUR()
    # Read title
    if nodename == "Title":
      Title = tmpstring
      print "  Title:",Title
    # Read comment
    if nodename == "Comment":
      Comment = tmpstring
      print "  Comment:",Comment
    

  # Reading of mesh coordinates
  #-----------------------------
  # Reading coordinates as single dimension arrays
  xcgn = CGNS.doubleArray(iDimL*jDimL*kDimL)
  ycgn = CGNS.doubleArray(iDimL*jDimL*kDimL)
  zcgn = CGNS.doubleArray(iDimL*jDimL*kDimL)

  ptr_dtype = CGNS.DataType_tp()
  ier,coordxname = CGNS.cg_coord_info(F_index,index_base,index_zone,1,ptr_dtype.cast())
  if ier != 0: raise ERREUR()
  ier,coordyname = CGNS.cg_coord_info(F_index,index_base,index_zone,2,ptr_dtype.cast())
  if ier != 0: raise ERREUR()
  ier,coordzname = CGNS.cg_coord_info(F_index,index_base,index_zone,3,ptr_dtype.cast())
  if ier != 0: raise ERREUR()
  
  rmin = CGNS.intArray(3)
  rmin[0] = 1
  rmin[1] = 1
  rmin[2] = 1
  rmax = CGNS.intArray(3)
  rmax[0] = iDimL
  rmax[1] = jDimL
  rmax[2] = kDimL
  ier = CGNS.cg_coord_read(F_index,index_base,index_zone,coordxname,
                           CGNS.RealDouble,rmin,rmax,xcgn)
  if ier != 0: raise ERREUR()
  ier = CGNS.cg_coord_read(F_index,index_base,index_zone,coordyname,
                           CGNS.RealDouble,rmin,rmax,ycgn)
  if ier != 0: raise ERREUR()
  ier = CGNS.cg_coord_read(F_index,index_base,index_zone,coordzname,
                           CGNS.RealDouble,rmin,rmax,zcgn)
  if ier != 0: raise ERREUR()

  # Store coordinates in 3 dimensions numpy arrays
  x = numpy.zeros((iDim,jDim,kDim))
  y = numpy.zeros((iDim,jDim,kDim))
  z = numpy.zeros((iDim,jDim,kDim))
  PasI = numpy.zeros((3))
  PasJ = numpy.zeros((3))
  PasK = numpy.zeros((3))
  NorI = numpy.zeros((3,2))
  NorJ = numpy.zeros((3,2))
  NorK = numpy.zeros((3,2))

  #----
  # 3D
  #----
  if physdim == 3:
    for k in range(0,kDim):
      for j in range(0,jDim):
        for i in range(0,iDim):
          n=i+j*iDimL+k*iDimL*jDimL
          x[i][j][k] = xcgn[n]
          y[i][j][k] = ycgn[n]
          z[i][j][k] = zcgn[n]
    for i in range(0,3):
      PasI[i] = PasNCGN[i+0*3]
      PasJ[i] = PasNCGN[i+1*3]
      PasK[i] = PasNCGN[i+2*3]
      for j in range(0,2):
        NorI[i][j] = NorNCGN[i+0*6+j*3]
        NorJ[i][j] = NorNCGN[i+1*6+j*3]
        NorK[i][j] = NorNCGN[i+2*6+j*3]

  #----
  # 2D
  #----
  if physdim == 2:

    # I<->K
    #-------
    if iDim == 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=k+j*iDimL+i*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasK[i] = PasNCGN[i+0*3]
        PasJ[i] = PasNCGN[i+1*3]
        PasI[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorK[i][j] = NorNCGN[i+0*6+j*3]
          NorJ[i][j] = NorNCGN[i+1*6+j*3]
          NorI[i][j] = NorNCGN[i+2*6+j*3]

    # J<->K
    #-------
    if jDim == 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=i+k*iDimL+j*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasI[i] = PasNCGN[i+0*3]
        PasK[i] = PasNCGN[i+1*3]
        PasJ[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorI[i][j] = NorNCGN[i+0*6+j*3]
          NorK[i][j] = NorNCGN[i+1*6+j*3]
          NorJ[i][j] = NorNCGN[i+2*6+j*3]

    # Nothing
    #---------
    if kDim == 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=i+j*iDimL+k*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasI[i] = PasNCGN[i+0*3]
        PasJ[i] = PasNCGN[i+1*3]
        PasK[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorI[i][j] = NorNCGN[i+0*6+j*3]
          NorJ[i][j] = NorNCGN[i+1*6+j*3]
          NorK[i][j] = NorNCGN[i+2*6+j*3]
  #----
  # 1D
  #----
  if physdim == 1:

    # Nothing
    #---------
    if iDim != 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=i+j*iDimL+k*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasI[i] = PasNCGN[i+0*3]
        PasJ[i] = PasNCGN[i+1*3]
        PasK[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorI[i][j] = NorNCGN[i+0*6+j*3]
          NorJ[i][j] = NorNCGN[i+1*6+j*3]
          NorK[i][j] = NorNCGN[i+2*6+j*3]

    # J<->I
    #-------
    if jDim != 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=j+i*iDimL+k*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasJ[i] = PasNCGN[i+0*3]
        PasI[i] = PasNCGN[i+1*3]
        PasK[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorJ[i][j] = NorNCGN[i+0*6+j*3]
          NorI[i][j] = NorNCGN[i+1*6+j*3]
          NorK[i][j] = NorNCGN[i+2*6+j*3]

    # K<->I
    #-------
    if kDim != 1:
      for k in range(0,kDim):
        for j in range(0,jDim):
          for i in range(0,iDim):
            n=j+i*iDimL+k*iDimL*jDimL
            x[i][j][k] = xcgn[n]
            y[i][j][k] = ycgn[n]
            z[i][j][k] = zcgn[n]
      for i in range(0,3):
        PasK[i] = PasNCGN[i+0*3]
        PasJ[i] = PasNCGN[i+1*3]
        PasI[i] = PasNCGN[i+2*3]
        for j in range(0,2):
          NorK[i][j] = NorNCGN[i+0*6+j*3]
          NorJ[i][j] = NorNCGN[i+1*6+j*3]
          NorI[i][j] = NorNCGN[i+2*6+j*3]

  print "Done reading coordinates"


  # Debugging output (equivalent to SPIDER Flow'Design ASCII format)
  #------------------
  if DebugLevel >= 1:
    for k in range(0,kDim):
      for j in range(0,jDim):
        for i in range(0,iDim):
          xstmp="%22.15E" %x[i][j][k]
          ystmp="%22.15E" %y[i][j][k]
          zstmp="%22.15E" %z[i][j][k]
          print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorI[0][0]
    ystmp="%22.15E" %NorI[1][0]
    zstmp="%22.15E" %NorI[2][0]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorI[0][1]
    ystmp="%22.15E" %NorI[1][1]
    zstmp="%22.15E" %NorI[2][1]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorJ[0][0]
    ystmp="%22.15E" %NorJ[1][0]
    zstmp="%22.15E" %NorJ[2][0]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorJ[0][1]
    ystmp="%22.15E" %NorJ[1][1]
    zstmp="%22.15E" %NorJ[2][1]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorK[0][0]
    ystmp="%22.15E" %NorK[1][0]
    zstmp="%22.15E" %NorK[2][0]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %NorK[0][1]
    ystmp="%22.15E" %NorK[1][1]
    zstmp="%22.15E" %NorK[2][1]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %PasI[0]
    ystmp="%22.15E" %PasI[1]
    zstmp="%22.15E" %PasI[2]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %PasJ[0]
    ystmp="%22.15E" %PasJ[1]
    zstmp="%22.15E" %PasJ[2]
    print xstmp,ystmp,zstmp
    xstmp="%22.15E" %PasK[0]
    ystmp="%22.15E" %PasK[1]
    zstmp="%22.15E" %PasK[2]
    print xstmp,ystmp,zstmp

  # Close file and exit
  #---------------------
  ier = CGNS.cg_close(F_index)
  if ier != 0: raise ERREUR()
  sys.exit(0)

# Error handling
#----------------
except ERREUR:
  print "ERROR: Error reading file <",FileName,">"
  ier = CGNS.cg_close(F_index)
  sys.exit(1)
