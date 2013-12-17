# makefile for twc library and executable
#
# Author: Jared Clifford and Rowan Gollan
# Date: 19-Nov-2013

TARGET ?= for_gnu
OPT ?= -O2
include ../../../lib/util/source/systems.mk

#----------------------------------
# Directory variables

# Local to this area
SRC             := ../source

# Other useful directory variables
UTIL_SRC        := ../../../util/source
UTIL            := ../../../util
NM              := ../../../nm

#---------------------------------------
# File name variables

LIBTWC := libtwc.a

TWC_OBJECTS = adv.o \
	conjugateHeatUtil.o \
	e3con.o \
	init.o

build_new : $(LIBTWC)

clean :
	- rm *.o *.a

#----------------------------------------
# Building the static library

$(LIBTWC) : $(TWC_OBJECTS)
	$(AR) rv $(LIBTWC) $(TWC_OBJECTS)

adv.o : $(SRC)/adv.hh $(SRC)/adv.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/adv.cxx

conjugateHeatUtil.o : $(SRC)/conjugateHeatUtil.hh $(SRC)/conjugateHeatUtil.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/conjugateHeatUtil.cxx

e3con.o : $(SRC)/e3con.hh $(SRC)/e3con.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/e3con.cxx

init.o : $(SRC)/init.hh $(SRC)/init.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/init.cxx



