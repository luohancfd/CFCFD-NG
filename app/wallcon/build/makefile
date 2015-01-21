# makefile for wallcon library and executable
#
# Author: Justin Beri and Rowan Gollan
# Date: 09-May-2014

TARGET ?= for_gnu
OPT ?= -O2
INSTALL_DIR ?= $(HOME)/e3bin
include ../../../lib/util/source/systems.mk

#----------------------------------
# Directory variables

# Local to this area
SRC             := ../source

# Other useful directory variables
NM_SRC          := ../../../lib/nm/source
NM              := ../../../lib/nm
LIBNM           := $(NM)/build/libnm.a
GEOM2           := ../../../lib/geometry2
GEOM2_SRC       := $(GEOM2)/source
LIBGEOM2        := $(GEOM2)/build/libgeom2.a
LUA             := ../../../extern/lua-5.1.4
LUA_INCLUDE_DIR := $(LUA)/include

#---------------------------------------
# File name variables

WALLCON_EXE := wallcon

LIBWALLCON := libwallcon.a

WC_OBJECTS = bc_adiabatic.o \
	bc_convection.o \
	bc_conv_rad.o\
	bc_e3connection.o \
	bc_fixed_q.o \
	bc_fixed_t.o \
	bc_radiation.o \
	bc_timevarying_q.o \
	bc_timevarying.o \
	bc_userdef_q.o \
	bc_userdef_t.o \
	e3conn.o \
	initialise.o \
	output.o \
	solid_block.o \
	solid_cell.o \
	space_update.o \
	useful_funcs.o

build_new : $(LIBWALLCON)

install : $(INSTALL_DIR) stand_alone
	@echo Copying wallcon executable to $(INSTALL_DIR)
	cp -f $(WALLCON_EXE) $(INSTALL_DIR)

stand_alone: $(LIBWALLCON) $(LIBNM) $(LIBGEOM2) $(SRC)/main.cxx
	$(CXXLINK) -std=c++0x $(LFLAG) $(SRC)/main.cxx $(LIBWALLCON) $(LIBNM) $(LIBGEOM2) -o $(WALLCON_EXE)

clean :
	- rm *.o *.a $(WALLCON_EXE)

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)
#----------------------------------------
# Building the static library

$(LIBWALLCON) : $(WC_OBJECTS)
	$(AR) rv $(LIBWALLCON) $(WC_OBJECTS)

bc_adiabatic.o : $(SRC)/bc_adiabatic.hh $(SRC)/bc_adiabatic.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_adiabatic.cxx

bc_convection.o : $(SRC)/bc_convection.hh $(SRC)/bc_convection.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_convection.cxx

bc_conv_rad.o : $(SRC)/bc_conv_rad.hh $(SRC)/bc_conv_rad.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_conv_rad.cxx

bc_e3connection.o : $(SRC)/bc_e3connection.hh $(SRC)/bc_e3connection.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_e3connection.cxx

bc_fixed_q.o : $(SRC)/bc_fixed_q.hh $(SRC)/bc_fixed_q.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_fixed_q.cxx

bc_fixed_t.o : $(SRC)/bc_fixed_t.hh $(SRC)/bc_fixed_t.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_fixed_t.cxx

bc_radiation.o : $(SRC)/bc_radiation.hh $(SRC)/bc_radiation.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_radiation.cxx

bc_timevarying_q.o : $(SRC)/bc_timevarying_q.hh $(SRC)/bc_timevarying_q.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_timevarying_q.cxx

bc_timevarying.o : $(SRC)/bc_timevarying.hh $(SRC)/bc_timevarying.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_timevarying.cxx

bc_userdef_q.o : $(SRC)/bc_userdef_q.hh $(SRC)/bc_userdef_q.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_userdef_q.cxx

bc_userdef_t.o : $(SRC)/bc_userdef_t.hh $(SRC)/bc_userdef_t.cxx 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bc_userdef_t.cxx

e3conn.o : $(SRC)/e3conn.hh $(SRC)/e3conn.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/e3conn.cxx -I$(LUA_INCLUDE_DIR)

initialise.o : $(SRC)/initialise.hh $(SRC)/initialise.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/initialise.cxx

output.o : $(SRC)/output.hh $(SRC)/output.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/output.cxx

solid_block.o : $(SRC)/solid_block.hh $(SRC)/solid_block.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/solid_block.cxx

solid_cell.o : $(SRC)/solid_cell.hh $(SRC)/solid_cell.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/solid_cell.cxx

space_update.o : $(SRC)/solid_block.hh $(SRC)/space_update.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/space_update.cxx

useful_funcs.o : $(SRC)/useful_funcs.hh $(SRC)/useful_funcs.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/useful_funcs.cxx

#--------------------- external dependencies -------------------------
$(LIBNM) :
	cd $(NM)/build; make TARGET=$(TARGET) OPT=$(OPT) libnm.a

$(LIBGEOM2):
	cd $(GEOM2)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new
