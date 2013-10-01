## \file makefile
## \ingroup psr
## \brief Makefile for the PostShockRelax stuff
## \author RJG, DFP
## \date 29-Mar-2010
## \version 29-Mar-2010
##

TARGET ?= for_gnu
OPT ?= -O2
include ../../../lib/util/source/systems.mk

SRC             := ../source
INSTALL_DIR     := $(HOME)/e3bin
GAS             := ../../../lib/gas
NM              := ../../../lib/nm
UTIL_SRC        := ../../../lib/util/source
UTIL            := ../../../lib/util
RAD             := ../../../lib/radiation
RAD_SRC         := ../../../lib/radiation/source
GEOM            := ../../../lib/geometry2
GEOM_SRC        := ../../../lib/geometry2/source

LIBPSR          := libpsr.a
LIBGAS          := $(GAS)/build/libgas.a
LIBNM           := $(NM)/build/libnm.a
LIBUTIL         := $(UTIL)/build/libutil.a
LIBRAD         := $(RAD)/build/librad.a
LIBGEOM         := $(GEOM)/build/libgeom2.a
LIBINIPARSER    := $(UTIL)/iniparser/libiniparser.a

LUA             := ../../../extern/lua-5.1.4
LIBLUA          := $(LUA)/lib/liblua.a
LUA_INCLUDE_DIR := $(LUA)/include
LUALINK         := -ldl
ZLIB            := ../../../extern/zlib-1.2.3
LIBZLIB         := $(ZLIB)/libz.a

PSR_OBJECTS := conservation_systems.o \
	flow_state.o \
	post_shock_flow.o \
	poshax_radiation_transport.o

PYTHON_SCRIPTS := p3post.py

PSR_EXECUTABLE := poshax3.x

all: install
	@echo "--------------------------------"
	@echo "Poshax3 should now be up to date."
	@echo "--------------------------------"

test : $(TEST_EXE_FILES)

clean :
	- rm *.x *.o *.a *.py
	- cd $(RAD)/build; make clean

install : $(PSR_EXECUTABLE) $(PYTHON_SCRIPTS) $(INSTALL_DIR)
	- cp -f $(PSR_EXECUTABLE) $(PYTHON_SCRIPTS) $(INSTALL_DIR)
	cd $(GAS)/build; make install
	cd $(RAD)/build; make install

$(INSTALL_DIR) :
	mkdir $(INSTALL_DIR)

#----------- External libraries
$(LIBGAS) :
	cd $(GAS)/build; make libgas.a
#	cp $(GAS)/build/libgas.a .

$(LIBNM) :
	cd $(NM)/build; make libnm.a
#	cp $(NM)/build/libnm.a .

$(LIBUTIL) :
	cd $(UTIL)/build; make TARGET=$(TARGET) OPT=$(OPT) libutil.a
#	cp $(UTIL)/build/libutil.a .

$(LIBRAD) :
	cd $(RAD)/build; make librad.a
#	cp $(RAD)/build/librad.a .

$(LIBINIPARSER):
	cd $(UTIL)/iniparser; make

$(LIBGEOM) :
	cd $(GEOM)/build; make TARGET=$(TARGET) OPT=$(OPT) libgeom2.a

# Assume linux
PLAT := linux
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    PLAT := macosx
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    PLAT := cygwin
endif
# Even if we are using the Mingw compiler, assume that we are on cygwin.
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    PLAT := cygwin
endif

$(LIBLUA):
	cd $(LUA); make $(PLAT) local

$(LUA_INCLUDE_DIR): $(LIBLUA)

$(LIBZLIB):
	cd $(ZLIB); make libz.a

#----------- Python module (using swig) -------

$(LIBPSR) : $(PSR_OBJECTS)
	- rm $(LIBPSR)
	$(AR) rv $(LIBPSR) $(PSR_OBJECTS)
	- ranlib $(LIBPSR)

#----------- Python libraries ------------

p3post.py : $(SRC)/p3post.py
	cp $(SRC)/p3post.py .

#----------- executables -----------------
poshax3.x : poshax.o $(LIBPSR) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) \
                     $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(FLINK) -o poshax3.x poshax.o \
		$(LIBPSR) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB) $(LIBGEOM) $(LUALINK)

#----------- objects ---------
poshax.o : $(SRC)/poshax.cxx $(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/poshax.cxx -I$(LUA_INCLUDE_DIR)

conservation_systems.o : $(SRC)/conservation_systems.cxx \
                         $(SRC)/conservation_systems.hh $(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/conservation_systems.cxx -I$(LUA_INCLUDE_DIR)

flow_state.o : $(SRC)/flow_state.cxx $(SRC)/flow_state.hh $(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/flow_state.cxx -I$(LUA_INCLUDE_DIR)

post_shock_flow.o : $(SRC)/post_shock_flow.cxx $(SRC)/post_shock_flow.hh \
                    $(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/post_shock_flow.cxx -I$(LUA_INCLUDE_DIR)

poshax_radiation_transport.o : $(SRC)/poshax_radiation_transport.cxx \
		$(SRC)/poshax_radiation_transport.hh $(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/poshax_radiation_transport.cxx -I$(LUA_INCLUDE_DIR)

