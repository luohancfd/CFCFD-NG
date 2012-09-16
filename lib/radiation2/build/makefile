# makefile for radiation2 library
#
# Author: Daniel F. Potter
# Version: 02-July-2009

#----------------------------------------
# Some default compile options

TARGET ?= for_gnu
OPT ?= -O2
INSTALL_DIR ?= $(HOME)/e3bin
WITH_SPRADIAN ?= 0
SPRADIAN_SRC ?= ../../../extern/spradian07
include ../../util/source/systems.mk

OMPFLAG ?= 
ifeq ($(TARGET), for_gnu_openmp)
    OMPFLAG := -fopenmp
endif

ifeq ($(TARGET), for_intel_openmp)
    OMPFLAG := -openmp
endif

#----------------------------------
# Directory variables

# Local to this area
SRC             := ../source

# Other useful directory variables
GAS             := ../../gas
LIBGAS          := $(GAS)/build/libgas.a
GASPY           := $(GAS)/build/gaspy.py
NM              := ../../nm
LIBNM           := $(NM)/build/libnm.a
LUA             := ../../../extern/lua-5.1.4
LIBLUA          := $(LUA)/lib/liblua.a
LUA_INCLUDE_DIR := $(LUA)/include
LUALINK         := -ldl
UTIL_SRC        := ../../util/source
UTIL            := ../../util
GEOM            := ../../geometry2
GEOM_SRC        := ../../geometry2/source
LIBGEOM         := $(GEOM)/build/libgeom2.a
LIBUTIL         := $(UTIL)/build/libutil.a
ZLIB            := ../../../extern/zlib-1.2.3
LIBZLIB         := $(ZLIB)/libz.a
CEA2            := ../../../extern/cea2
CFPYLIB         := ../../cfpylib

ifeq ($(shell [ -d ../../../extern/cea2 ] && echo "true"), true)
  WITH_CEA2 := 1
else
  WITH_CEA2 := 0
endif

#---------------------------------------
# File name variables

LIBRAD2 := librad2.a

RAD_OBJECTS := spectral_model.o \
	equilibrium_air.o \
	photaura.o \
	lua_service.o \
	radiator.o \
	diatomic_radiator.o \
	diatomic_system.o \
	diatomic_band.o \
	atomic_radiator.o \
	atomic_line.o \
	electron_radiator.o \
	planck_radiator.o \
	spectra_pieces.o \
	LOS_pieces.o \
	photoionisation.o \
	cr_reactions.o \
	cr_rr_coeffs.o \
	parade.o \
	parade_radiator.o \
	mersenne.o
	
ifeq ($(WITH_SPRADIAN), 1)
    RAD_OBJECTS += spradian.o spradian_radiator.o radipac6.o
endif

TEST_EXE := photaura_test.x

ifeq ($(WITH_SPRADIAN), 1)
    TEST_EXE += spradian_test.x run_spradian.x
endif

TEST_PY := photaura_test.py radGUI.py

ifeq ($(WITH_CEA2), 1)
	TEST_EXE += cea2
endif

# Assume Unix/Linux
LOADABLE_MODULE := _librad2.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _librad2.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _librad2.dll
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _librad2.dylib
endif

#---------------------------------------
# main targets

build_new: $(LIBRAD2)
	@echo "------------------------------------------------"
	@echo "Radiation library (v2) should now be up to date."
	@echo "------------------------------------------------"

install : build_new $(LOADABLE_MODULE) librad2.py
	@echo "Install python pieces from the librad2 functions."
	- cp $(LOADABLE_MODULE) librad2.py $(INSTALL_DIR)
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
	- cp $(INSTALL_DIR)/_librad2.dylib $(INSTALL_DIR)/_librad2.so
endif
	- cp $(SRC)/radiator_library2.py $(SRC)/parade_radiator_library.py \
		$(SRC)/script_rad2.py $(SRC)/EQ-spectra.py $(INSTALL_DIR)
	
test: install $(TEST_EXE) $(TEST_PY) $(GASPY)
	cd ../../; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)
	- cp $(TEST_PY) $(INSTALL_DIR)
	- cp $(TEST_EXE) $(INSTALL_DIR)
ifeq ($(WITH_CEA2), 1)
	- cp $(CEA2)/thermo.inp $(CEA2)/trans.inp $(INSTALL_DIR)
	@echo "---------------------------------------------------"
	@echo "Installation of radiation test components complete."
	@echo "CEA2 has been included in the build."
	@echo "---------------------------------------------------"
else
	@echo "--------------------------------------------------------------------"
	@echo "Installation of radiation test components complete."
	@echo "The cea2_gas python module will NOT be functional"
	@echo "because the CEA2 program has not been found in $(CEA2)."
	@echo "Contact the NASA Glenn Research Center to request the software:"
	@echo "http://www.grc.nasa.gov/WWW/CEAWeb/"
	@echo "--------------------------------------------------------------------"
endif
	
clean :
	- rm *.o *.py $(LIBRAD2) $(LOADABLE_MODULE) *.x *.a *.mod $(TEST_EXE)
	- rm $(SRC)/librad2.py $(SRC)/librad2_wrap.cxx
	- cd $(GAS)/build; make clean
	
clean_local :
	- rm *.o *.py $(LIBRAD2) $(LOADABLE_MODULE) *.x *.a *.mod $(TEST_EXE)
	- rm $(SRC)/librad2.py $(SRC)/librad2_wrap.cxx
	
clean_swig :
	- rm librad2.py $(LOADABLE_MODULE) librad2_wrap.o
	- rm $(SRC)/librad2.py $(SRC)/librad2_wrap.cxx

#----------- External libraries
$(LIBGAS) :
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) libgas.a
	
$(GASPY) :
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) install

$(LIBNM) :
	cd $(NM)/build; make TARGET=$(TARGET) OPT=$(OPT) libnm.a

$(LIBUTIL) :
	cd $(UTIL)/build; make TARGET=$(TARGET) OPT=$(OPT) libutil.a
	
$(LIBGEOM) :
	cd $(GEOM)/build; make TARGET=$(TARGET) OPT=$(OPT) libgeom2.a
	
# Assume linux
PLAT := linux
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    PLAT := macosx
endif
$(LIBLUA):
	cd $(LUA); make $(PLAT) local

$(LUA_INCLUDE_DIR): $(LIBLUA)

$(LIBZLIB):
	cd $(ZLIB); make libz.a

# -----------  Static library --------------------------------
$(LIBRAD2) : $(RAD_OBJECTS)
ifneq ($(strip $(wildcard $(LIBRAD2))),)
	- rm $(LIBRAD2)
endif
	$(AR) rv $(LIBRAD2) $(RAD_OBJECTS)
	- ranlib $(LIBRAD2)

# ----------  Loadable library for Python under Linux --------
_librad2.so : librad2_wrap.o $(LIBRAD2) $(LIBNM) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LDFLAG) $(FLINK) $(OMPFLAG) -shared -o _librad2.so librad2_wrap.o \
		$(LIBRAD2) $(LIBNM) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)

# ----------  Loadable library for Python under Mac OS X --------
_librad2.dylib : librad2_wrap.o $(LIBRAD2) $(LIBNM) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LDFLAG) $(FLINK) $(OMPFLAG) -o _librad2.dylib librad2_wrap.o \
		$(LIBRAD2) $(LIBNM) $(LIBLUA) $(LIBZLIB) $(LIBGEOM) \
		-framework Python -bundle -bundle_loader $(PYTHON_BIN_DIR)/python 
		 
# ----------  Loadable library for Python under MS-Windows --------
_librad2.dll : librad2_wrap.o $(LIBRAD2) $(LIBNM) $(LIBPYTHON) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LDFLAG) $(FLINK) $(OMPFLAG) -shared -o _librad2.dll librad2_wrap.o \
		-L$(PYTHON_DIR)/libs \
		$(LIBRAD2) $(LIBNM) $(LIBPYTHON) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)

# The SWIG interface file
$(SRC)/librad2_wrap.cxx $(SRC)/librad2.py: $(SRC)/librad2.i $(SRC)/spectral_model.hh
	swig -python -c++ $(SRC)/librad2.i

librad2.py : $(SRC)/librad2.py
	cp $(SRC)/librad2.py .

ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    CXXFLAG3 := -c
else
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CXXFLAG3 := -c
    else
        # Unix/Linux and Darwin on OS X
        CXXFLAG3 := -c -fPIC
    endif
endif

# -----------  Swig wrapper -----------

librad2_wrap.o : $(SRC)/librad2_wrap.cxx $(LIBLUA)
	$(CXX) $(CXXFLAG3) $(SRC)/librad2_wrap.cxx \
	-I$(PYTHON_INCLUDE_DIR) -I$(LUA_INCLUDE_DIR)

# -----------  Object files -----------

spectral_model.o : $(SRC)/spectral_model.cxx $(SRC)/spectral_model.hh \
		$(GAS)/models/gas_data.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DWITH_SPRADIAN=$(WITH_SPRADIAN) $(SRC)/spectral_model.cxx \
		-I$(LUA_INCLUDE_DIR)

equilibrium_air.o : $(SRC)/equilibrium_air.cxx $(SRC)/equilibrium_air.hh \
		$(SRC)/spectral_model.hh $(SRC)/radiation_constants.hh $(GAS)/models/gas_data.hh
	$(CXX) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/equilibrium_air.cxx
	
lua_service.o : $(UTIL_SRC)/lua_service.hh $(UTIL_SRC)/lua_service.cxx $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) $(UTIL_SRC)/lua_service.cxx -I$(LUA_INCLUDE_DIR) -I$(ZLIB)
	
photaura.o : $(SRC)/photaura.cxx $(SRC)/photaura.hh \
		$(SRC)/photaura.hh $(SRC)/radiation_constants.hh \
		$(SRC)/spectral_model.hh $(GAS)/models/gas_data.hh
	$(CXX) $(CXXFLAG) $(SRC)/photaura.cxx -I$(LUA_INCLUDE_DIR)
	
radiator.o : $(SRC)/radiator.cxx $(SRC)/radiator.hh $(SRC)/photoionisation.hh
	$(CXX) $(CXXFLAG) $(SRC)/radiator.cxx -I$(LUA_INCLUDE_DIR)
	
diatomic_radiator.o : $(SRC)/diatomic_radiator.cxx $(SRC)/diatomic_radiator.hh \
		$(SRC)/radiator.hh $(SRC)/spectral_model.hh
	$(CXX) $(CXXFLAG) $(SRC)/diatomic_radiator.cxx -I$(LUA_INCLUDE_DIR)
	
diatomic_system.o : $(SRC)/diatomic_system.cxx $(SRC)/diatomic_system.hh \
		$(SRC)/diatomic_radiator.hh
	$(CXX) $(CXXFLAG) $(SRC)/diatomic_system.cxx -I$(LUA_INCLUDE_DIR)
	
diatomic_band.o : $(SRC)/diatomic_band.cxx $(SRC)/diatomic_band.hh
	$(CXX) $(CXXFLAG) $(SRC)/diatomic_band.cxx -I$(LUA_INCLUDE_DIR)

atomic_radiator.o : $(SRC)/atomic_radiator.cxx $(SRC)/atomic_radiator.hh \
		$(SRC)/radiator.hh $(SRC)/spectral_model.hh
	$(CXX) $(CXXFLAG) $(SRC)/atomic_radiator.cxx -I$(LUA_INCLUDE_DIR)
	
atomic_line.o : $(SRC)/atomic_line.cxx $(SRC)/atomic_line.hh
	$(CXX) $(CXXFLAG) $(SRC)/atomic_line.cxx -I$(LUA_INCLUDE_DIR)

electron_radiator.o : $(SRC)/electron_radiator.cxx $(SRC)/electron_radiator.hh \
		$(SRC)/radiator.hh
	$(CXX) $(CXXFLAG) $(SRC)/electron_radiator.cxx -I$(LUA_INCLUDE_DIR)

planck_radiator.o : $(SRC)/planck_radiator.cxx $(SRC)/planck_radiator.hh \
		$(SRC)/radiator.hh
	$(CXX) $(CXXFLAG) $(SRC)/planck_radiator.cxx -I$(LUA_INCLUDE_DIR)
	
spectra_pieces.o : $(SRC)/spectra_pieces.cxx $(SRC)/spectra_pieces.hh 
	$(CXX) $(CXXFLAG) $(SRC)/spectra_pieces.cxx  -I$(LUA_INCLUDE_DIR)
	
LOS_pieces.o : $(SRC)/LOS_pieces.cxx $(SRC)/LOS_pieces.hh 
	$(CXX) $(CXXFLAG) $(SRC)/LOS_pieces.cxx  -I$(LUA_INCLUDE_DIR)
	
photoionisation.o : $(SRC)/photoionisation.cxx $(SRC)/photoionisation.hh
	$(CXX) $(CXXFLAG) $(SRC)/photoionisation.cxx -I$(LUA_INCLUDE_DIR)
	
cr_reactions.o : $(SRC)/cr_reactions.cxx $(SRC)/cr_reactions.hh $(SRC)/cr_rr_coeffs.hh
	$(CXX) $(CXXFLAG) $(SRC)/cr_reactions.cxx -I$(LUA_INCLUDE_DIR)
	
cr_rr_coeffs.o : $(SRC)/cr_rr_coeffs.cxx $(SRC)/cr_rr_coeffs.hh
	$(CXX) $(CXXFLAG) $(SRC)/cr_rr_coeffs.cxx -I$(LUA_INCLUDE_DIR)
	
parade.o : $(SRC)/parade.cxx $(SRC)/parade.hh \
		$(SRC)/spectral_model.hh $(SRC)/radiation_constants.hh $(SRC)/parade_radiator.hh \
		$(GAS)/models/gas_data.hh
	$(CXX) $(CXXFLAG) $(OMPFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/parade.cxx
	
parade_radiator.o : $(SRC)/parade_radiator.cxx $(SRC)/parade_radiator.hh
	$(CXX) $(CXXFLAG) $(SRC)/parade_radiator.cxx -I$(LUA_INCLUDE_DIR)
	
spradian.o : $(SRC)/spradian.cxx $(SRC)/spradian.hh \
		$(SRC)/spectral_model.hh $(SRC)/radiation_constants.hh $(SRC)/spradian_radiator.hh \
		$(GAS)/models/gas_data.hh
	$(CXX) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/spradian.cxx
	
spradian_radiator.o : $(SRC)/spradian_radiator.cxx $(SRC)/spradian_radiator.hh
	$(CXX) $(CXXFLAG) $(SRC)/spradian_radiator.cxx -I$(LUA_INCLUDE_DIR)
	
radipac6.o : $(SPRADIAN_SRC)/radipac6.f90
	$(F90) $(F90FLAG) $(SPRADIAN_SRC)/radipac6.f90

mersenne.o : $(UTIL_SRC)/mersenne.cpp $(UTIL_SRC)/randomc.h
	$(CXX) $(CXXFLAG) $(UTIL_SRC)/mersenne.cpp
	
# ------------ Test python libraries and executables -------
	
photaura_test.py : $(SRC)/photaura_test.py
	- cp $(SRC)/photaura_test.py .
	
radGUI.py : $(SRC)/radGUI.py
	- cp $(SRC)/radGUI.py .

# ------------ Test executables ------------

cea2 : $(CEA2)/cea2.f
	gfortran -std=legacy -o cea2 $(CEA2)/cea2.f

photaura_test.x : photaura_test.o $(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LFLAG) $(FLINK) $(OMPFLAG) -o photaura_test.x photaura_test.o \
		$(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LUALINK) $(LIBGEOM)
		
photaura_test.o : $(SRC)/photaura_test.cxx 
	$(CXX) $(CXXFLAG) $(SRC)/photaura_test.cxx -I$(LUA_INCLUDE_DIR)
	
spradian_test.x : spradian_test.o $(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LFLAG) $(FLINK) $(OMPFLAG) -o spradian_test.x spradian_test.o \
		$(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LUALINK) $(LIBGEOM)
		
spradian_test.o : $(SRC)/spradian_test.cxx $(SRC)/spradian.cxx 
	$(CXX) $(CXXFLAG) $(SRC)/spradian_test.cxx -I$(LUA_INCLUDE_DIR)
	
run_spradian.x : run_spradian.o $(LIBRAD2) $(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LIBGEOM)
	$(CXXLINK) $(LFLAG) $(FLINK) -o run_spradian.x run_spradian.o \
		$(LIBRAD2) $(LIBNM) $(LIBUTIL) $(LIBGAS) $(LIBLUA) $(LIBZLIB) $(LUALINK) $(LIBGEOM) 
	
run_spradian.o : $(SRC)/run_spradian.cxx 
	$(CXX) $(CXXFLAG) $(SRC)/run_spradian.cxx -I$(LUA_INCLUDE_DIR)
	
