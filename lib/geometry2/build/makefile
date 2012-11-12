## \file makefile
## \ingroup geom2
## \brief Makefile for the experimental libgeom2 module.
## \author PJ
## \version 27-Dec-2005 start of the C++ experiment
## \version 11-Jan-2006 Genaralized to work on Ubuntu and 
##     different versions of Python. 
##     Presumably works with Cygwin and Mac-OSX.
##     g++ is the hard-coded compiler.
##
## \usage   make
##          make TARGET=for_gnu
##          make install
##          make clean
#----------------------------------------------------------------------

WITH_PYTHON_CALLBACKS ?= 1
TARGET ?= for_gnu
OPT ?= -O2
INSTALL_DIR ?= $(HOME)/e3bin
include ../../util/source/systems.mk

UTIL_SRC        := ../../util/source
NM              := ../../nm
NM_SRC          := $(NM)/source
SRC             := ../source
LUA             := ../../../extern/lua-5.1.4
LIBLUA          := $(LUA)/lib/liblua.a
LUA_INCLUDE_DIR := $(LUA)/include

#----------------------------------------------------------------------

EXE_FILES := geom_test.x gpath_test.x surface_test.x \
	secant_test.x geom_test.py gpath_test.py gpath_test_bezier.py \
	fobject_test.x surface_test.py volume_test.py meshvolume_test.py \
	triangle_patch_test.x triangle_patch_test.py triangle_join_test.py \
	revolved_surface_test.x revolved_surface_test.py \
	bezier_patch_test.py fit_spline.x helix_test.py polar_test.py \
	lgeom

# Assume Unix/Linux
LOADABLE_MODULE := _libgeom2.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _libgeom2.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _libgeom2.dll
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _libgeom2.dylib
endif

PYTHON_MODULE_FILES := libgeom2.py blockgrid2d.py blockgrid3d.py

LIBGEOM2 := libgeom2.a

GEOM_OBJECTS :=	geom.o gpath.o secant.o surface.o fobject.o volume.o nelmin.o \
	gpath_utils.o no_fuss_linear_algebra.o nurbs.o golden_section_search.o \
	nurbs_utils.o zero_system.o zero_finders.o
ifeq ($(WITH_PYTHON_CALLBACKS), 1)
    GEOM_OBJECTS += pypath.o pysurface.o pyvolume.o
endif

LUA_OBJECTS := l_geom.o l_gpath.o l_gpath_utils.o l_nurbs_utils.o l_surface.o
LIBLUAGEOM2 := libluageom2.a
LUAGEOM2_MODULE := geometry.so
LUA_TEST_FILES := geom_test.lua gpath_test.lua

#------------------------------------------------------------

build_new: $(EXE_FILES) $(LIBGEOM2) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES)
	@echo "Geometry-library-2 should be up to date."

clean:
	-rm $(EXE_FILES) $(LIBGEOM2) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES)
	-rm *.o
	-rm *.pyc
	-rm *.wrl
	-rm *.vtk
	-rm *.a
	-rm $(SRC)/libgeom2.py
	-rm $(SRC)/libgeom2_wrap.cxx
	-rm $(LUAGEOM2_MODULE) $(LUA_TEST_FILES) 
	- cd $(LUA); make clean

install: $(EXE_FILES) $(LIBGEOM2) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES) \
		$(INSTALL_DIR)
	@echo "Install python pieces from the geometry-2 functions."
	- cp -f $(PYTHON_MODULE_FILES) $(INSTALL_DIR)
	- cp -f $(LOADABLE_MODULE) $(INSTALL_DIR)
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
	- cp $(INSTALL_DIR)/_libgeom2.dylib $(INSTALL_DIR)/_libgeom2.so
endif

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

lua_shared_library: $(LUA_SHARED_LIBRARY)

#--------------- EXE Files ----------------------------------

geom_test.x : geom_test.o $(LIBGEOM2)
	$(CXXLINK) -o geom_test.x geom_test.o $(LIBGEOM2)

geom_test.py : $(SRC)/geom_test.py
	cp $(SRC)/geom_test.py .

gpath_test.py : $(SRC)/gpath_test.py
	cp $(SRC)/gpath_test.py .

gpath_test_bezier.py : $(SRC)/gpath_test_bezier.py
	cp $(SRC)/gpath_test_bezier.py .

helix_test.py : $(SRC)/helix_test.py
	cp $(SRC)/helix_test.py .

polar_test.py : $(SRC)/polar_test.py
	cp $(SRC)/polar_test.py .

gpath_test.x : gpath_test.o $(LIBGEOM2)
	$(CXXLINK) -o gpath_test.x gpath_test.o $(LIBGEOM2)

test_create_cubic.x : test_create_cubic.o $(LIBGEOM2) no_fuss_linear_algebra.o
	$(CXXLINK) -o test_create_cubic.x test_create_cubic.o $(LIBGEOM2) no_fuss_linear_algebra.o


secant_test.x : secant_test.o secant.o
	$(CXXLINK) -o secant_test.x secant_test.o secant.o

surface_test.x : surface_test.o $(LIBGEOM2)
	$(CXXLINK) -o surface_test.x surface_test.o $(LIBGEOM2)

triangle_patch_test.x : triangle_patch_test.o $(LIBGEOM2)
	$(CXXLINK) -o triangle_patch_test.x triangle_patch_test.o $(LIBGEOM2)

revolved_surface_test.x : revolved_surface_test.o $(LIBGEOM2)
	$(CXXLINK) -o revolved_surface_test.x revolved_surface_test.o $(LIBGEOM2)

fobject_test.x : fobject_test.o $(LIBGEOM2)
	$(CXXLINK) -o fobject_test.x fobject_test.o $(LIBGEOM2)

surface_test.py : $(SRC)/surface_test.py
	cp $(SRC)/surface_test.py .

triangle_patch_test.py : $(SRC)/triangle_patch_test.py
	cp $(SRC)/triangle_patch_test.py .

bezier_patch_test.py : $(SRC)/bezier_patch_test.py
	cp $(SRC)/bezier_patch_test.py .

triangle_join_test.py : $(SRC)/triangle_join_test.py
	cp $(SRC)/triangle_join_test.py .

revolved_surface_test.py : $(SRC)/revolved_surface_test.py
	cp $(SRC)/revolved_surface_test.py .

volume_test.py : $(SRC)/volume_test.py
	cp $(SRC)/volume_test.py .

meshvolume_test.py : $(SRC)/meshvolume_test.py
	cp $(SRC)/meshvolume_test.py .

blockgrid2d.py : $(SRC)/blockgrid2d.py
	cp $(SRC)/blockgrid2d.py .

blockgrid3d.py : $(SRC)/blockgrid3d.py
	cp $(SRC)/blockgrid3d.py .

fit_spline.x : fit_spline.o best_fit_spline.o $(LIBGEOM2)
	$(CXXLINK) -o fit_spline.x fit_spline.o best_fit_spline.o $(LIBGEOM2)

# A Lua interpreter with geometry functions built in.
lgeom : $(SRC)/lgeom.cxx $(LIBGEOM2) $(LIBLUAGEOM2) $(LIBLUA)
	$(CXXLINK) $(LDFLAG) -o lgeom $(SRC)/lgeom.cxx $(LIBLUAGEOM2) $(LIBGEOM2) $(LIBLUA) $(LLIB) \
	-I$(LUA_INCLUDE_DIR) -ldl -lreadline -lhistory -lncurses

# A library module for a Lua interpreter -- lua> require "geometry"
$(LUAGEOM2_MODULE) : $(LUA_OBJECTS) $(LIBGEOM2) $(LIBLUA)
	$(CXXLINK) $(LDFLAG) -shared -o geometry.so $(LUA_OBJECTS) $(LIBGEOM2) $(LIBLUA) $(LLIB)

geom_test.lua : $(SRC)/geom_test.lua
	cp $(SRC)/geom_test.lua .

gpath_test.lua : $(SRC)/gpath_test.lua
	cp $(SRC)/gpath_test.lua .

#------------------------------------------------------------
# Static Library to link into other C++ programs.

$(LIBGEOM2) : $(GEOM_OBJECTS)
ifneq ($(strip $(wildcard $(LIBGEOM2))),)
	- rm $(LIBGEOM2)
endif
	$(AR) rv $(LIBGEOM2) $(GEOM_OBJECTS)
	- ranlib $(LIBGEOM2)

# Static library to link with lua-enabled C++ programs which require
# access to geometry objects/functions
$(LIBLUAGEOM2) : $(LIBGEOM2) $(LUA_OBJECTS)
	$(AR) rv $(LIBLUAGEOM2) $(LUA_OBJECTS)
	- ranlib $(LIBLUAGEOM2)


# Assume linux
PLAT := linux
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    PLAT := macosx
endif
$(LIBLUA):
	cd $(LUA); make $(PLAT) local

#-----------------------------------------------------------------------
# Loadable library for Python.

# Dynamic link library for MINGW32 environment on MS-Windows
_libgeom2.dll : libgeom2_wrap.o $(LIBGEOM2) $(LIBPYTHON)
	$(CXXLINK) $(LDFLAG) -shared -o _libgeom2.dll libgeom2_wrap.o \
		-L$(PYTHON_DIR)/libs \
		$(LIBGEOM2) $(LIBPYTHON) $(LLIB)

# Shared object for Unix/Linux
_libgeom2.so : libgeom2_wrap.o $(LIBGEOM2)
	$(CXXLINK) $(LDFLAG) -shared -o _libgeom2.so libgeom2_wrap.o \
		$(LIBGEOM2) $(LLIB)

# Dynamic library for Mac OS X
_libgeom2.dylib : libgeom2_wrap.o $(LIBGEOM2)
	$(CXXLINK) $(LDFLAG) -o _libgeom2.dylib libgeom2_wrap.o \
		$(LIBGEOM2) $(LLIB) -framework Python \
		-bundle -bundle_loader $(PYTHON_BIN_DIR)/python

# The SWIG interface file and corresponding documentation
# Normally, we would not try to rebuild this file because 
# some computers will not have swig.

$(SRC)/libgeom2_wrap.cxx $(SRC)/libgeom2.py : $(SRC)/libgeom2.i \
		$(SRC)/gpath.hh $(SRC)/geom.hh $(SRC)/pypath.hh \
		$(SRC)/surface.hh $(NM_SRC)/secant.hh $(SRC)/volume.hh \
		$(NM_SRC)/fobject.hh $(SRC)/pysurface.hh $(SRC)/pyvolume.hh
	swig -python -c++ -w512 \
	     -DWITH_PYTHON_CALLBACKS=$(WITH_PYTHON_CALLBACKS) \
	     $(SRC)/libgeom2.i

libgeom2.py : $(SRC)/libgeom2.py
	cp $(SRC)/libgeom2.py .

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
libgeom2_wrap.o : $(SRC)/libgeom2_wrap.cxx
	$(CXX) $(CXXFLAG3) -I$(PYTHON_INCLUDE_DIR) \
		$(SRC)/libgeom2_wrap.cxx

#-----------------------------------------------------------------

geom.o : $(SRC)/geom.cxx $(SRC)/geom.hh
	$(CXX) $(CXXFLAG) $(SRC)/geom.cxx

geom_test.o : $(SRC)/geom_test.cxx $(SRC)/geom.hh
	$(CXX) $(CXXFLAG) $(SRC)/geom_test.cxx

test_create_cubic.o : $(SRC)/test_create_cubic.cxx $(SRC)/gpath_utils.hh $(SRC)/geom.hh
	$(CXX) $(CXXFLAG) $(SRC)/test_create_cubic.cxx

golden_section_search.o : $(NM_SRC)/golden_section_search.cxx $(NM_SRC)/golden_section_search.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/golden_section_search.cxx

gpath_test.o : $(SRC)/gpath_test.cxx $(SRC)/gpath.hh $(SRC)/geom.hh
	$(CXX) $(CXXFLAG) $(SRC)/gpath_test.cxx

gpath.o : $(SRC)/gpath.cxx $(SRC)/gpath.hh $(SRC)/geom.hh $(NM_SRC)/secant.hh \
		$(NM_SRC)/fobject.hh
	$(CXX) $(CXXFLAG) $(SRC)/gpath.cxx

pypath.o : $(SRC)/pypath.cxx $(SRC)/pypath.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh $(NM_SRC)/secant.hh
ifeq ($(TARGET), for_pgi)
	$(CXX) $(CXXFLAG) -I$(PYTHON_INCLUDE_DIR) $(SRC)/pypath.cxx
else
	$(CXX) $(CXXFLAG) -Wno-long-long -I$(PYTHON_INCLUDE_DIR) $(SRC)/pypath.cxx
endif


gpath_utils.o : $(SRC)/gpath_utils.cxx $(SRC)/gpath_utils.hh $(SRC)/gpath.hh $(SRC)/geom.hh \
		$(NM_SRC)/no_fuss_linear_algebra.hh
	$(CXX) $(CXXFLAG) $(SRC)/gpath_utils.cxx

no_fuss_linear_algebra.o : $(NM_SRC)/no_fuss_linear_algebra.cxx $(NM_SRC)/no_fuss_linear_algebra.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/no_fuss_linear_algebra.cxx

nurbs.o : $(SRC)/nurbs.cxx $(SRC)/nurbs.hh
	$(CXX) $(CXXFLAG) $(SRC)/nurbs.cxx

nurbs_utils.o : $(SRC)/nurbs_utils.cxx $(SRC)/nurbs_utils.hh $(SRC)/gpath.hh $(SRC)/geom.hh \
		$(SRC)/nurbs.hh
	$(CXX) $(CXXFLAG) $(SRC)/nurbs_utils.cxx

secant_test.o : $(NM_SRC)/secant_test.cxx $(NM_SRC)/secant.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/secant_test.cxx

secant.o : $(NM_SRC)/secant.cxx $(NM_SRC)/secant.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/secant.cxx

nelmin.o : $(NM_SRC)/nelmin.cxx $(NM_SRC)/nelmin.hh $(NM_SRC)/fobject.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/nelmin.cxx

surface_test.o : $(SRC)/surface_test.cxx $(SRC)/surface.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh
	$(CXX) $(CXXFLAG) $(SRC)/surface_test.cxx

triangle_patch_test.o : $(SRC)/triangle_patch_test.cxx $(SRC)/surface.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh
	$(CXX) $(CXXFLAG) $(SRC)/triangle_patch_test.cxx

revolved_surface_test.o : $(SRC)/revolved_surface_test.cxx $(SRC)/surface.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh
	$(CXX) $(CXXFLAG) $(SRC)/revolved_surface_test.cxx

surface.o : $(SRC)/surface.cxx $(SRC)/surface.hh $(SRC)/geom.hh \
		$(SRC)/gpath.hh $(NM_SRC)/secant.hh
	$(CXX) $(CXXFLAG) $(SRC)/surface.cxx

pysurface.o : $(SRC)/pysurface.cxx $(SRC)/pysurface.hh $(SRC)/surface.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh $(NM_SRC)/secant.hh
ifeq ($(TARGET), for_pgi)
	$(CXX) $(CXXFLAG) -I$(PYTHON_INCLUDE_DIR) $(SRC)/pysurface.cxx
else
	$(CXX) $(CXXFLAG) -Wno-long-long -I$(PYTHON_INCLUDE_DIR) $(SRC)/pysurface.cxx
endif
fobject.o : $(NM_SRC)/fobject.cxx $(NM_SRC)/fobject.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/fobject.cxx

fobject_test.o : $(NM_SRC)/fobject_test.cxx $(NM_SRC)/fobject.hh
	$(CXX) $(CXXFLAG) $(NM_SRC)/fobject_test.cxx

volume.o : $(SRC)/volume.cxx $(SRC)/volume.hh $(SRC)/surface.hh $(SRC)/geom.hh \
		$(SRC)/gpath.hh $(NM_SRC)/secant.hh
	$(CXX) $(CXXFLAG) $(SRC)/volume.cxx

pyvolume.o : $(SRC)/pyvolume.cxx $(SRC)/pyvolume.hh $(SRC)/volume.hh $(SRC)/surface.hh \
		$(SRC)/geom.hh $(SRC)/gpath.hh $(NM_SRC)/secant.hh
ifeq ($(TARGET), for_pgi)
	$(CXX) $(CXXFLAG) -I$(PYTHON_INCLUDE_DIR) $(SRC)/pyvolume.cxx
else
	$(CXX) $(CXXFLAG) -Wno-long-long -I$(PYTHON_INCLUDE_DIR) $(SRC)/pyvolume.cxx
endif

best_fit_spline.o : $(SRC)/best_fit_spline.cxx $(SRC)/best_fit_spline.hh $(SRC)/gpath.hh \
		$(NM_SRC)/fobject.hh $(NM_SRC)/nelmin.hh
	$(CXX) $(CXXFLAG) $(SRC)/best_fit_spline.cxx

fit_spline.o : $(SRC)/fit_spline.cxx
	$(CXX) $(CXXFLAG) $(SRC)/fit_spline.cxx

l_geom.o : $(SRC)/l_geom.cxx $(SRC)/l_gpath.hh $(SRC)/geom.hh $(LIBLUA)
	$(CXX) $(CXXFLAG) $(SRC)/l_geom.cxx -I$(LUA_INCLUDE_DIR)

l_gpath.o : $(SRC)/l_gpath.cxx $(SRC)/l_gpath.hh $(SRC)/gpath.hh $(SRC)/geom.hh $(SRC)/l_geom.hh $(LIBLUA)
	$(CXX) $(CXXFLAG) $(SRC)/l_gpath.cxx -I$(LUA_INCLUDE_DIR)

l_gpath_utils.o : $(SRC)/l_gpath_utils.cxx $(SRC)/l_gpath_utils.hh $(SRC)/gpath_utils.hh \
	$(SRC)/gpath.hh $(SRC)/geom.hh $(SRC)/l_geom.hh $(LIBLUA)
	$(CXX) $(CXXFLAG) $(SRC)/l_gpath_utils.cxx -I$(LUA_INCLUDE_DIR)

l_nurbs_utils.o : $(SRC)/l_nurbs_utils.cxx $(SRC)/l_nurbs_utils.hh $(SRC)/l_gpath_utils.hh $(SRC)/gpath_utils.hh \
	$(SRC)/gpath.hh $(SRC)/geom.hh $(SRC)/l_geom.hh $(UTIL_SRC)/lunar.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) $(SRC)/l_nurbs_utils.cxx -I$(LUA_INCLUDE_DIR)

l_surface.o : $(SRC)/l_surface.cxx $(SRC)/l_surface.hh $(SRC)/surface.cxx $(SRC)/surface.hh \
	$(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) $(SRC)/l_surface.cxx -I$(LUA_INCLUDE_DIR)

zero_system.o : $(NM_SRC)/zero_system.cxx $(NM_SRC)/zero_system.hh $(NM_SRC)/no_fuss_linear_algebra.hh
	$(CXX) $(CXXFLAG) $(ARRAY_SIZES) $(NM_SRC)/zero_system.cxx
	
zero_finders.o : $(NM_SRC)/zero_finders.cxx $(NM_SRC)/zero_finders.hh \
		$(NM_SRC)/zero_system.hh $(NM_SRC)/no_fuss_linear_algebra.hh
	$(CXX) $(CXXFLAG) $(ARRAY_SIZES) $(NM_SRC)/zero_finders.cxx
