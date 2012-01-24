# \file    makefile
# \ingroup geom
# \brief   Makefile for geometry library used by the CFCFD codes.
# \author  PA Jacobs
# \version 1.0 adapted from gas_models and plot makefiles
#
# \usage   make
#          make TARGET=for_gnu
#          make test
#          make install
#          make clean
#----------------------------------------------------------------------

TARGET ?= for_gnu
include ../../util/source/systems.mk

UTIL_SRC        := ../../util/source
NM_SRC          := ../../nm/source
PY_NM_SRC       := ../../cfpylib/nm
SRC             := ../source
INSTALL_DIR     := $(HOME)/cfd_bin

#----------------------------------------------------------------------

EXE_FILES := 

# Assume Unix/Linux
LOADABLE_MODULE := _libgeom.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _libgeom.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _libgeom.dll
endif
ifeq ($(findstring darwin, $(SYSTEM)), darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _libgeom.dylib
endif

PYTHON_MODULE_FILES := libgeom.py geom.py gpath.py bezier.py \
	bezier_spline.py simple_boxes.py zero_solvers.py \
	block3d_grid.py roberts.py geom_elmer.py geom_mbcns.py \
	hyptan.py cluster_points.py valliammai.py

LIBGEOM := libgeom.a

GEOM_OBJECTS :=	geom.o bezier.o gpath.o


#------------------------------------------------------------

build_new: $(EXE_FILES) $(LIBGEOM) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES)
	@echo "Have now built the geometry library."

clean:
	-rm $(EXE_FILES) $(LIBGEOM) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES)
	-rm *.o
	-rm *.pyc

install: $(EXE_FILES) $(LIBGEOM) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES) \
		$(INSTALL_DIR)
	@echo "Install python pieces from the geometry functions."
	- cp -f $(PYTHON_MODULE_FILES) $(INSTALL_DIR)
	- cp -f $(LOADABLE_MODULE) $(INSTALL_DIR)

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

#--------------- EXE Files ----------------------------------
# None, presently.

#------------------------------------------------------------
# Static Library to link into other C programs.

$(LIBGEOM) : $(GEOM_OBJECTS)
	- rm $(LIBGEOM)
	ar rv $(LIBGEOM) $(GEOM_OBJECTS)
	- ranlib $(LIBGEOM)

#-----------------------------------------------------------------------
# Loadable library for Python.

# Dynamic link library for MINGW32 environment on MS-Windows
_libgeom.dll : libgeom_wrap.o $(LIBGEOM) $(LIBPYTHON)
	$(LINK) $(LDFLAG) -shared -o _libgeom.dll libgeom_wrap.o \
		-L$(PYTHON_DIR)/libs \
		$(LIBGEOM) $(LIBPYTHON) $(LLIB)

# Shared object for Unix/Linux
_libgeom.so : libgeom_wrap.o $(LIBGEOM)
	$(LINK) $(LDFLAG) -shared -o _libgeom.so libgeom_wrap.o $(LIBGEOM) $(LLIB)

# Dynamic library for Mac OS X
_libgeom.dylib : libgeom_wrap.o $(LIBGEOM)
	gcc $(LDFLAG) -dynamiclib -o _libgeom.dylib libgeom_wrap.o \
		$(LIBGEOM) $(LLIB) -framework Python

# The SWIG interface file and corresponding documentation
# Normally, we would not try to rebuild this file because 
# some computers will not have swig.

$(SRC)/libgeom_wrap.c : $(SRC)/libgeom.i $(SRC)/gpath.h $(SRC)/geom.h $(SRC)/bezier.h
	swig -python $(SRC)/libgeom.i

libgeom.py : $(SRC)/libgeom.py
	cp $(SRC)/libgeom.py .

$(SRC)/libgeom.py : $(SRC)/libgeom.i $(SRC)/gpath.h $(SRC)/geom.h $(SRC)/bezier.h
	swig -python $(SRC)/libgeom.i

geom.py : $(SRC)/geom.py
	cp $(SRC)/geom.py .

gpath.py : $(SRC)/gpath.py
	cp $(SRC)/gpath.py .

geom_elmer.py : $(SRC)/geom_elmer.py
	cp $(SRC)/geom_elmer.py .

geom_mbcns.py : $(SRC)/geom_mbcns.py
	cp $(SRC)/geom_mbcns.py .

simple_boxes.py : $(SRC)/simple_boxes.py
	cp $(SRC)/simple_boxes.py .

bezier.py : $(SRC)/bezier.py
	cp $(SRC)/bezier.py .

bezier_spline.py : $(SRC)/bezier_spline.py
	cp $(SRC)/bezier_spline.py .

zero_solvers.py: $(PY_NM_SRC)/zero_solvers.py
	cp $(PY_NM_SRC)/zero_solvers.py .

roberts.py: $(PY_NM_SRC)/roberts.py
	cp $(PY_NM_SRC)/roberts.py .

hyptan.py: $(PY_NM_SRC)/hyptan.py
	cp $(PY_NM_SRC)/hyptan.py .

valliammai.py: $(PY_NM_SRC)/valliammai.py
	cp $(PY_NM_SRC)/valliammai.py .

cluster_points.py : $(SRC)/cluster_points.py
	cp $(SRC)/cluster_points.py .

block3d_grid.py : $(SRC)/block3d_grid.py
	cp $(SRC)/block3d_grid.py .

#---------------------------------------------------------------
# Individual Object Files

geom.o : $(SRC)/geom.c $(SRC)/geom.h $(UTIL_SRC)/logfile.h
	$(COMPILE) $(CFLAG) $(SRC)/geom.c -o geom.o

bezier.o : $(SRC)/bezier.c $(SRC)/bezier.h $(SRC)/geom.h
	$(COMPILE) $(CFLAG) $(SRC)/bezier.c -o bezier.o

gpath.o : $(SRC)/gpath.c $(SRC)/gpath.h $(SRC)/bezier.h $(SRC)/geom.h
	$(COMPILE) $(CFLAG) $(SRC)/gpath.c -o gpath.o

ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    CFLAG3 := -c
else
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG3 := -c
    else
        # Unix/Linux and Darwin on OS X
        CFLAG3 := -c -fPIC
    endif
endif
libgeom_wrap.o : $(SRC)/libgeom_wrap.c
	$(COMPILE) $(CFLAG3) -I$(PYTHON_INCLUDE_DIR) \
		$(SRC)/libgeom_wrap.c

#--------------------------------------------------------------------
#
test:
	gcc -Wall -DTEST_GEOM ../source/geom.c -o geom_test.exe -lm
	gcc -Wall -o bezierd.exe ../source/bezierd.c ../source/bezier.c -lm
	gcc -Wall -o bezplot.exe ../source/bezplot.c ../source/bezier.c -lm
	gcc -Wall -o gpathd.exe ../source/gpathd.c ../source/gpath.c \
		 ../source/bezier.c ../source/geom.c -lm
	echo "Test executables built."

#------------------------ end of makefile ----------------------
