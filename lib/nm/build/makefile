## \file makefile
## \ingroup nm
## \brief Makefile for numerical methods
## \author RJG
## \date 26-Apr-2006
## \version 26-Apr-2006 : factored certain parts into a library
##                        testing programs are a separate target

TARGET ?= for_gnu
OPT ?= -O2
INSTALL_DIR ?= $(HOME)/cfd_bin
INSTALL_DIR2 ?= $(HOME)/e3bin

include ../../util/source/systems.mk

SRC		:= ../source
LIBNM		:= libnm.a

UTIL		:= ../../util
LIBUTIL		:= $(UTIL)/build/libutil.a

PYTHON_MODULE_FILES := libnm.py

NM_OBJECTS := no_fuss_linear_algebra.o \
	zero_system.o \
	zero_finders.o \
	ode_system.o \
	ode_solver.o \
	ode_step.o \
	curve_fit.o \
	segmented_curve_fit.o \
	exponential_integrals.o \
	quadrature.o \
	linear_interpolation.o \
	lu_decomp.o \
	ridder.o

TEST_EXE_FILES := fobject_test.x \
	nelmin_test.x \
	ode_test.x \
	no_fuss_test.x \
	lu_decomp_test.x \
	ridder_test.x
#	zero_finders_test.x

LOADABLE_MODULE := _libnm.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _libnm.dll 
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _libnm.dll 
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _libnm.dylib
endif

build_new: $(LIBNM) $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES)

python_module : $(PYTHON_MODULE_FILES) $(LOADABLE_MODULE) install

test : $(TEST_EXE_FILES) 

clean :
	- rm *.x *.o *.so *.a *.py
	- rm $(SRC)/libnm_wrap.cxx $(SRC)/libnm.py
	- rm $(LOADABLE_MODULE)
	- cd ../../util/build; make clean

install: $(LOADABLE_MODULE) $(PYTHON_MODULE_FILES) $(INSTALL_DIR)  $(INSTALL_DIR2)
	@echo "Install python pieces from the nm functions."
	- cp -f $(PYTHON_MODULE_FILES) $(INSTALL_DIR)
	- cp -f $(LOADABLE_MODULE) $(INSTALL_DIR)
	- cp -f $(PYTHON_MODULE_FILES) $(INSTALL_DIR2)
	- cp -f $(LOADABLE_MODULE) $(INSTALL_DIR2)

ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
	- cp $(INSTALL_DIR)/_libnm.dylib $(INSTALL_DIR)/_libnm.so
endif

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

$(INSTALL_DIR2):	
	mkdir $(INSTALL_DIR2)

#-------- Python module (using SWIG) ---------

$(LIBNM) : $(NM_OBJECTS)
ifneq ($(strip $(wildcard $(LIBNM))),)
	- rm $(LIBNM)
endif
	$(AR) rv $(LIBNM) $(NM_OBJECTS)
	- ranlib $(LIBNM)

libnm.py : $(SRC)/libnm.py
	cp $(SRC)/libnm.py .

$(LIBUTIL) :
	cd $(UTIL)/build; make TARGET=$(TARGET) OPT=$(OPT) libutil.a

$(SRC)/libnm.py : $(SRC)/libnm.i $(SRC)/no_fuss_linear_algebra.hh \
	$(SRC)/zero_system.hh $(SRC)/zero_finders.hh
	swig -python -c++ -w512 -Wall $(SRC)/libnm.i

$(SRC)/libnm_wrap.cxx : $(SRC)/libnm.i $(SRC)/no_fuss_linear_algebra.hh \
	$(SRC)/zero_system.hh $(SRC)/zero_finders.hh
	swig -python -c++ -w512 -Wall $(SRC)/libnm.i

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
libnm_wrap.o : $(SRC)/libnm_wrap.cxx
	$(CXXCOMPILE) $(CXXFLAG3) -I$(PYTHON_INCLUDE_DIR) \
		$(SRC)/libnm_wrap.cxx

# Dynamic link library for MINGW32 environment on MS-Windows
_libnm.dll : libnm_wrap.o $(LIBUTIL) $(LIBNM) $(LIBPYTHON)
	$(CXXLINK) $(LDFLAG) -shared -o _libnm.dll libnm_wrap.o \
		-L$(PYTHON_DIR)/libs \
		$(LIBUTIL) $(LIBNM) $(LIBPYTHON) $(LLIB)

# Shared object for Unix/Linux
_libnm.so : libnm_wrap.o $(LIBUTIL) $(LIBNM)
	$(CXXLINK) $(LDFLAG) -shared -o _libnm.so libnm_wrap.o \
		$(LIBUTIL) $(LIBNM) $(LLIB)

# Dynamic library for Mac OS X
_libnm.dylib : libnm_wrap.o $(LIBUTIL) $(LIBNM)
	$(CXXLINK) $(LDFLAG) -o _libnm.dylib libnm_wrap.o \
		$(LIBUTIL) $(LIBNM) $(LLIB) -framework Python \
		-bundle -bundle_loader $(PYTHON_BIN_DIR)/python

#-------- executable -----------

nelmin_test.x : nelmin_test.o nelmin.o fobject.o $(LIBNM)
	$(CXXLINK) -o nelmin_test.x nelmin_test.o nelmin.o fobject.o $(LIBNM)

fobject_test.x : fobject_test.o fobject.o $(LIBNM)
	$(CXXLINK) -o fobject_test.x fobject_test.o fobject.o $(LIBNM)

ode_test.x : ode_test.o ode_step.o ode_system.o ode_solver.o no_fuss_linear_algebra.o
	$(CXXLINK) -o ode_test.x ode_test.o ode_step.o ode_system.o ode_solver.o no_fuss_linear_algebra.o

no_fuss_test.x : no_fuss_test.o no_fuss_linear_algebra.o
	$(CXXLINK) -o no_fuss_test.x no_fuss_test.o no_fuss_linear_algebra.o

lu_decomp_test.x : lu_decomp_test.o lu_decomp.o no_fuss_linear_algebra.o
	$(CXXLINK) -o lu_decomp_test.x lu_decomp_test.o lu_decomp.o no_fuss_linear_algebra.o

zero_finders_test.x : zero_finders_test.o zero_finders.o
	$(CXXLINK) -o zero_finders_test.x zero_finders_test.o zero_finders.o

Richardson_extrapolation_test.x : $(SRC)/Richardson_extrapolation_test.cxx \
	Richardson_extrapolation.o
	$(CXXLINK) -o Richardson_extrapolation_test.x $(SRC)/Richardson_extrapolation_test.cxx \
	Richardson_extrapolation.o

ridder_test.x : ridder_test.o $(LIBNM)
	$(CXXLINK) -o ridder_test.x ridder_test.o $(LIBNM)

#-------- objects for testing executables ---------
nelmin_test.o : $(SRC)/nelmin_test.cxx $(SRC)/nelmin.hh $(SRC)/fobject.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/nelmin_test.cxx

fobject_test.o : $(SRC)/fobject_test.cxx $(SRC)/fobject.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/fobject_test.cxx

ode_test.o : $(SRC)/ode_test.cxx $(SRC)/ode_solver.hh $(SRC)/ode_system.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/ode_test.cxx

no_fuss_test.o : $(SRC)/no_fuss_test.cxx $(SRC)/no_fuss_linear_algebra.hh 
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/no_fuss_test.cxx

lu_decomp_test.o : $(SRC)/lu_decomp_test.cxx $(SRC)/lu_decomp.cxx $(SRC)/lu_decomp.hh \
	$(SRC)/no_fuss_linear_algebra.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/lu_decomp_test.cxx

zero_finders_test.o : $(SRC)/zero_finders_test.cxx $(SRC)/zero_system.hh $(SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/zero_finders_test.cxx

ridder_test.o : $(SRC)/ridder_test.cxx $(SRC)/ridder.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/ridder_test.cxx


#-------- objects compilation -----------
nelmin.o : $(SRC)/nelmin.cxx $(SRC)/nelmin.hh $(SRC)/fobject.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/nelmin.cxx

fobject.o : $(SRC)/fobject.cxx $(SRC)/fobject.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/fobject.cxx

curve_fit.o : $(SRC)/curve_fit.cxx $(SRC)/curve_fit.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/curve_fit.cxx

segmented_curve_fit.o : $(SRC)/segmented_curve_fit.cxx $(SRC)/segmented_curve_fit.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/segmented_curve_fit.cxx

ode_system.o : $(SRC)/ode_system.cxx $(SRC)/ode_system.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/ode_system.cxx

ode_step.o : $(SRC)/ode_step.cxx $(SRC)/ode_step.hh $(SRC)/ode_system.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/ode_step.cxx

ode_solver.o : $(SRC)/ode_solver.cxx $(SRC)/ode_solver.hh $(SRC)/ode_system.hh \
		$(SRC)/ode_step.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/ode_solver.cxx

no_fuss_linear_algebra.o : $(SRC)/no_fuss_linear_algebra.cxx $(SRC)/no_fuss_linear_algebra.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/no_fuss_linear_algebra.cxx

zero_finders.o : $(SRC)/zero_finders.cxx $(SRC)/zero_finders.hh \
		$(SRC)/zero_system.hh $(SRC)/no_fuss_linear_algebra.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/zero_finders.cxx

zero_system.o : $(SRC)/zero_system.cxx $(SRC)/zero_system.hh $(SRC)/no_fuss_linear_algebra.hh
	$(CXXCOMPILE) $(CXXFLAG) $(ARRAY_SIZES) $(SRC)/zero_system.cxx

Richardson_extrapolation.o : $(SRC)/Richardson_extrapolation.cxx $(SRC)/Richardson_extrapolation.hh \
	$(SRC)/functor.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/Richardson_extrapolation.cxx

exponential_integrals.o : $(SRC)/exponential_integrals.cxx $(SRC)/exponential_integrals.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/exponential_integrals.cxx

quadrature.o : $(SRC)/quadrature.cxx $(SRC)/quadrature.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/quadrature.cxx

linear_interpolation.o : $(SRC)/linear_interpolation.cxx $(SRC)/linear_interpolation.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/linear_interpolation.cxx

lu_decomp.o : $(SRC)/lu_decomp.cxx $(SRC)/lu_decomp.hh $(SRC)/no_fuss_linear_algebra.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/lu_decomp.cxx

ridder.o : $(SRC)/ridder.cxx $(SRC)/ridder.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/ridder.cxx
