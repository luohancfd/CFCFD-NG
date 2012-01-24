# basic makefile for pure-fluid thermo functions

# Some default compile options
TARGET ?= for_gnu
OPT ?= -O2
include ../../util/source/systems.mk

# Assume Unix/Linux
PYTHONFLUID_MODULE := libfluid.py _libfluid.so
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
   # Darwin environment on a Mac OS X
  PYTHONFLUID_MODULE := libfluid.py _libfluid.dylib
endif
CXXFLAG3 := -c -fPIC

all : fluid_thermo_test.exe co2_cycle.exe lab_co2_cycle.exe $(PYTHONFLUID_MODULE)
	echo "Should be up to date."

clean:
	- rm *.o *.exe $(PYTHONFLUID_MODULE) libfluid_wrap.cxx *.pyc


libfluid_wrap.cxx libfluid.py: libfluid.i fluid_thermo.hh
	swig -python -c++ -w512 libfluid.i

libfluid_wrap.o : libfluid_wrap.cxx
	$(CXX) $(CXXFLAG3)  libfluid_wrap.cxx -I$(PYTHON_INCLUDE_DIR)

# Shared object for Unix/Linux
_libfluid.so : libfluid.py libfluid_wrap.o fluid_thermo.o
	$(CXXLINK) $(LDFLAG) -shared -o _libfluid.so \
	libfluid_wrap.o fluid_thermo.o $(LLIB)

# Dynamic library for Mac OS X
_libfluid.dylib : libfluid.py libfluid_wrap.o fluid_thermo.o
	$(CXXLINK) $(LDFLAG) -o _libfluid.dylib \
	libfluid_wrap.o fluid_thermo.o \
	$(LLIB) -framework Python \
	-bundle -bundle_loader $(PYTHON_BIN_DIR)/python

fluid_thermo_test.exe : fluid_thermo_test.o fluid_thermo.o
	g++ -o fluid_thermo_test.exe fluid_thermo_test.o fluid_thermo.o -lm

fluid_thermo_test.o : fluid_thermo_test.cxx fluid_thermo.hh
	g++ -c -Wall -pedantic fluid_thermo_test.cxx

co2_cycle.exe : co2_cycle.o fluid_thermo.o
	g++ -o co2_cycle.exe co2_cycle.o fluid_thermo.o -lm

co2_cycle.o : co2_cycle.cxx fluid_thermo.hh
	g++ -c -Wall -pedantic co2_cycle.cxx

lab_co2_cycle.exe : lab_co2_cycle.o fluid_thermo.o
	g++ -o lab_co2_cycle.exe lab_co2_cycle.o fluid_thermo.o -lm

lab_co2_cycle.o : lab_co2_cycle.cxx fluid_thermo.hh
	g++ -c -Wall -pedantic lab_co2_cycle.cxx

fluid_thermo.o : fluid_thermo.cxx fluid_thermo.hh
	g++ -c -Wall -pedantic fluid_thermo.cxx
