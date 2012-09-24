# \file    makefile
# \ingroup L1d3
# \brief   Makefile for L1d3 and friends.
# \author  PA Jacobs
# \version 1.0 adapted from L1d2
#
# \usage   make
#          make TARGET=for_gnu
#          make install
#          make clean
#          make tags
#          make tarfile
#----------------------------------------------------------------------

TARGET ?= for_gnu
WITH_PYTHON_BITS ?= 1
INSTALL_DIR ?= $(HOME)/e3bin

# We seem to be having trouble with optimizing compilers.
# This optimization level may be over-ridden at the command line.
OPT ?= -O2

NSPECD ?= 30
NDIM   ?= 2000
include ../../../lib/util/source/systems.mk
LLIB += -lpopt -ldl

SRC             := ../source
EILMER3         := ../../../app/eilmer3
PLOT            := ../../../lib/plot
UTIL            := ../../../lib/util
UTIL_SRC        := $(UTIL)/source
LIBUTIL         := $(UTIL)/build/libutil.a
GAS             := ../../../lib/gas
LIBGAS          := $(GAS)/build/libgas.a
NM              := ../../../lib/nm
NM_SRC          := $(NM)/source
LIBNM           := $(NM)/build/libnm.a
LIBINIPARSER    := $(UTIL)/iniparser/libiniparser.a
LUA             := ../../../extern/lua-5.1.4
LIBLUA          := $(LUA)/lib/liblua.a
LUA_INCLUDE_DIR := $(LUA)/include
ZLIB            := ../../../extern/zlib-1.2.3
LIBZLIB         := $(ZLIB)/libz.a

EXEFILES = l1d.exe \
	l_prep.exe \
	l_post.exe \
	l_hist.exe \
	sptime.exe \
	piston.exe \
	l_trim.exe \
	$(PLOT)/build/mb_cont.exe

OBJECTS = l_misc.o \
	l_kernel.o \
	l_adapt.o \
	l_tstep.o \
	l_bc.o \
	l_io.o \
	l_rivp.o \
	roberts.o

PY_FILES = l_script.py gaspy.py

# Assume Unix/Linux
LOADABLE_MODULE := _gaspy.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _gaspy.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _gaspy.dll
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _gaspy.dylib
endif
PY_FILES += $(LOADABLE_MODULE)

#-----------------------------------------------------------

build_all : $(EXEFILES) $(PY_FILES)
	@echo "------------------------------"
	@echo "L1d3 should be up-to-date now."
	@echo "------------------------------"

install: $(INSTALL_DIR) $(EXEFILES) $(PY_FILES)
	@echo Put executable files into $(INSTALL_DIR)
	cp -f $(EXEFILES) $(INSTALL_DIR)
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd $(NM)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) INSTALL_DIR2=$(INSTALL_DIR) install
	cd $(PLOT)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cp -f $(PY_FILES) $(INSTALL_DIR)

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

clean:
	@echo Removing library file
	- rm l1d3.a
	@echo Removing object files
	- rm $(OBJECTS) l1d.o l_prep.o l_post.o \
		l_hist.o sptime.o piston.o l_trim.o
	@echo Removing executable files
	- rm $(EXEFILES)
	- rm $(PY_FILES) *.pyc
	@echo Cleaning up Gas Models
	- cd $(GAS)/build; make clean
	- cd $(NM)/build; make clean
	- cd $(PLOT)/build; make clean
	- cd $(UTIL)/build; make clean
	- cd $(UTIL)/iniparser; make clean
	- cd $(LUA); make clean

doc: html_doc latex_doc
	@echo "Done both HTML and Latex doc."

html_doc:
	epydoc --html -o ../doc/l_script \
		--name "Using l_script to define L1d simulations." \
		$(SRC)/l_script.py

latex_doc:
	epydoc --latex -o ../doc/l_script \
		--name "Using l_script to define L1d simulations." \
		$(SRC)/l_script.py
	cd ../doc/l_script; \
		latex l_script; \
		makeindex l_script; \
		latex l_script; \
		dvips -o l_script.ps l_script; \
		ps2pdf l_script.ps l_script.pdf

tags:
	cd ../source; \
	etags ./*.cxx ./*.hh \
		$(GAS)/kinetics/*.cxx $(GAS)/kinetics/*.hh \
		$(GAS)/models/*.cxx $(GAS)/models/*.hh \
		$(UTIL_SRC)/*.cxx $(UTIL_SRC)/*.hh $(UTIL_SRC)/*.h \
		$(NM_SRC)/*.cxx $(NM_SRC)/*.hh $(NM_SRC)/*.c $(NM_SRC)/*.h

tarfile:
	@echo Build a compressed tar-file in the $HOME directory
	- cd $(HOME); mv L1d3.tar.gz L1d3.tar.gz.old
	cd $(HOME); tar --exclude ".svn" --exclude "*~" -cvf L1d2.tar cfcfd/app/L1d3
	cd $(HOME); gzip L1d3.tar

#-----------------------------------------------------------

l1d3.a : $(OBJECTS)
ifneq ($(strip $(wildcard l1d3.a)),)
	- rm l1d3.a
endif
	$(AR) rv l1d3.a $(OBJECTS)
	- ranlib l1d3.a

$(LIBGAS) \
$(GAS)/build/_gaspy.so \
$(GAS)/build/_gaspy.dll \
$(GAS)/build/_gaspy.dylib \
$(GAS)/build/gaspy.py :
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new


gasfile:
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) \
			      INSTALL_DIR=$(INSTALL_DIR) install

$(LIBNM) \
$(NM)/build/_libnm.so \
$(NM)/build/_libnm.dll \
$(NM)/build/_libnm.dylib \
$(NM)/build/libnm.py:
	cd $(NM)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new

$(PLOT)/build/mb_cont.exe:
	cd $(PLOT)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR)

$(LIBINIPARSER):
	cd $(UTIL)/iniparser; make

$(LIBUTIL):
	cd $(UTIL)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new


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
$(LIBLUA) $(LUA_INCLUDE_DIR):
	cd $(LUA); make $(PLAT) local

$(LIBZLIB):
	cd $(ZLIB); make libz.a

#-----------------------------------------------------------

# -----------------
# Executable codes.
# -----------------

l1d.exe : l1d.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) l1d.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o l1d.exe

l_prep.exe : l_prep.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) l_prep.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o l_prep.exe

l_post.exe : l_post.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) l_post.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o l_post.exe

l_hist.exe : l_hist.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER)  $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) l_hist.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o l_hist.exe

sptime.exe : sptime.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) sptime.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o sptime.exe

piston.exe : piston.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) piston.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o piston.exe

l_trim.exe : l_trim.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) l_trim.o l1d3.a $(LIBGAS) $(LIBNM) $(LIBUTIL) $(LIBINIPARSER) \
		$(LIBLUA) $(LIBZLIB) $(LLIB) -o l_trim.exe

# Python files.

l_script.py : $(SRC)/l_script.py
	cp $(SRC)/l_script.py .
	chmod +x l_script.py

gaspy.py : $(GAS)/build/gaspy.py
	cp $(GAS)/build/gaspy.py .

_gaspy.so : $(GAS)/build/_gaspy.so
	cp $(GAS)/build/_gaspy.so .

_gaspy.dll : $(GAS)/build/_gaspy.dll
	cp $(GAS)/build/_gaspy.dll .

_gaspy.dylib : $(GAS)/build/_gaspy.dylib
	cp $(GAS)/build/_gaspy.dylib .

#----------------------------------------------------------

# -------------------------
# Object file dependencies.
# -------------------------

l1d.o : $(SRC)/l1d.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l1d.cxx -o l1d.o

l_kernel.o  : $(SRC)/l_kernel.cxx $(SRC)/l_kernel.hh $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_kernel.cxx -o l_kernel.o

l_bc.o    : $(SRC)/l_bc.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_bc.cxx -o l_bc.o

l_adapt.o    : $(SRC)/l_adapt.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_adapt.cxx -o l_adapt.o

l_io.o    : $(SRC)/l_io.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_io.cxx -o l_io.o

l_tstep.o : $(SRC)/l_tstep.cxx $(SRC)/l1d.hh $(NM_SRC)/qd_power.h $(NM_SRC)/qd_log10.h  $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_tstep.cxx -o l_tstep.o

l_misc.o  : $(SRC)/l_misc.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_misc.cxx -o l_misc.o

l_rivp.o  : $(SRC)/l_rivp.cxx $(SRC)/l1d.hh $(NM_SRC)/qd_power.h  $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_rivp.cxx -o l_rivp.o

roberts.o : $(NM_SRC)/roberts.c $(NM_SRC)/roberts.h
	$(COMPILE) $(CFLAG) $(NM_SRC)/roberts.c -o roberts.o

l_prep.o  : $(SRC)/l_prep.cxx $(SRC)/l1d.hh $(NM_SRC)/roberts.h $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_prep.cxx -o l_prep.o

l_post.o  : $(SRC)/l_post.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_post.cxx -o l_post.o

l_hist.o  : $(SRC)/l_hist.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -DNSPECD=$(NSPECD) -I$(LUA_INCLUDE_DIR) \
		$(SRC)/l_hist.cxx -o l_hist.o

sptime.o  : $(SRC)/sptime.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/sptime.cxx -o sptime.o

piston.o  : $(SRC)/piston.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/piston.cxx -o piston.o

l_trim.o  : $(SRC)/l_trim.cxx $(SRC)/l1d.hh $(LUA_INCLUDE_DIR)
	$(CXX) $(CXXFLAG) -DNDIM=$(NDIM) -I$(LUA_INCLUDE_DIR) $(SRC)/l_trim.cxx -o l_trim.o

# -------------------- end of makefile -------------------------------
