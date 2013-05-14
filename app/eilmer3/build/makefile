# \file    makefile
# \ingroup eilmer3
# \brief   Makefile for the Navier-Stokes program collection "eilmer3".
# \author  PA Jacobs
# \version 3.x ripped off from the mbcns2 makefile
#
# \usage   make
#          make TARGET=for_intel_mpi INSTALL_DIR=/work1/e4pjacob/e3bin
#          make INSTALL_DIR=/work1/e4pjacob/e3bin install
#          make clean
#          make tags
#          make tarfile
#----------------------------------------------------------------------

TARGET ?= for_gnu
MEMORY_SIZE ?= large
WITH_IMPLICIT ?= 0
INSTALL_DIR ?= $(HOME)/e3bin

ifeq ($(TARGET), for_openmpi)
    WITH_MPI := 1
endif

ifeq ($(TARGET), for_clang_openmpi)
    WITH_MPI := 1
endif

ifeq ($(TARGET), for_intel_mpi)
    WITH_MPI := 1
endif

ifeq ($(TARGET), for_intel_mpi_profile)
    WITH_MPI := 1
endif

ifeq ($(TARGET), for_intel_mpi_trace)
    WITH_MPI := 1
endif

ifeq ($(TARGET), for_mac_openmpi)
    WITH_MPI := 1
endif

# We seem to be having trouble with optimizing compilers.
# This optimization level may be over-ridden at the command line.
OPT ?= -O2

include ../../../lib/util/source/systems.mk
CFLAG += -DWITH_IMPLICIT=$(WITH_IMPLICIT)
CXXFLAG += -DWITH_IMPLICIT=$(WITH_IMPLICIT)
LLIB += -lpopt -ldl -L/opt/local/lib

# Test if we are asking for debugging information for a single cell
ifdef PRNT_DBG_CELL
    CXXFLAG += -DPRNT_DBG_CELL
endif
ifdef DBG_BLK
    CXXFLAG += -DDBG_BLK=$(DBG_BLK)
endif
ifdef DBG_CELL_I
    CXXFLAG += -DDBG_CELL_I=$(DBG_CELL_I)
endif
ifdef DBG_CELL_J
    CXXFLAG += -DDBG_CELL_J=$(DBG_CELL_J)
endif
ifdef DBG_CELL_K
    CXXFLAG += -DDBG_CELL_I=$(DBG_CELL_K)
endif

SRC             := ../source
LIBE3           := libe3.a
UTIL            := ../../../lib/util
UTIL_SRC        := $(UTIL)/source
LIBUTIL         := $(UTIL)/build/libutil.a
GEOM2           := ../../../lib/geometry2
GEOM2_SRC       := $(GEOM2)/source
LIBGEOM2        := $(GEOM2)/build/libgeom2.a
PLOT_SRC        := ../../../lib/plot/source
NM              := ../../../lib/nm
NM_SRC          := $(NM)/source
LIBNM           := $(NM)/build/libnm.a
GAS             := ../../../lib/gas
LIBGAS          := $(GAS)/build/libgas.a
RAD             := ../../../lib/radiation
RAD_SRC         := $(RAD)/source
LIBRAD          := $(RAD)/build/librad.a
LIBINIPARSER    := $(UTIL)/iniparser/libiniparser.a
LUA             := ../../../extern/lua-5.1.4
LIBLUA          := $(LUA)/lib/liblua.a
LUA_INCLUDE_DIR := $(LUA)/include
ZLIB            := ../../../extern/zlib-1.2.3
LIBZLIB         := $(ZLIB)/libz.a

#------------------------------------------------------------

E3_EXE_FILES := e3shared.exe e3console.tcl e3rad.exe
E3_MPI_EXE_FILES := e3mpi.exe
ifeq ($(WITH_MPI), 1)
    E3_EXE_FILES += $(E3_MPI_EXE_FILES)
endif
E3_UTILS := ../utilities/p2e.py ../utilities/import_grid.py \
            ../utilities/tecplot_feblock.py ../utilities/csv_util.py \
            ../utilities/laura2vtk.py ../utilities/import_laura_plot3d.py \
            ../utilities/shock1d.py ../utilities/e3loadbalance.py

E3_OBJECTS_1 := e3shared.o e3rad.o
E3_OBJECTS_1_MPI := e3mpi.o

E3_OBJECTS_2 := kernel.o \
	block.o \
	block_geometry.o \
	block_moving_grid.o \
	block_io.o \
	block_bgk.o \
	block_filter.o \
	exch2d.o \
	exch3d.o \
	bc.o \
	bc_catalytic.o \
	bc_adjacent.o \
	bc_supersonic_in.o \
	bc_extrapolate_out.o \
	bc_shock_fitting_in.o \
	bc_slip_wall.o \
	bc_adiabatic.o \
	bc_fixed_t.o \
	bc_sliding_t.o \
	bc_subsonic_in.o \
	bc_transient_uniform.o \
	bc_static_profile.o \
	bc_fixed_p_out.o \
	bc_ablating.o \
	bc_surface_energy_balance.o \
	bc_user_defined.o \
	bc_menter_correction.o \
	bc_fstc.o \
	init.o \
	block_invs.o \
	visc.o \
	visc3D.o \
	baldwin_lomax.o \
	diffusion.o \
	piston.o \
	rivp.o \
	efm.o \
	ausmdv.o \
	adaptive_flux.o \
        ausm_plus_up.o \
	hlle.o \
	flux_calc.o \
	ausm.o \
	cell.o \
	one_d_interp.o \
	flex_cell.o \
	time_to_go.o \
	c-flow-condition.o \
	lua_service.o \
	lua_service_for_e3.o \
	dbc_assert.o \
	radiation_transport.o \
	implicit.o \
	cell_finder.o \
	mersenne.o \
	ray_tracing_pieces.o \
	bgk.o 

ifeq ($(LES_DATA), 1)
    E3_OBJECTS_2 += les.o
endif

E3_OBJECTS_2_MPI := exch_mpi.o

PY_FILES = e3prep.py e3post.py turbo_post.py cgns_grid.py e3cgns.py e3history.py \
	e3_block.py e3_render.py e3_grid.py e3_flow.py bc_defs.py flux_dict.py \
	libprep3.py e3_defs.py e3_bc_util.py

# Assume Unix/Linux
LOADABLE_MODULE := _libprep3.so
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LOADABLE_MODULE := _libprep3.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LOADABLE_MODULE := _libprep3.dll
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LOADABLE_MODULE := _libprep3.dylib
endif
# The only place the loadable module is used is in e3prep.py.
PY_FILES += $(LOADABLE_MODULE)

#------------------------------------------------------------

build_new : $(E3_EXE_FILES) $(PY_FILES) 
	@echo "Built with C and C++ compilers "
	- $(COMPILE) --version
	- $(CXXCOMPILE) --version
	- @echo "CFLAG= $(CFLAG)"
	- @echo "CXXFLAG= $(CXXFLAG)"
	@echo "--------------------------------"
	@echo "Eilmer3 should be now up-to-date"
	@echo "--------------------------------"

install : $(INSTALL_DIR) $(E3_EXE_FILES) $(PY_FILES) $(E3_UTILS)
	@echo Copying executable files to $(INSTALL_DIR)
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd $(RAD)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd $(GEOM2)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	- cp -f $(E3_EXE_FILES) $(INSTALL_DIR)
	- cp -f $(PY_FILES) $(INSTALL_DIR)
	- cp -f $(E3_UTILS) $(INSTALL_DIR)
	- chmod +x $(INSTALL_DIR)/import_laura_plot3d.py
	- chmod +x $(INSTALL_DIR)/shock1d.py
	- chmod +x $(INSTALL_DIR)/e3loadbalance.py
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
	- cp $(INSTALL_DIR)/_libprep3.dylib $(INSTALL_DIR)/_libprep3.so
endif
	cd ../../../lib/; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

clean :
	@echo Deleting executable files
	- rm $(E3_EXE_FILES) $(E3_MPI_EXE_FILES)
	- rm $(PY_FILES)
	@echo Deleting object files
	- rm $(E3_OBJECTS_1) $(E3_OBJECTS_2) $(E3_OBJECTS_1_MPI) $(E3_OBJECTS_2_MPI)
	- rm *.pyc
	- rm libprep3_wrap.o
	- rm $(SRC)/libprep3.py
	- rm $(SRC)/libprep3_wrap.cxx
	@echo Deleting archive/library of object files
	- rm $(LIBE3)
	@echo Cleaning up other areas
	- cd $(GAS)/build; make clean
	- cd $(NM)/build; make clean
	- cd $(GEOM2)/build; make clean
	- cd $(UTIL)/build; make clean
	- cd $(UTIL)/iniparser; make clean
	- cd $(LUA); make clean
	- cd $(RAD)/build; make clean

tarfile :
	- cd $(HOME); mv eilmer3.tar.gz eilmer3.tar.gz.old
	cd $(HOME); tar --exclude ".svn" --exclude "*~" -cvf eilmer3.tar cfcfd/app/eilmer3
	cd $(HOME); gzip eilmer3.tar

tags:
	cd ../source; etags ./*.cxx ./*.hh \
		$(UTIL_SRC)/*.cxx $(UTIL_SRC)/*.c $(UTIL_SRC)/*.hh \
		$(NM_SRC)/*.cxx $(NM_SRC)/*.hh $(NM_SRC)/*.c $(NM_SRC)/*.h \
		$(GEOM2_SRC)/*.cxx $(GEOM2_SRC)/*.hh \
		$(GAS)/models/*.cxx $(GAS)/models/*.hh \
		$(GAS)/kinetics/*.cxx $(GAS)/models/*.hh \
		$(UTIL)/iniparser/src/*.c $(UTIL)/iniparser/src/*.h \
		$(RAD_SRC)/*.cxx $(RAD_SRC)/*.hh

e3rad : 
	@echo Installing a fresh copy of e3rad.exe with all dependent components rebuilt
	- rm $(RAD)/build/librad.a $(RAD)/build/parade.o
	make TARGET=$(TARGET) OPT=$(OPT) $(LIBRAD) 
	- rm e3rad.exe e3rad.o radiation_transport.o cell_finder.o
	make TARGET=$(TARGET) OPT=$(OPT) e3rad.exe
	- cp -f e3rad.exe $(INSTALL_DIR)


#---------------------------------------------------------------------

# -----------------------
# Library of code modules
# -----------------------

$(LIBE3) : $(E3_OBJECTS_2)
ifneq ($(strip $(wildcard $(LIBE3))),)
	- rm $(LIBE3)
endif
	$(AR) rv $(LIBE3) $(E3_OBJECTS_2)
	- ranlib $(LIBE3)

$(LIBRAD):
	cd $(RAD)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new

$(LIBGAS):
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new

gasfile:
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) \
			      INSTALL_DIR=$(INSTALL_DIR) install

$(LIBNM):
	cd $(NM)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new

$(LIBGEOM2):
	cd $(GEOM2)/build; make TARGET=$(TARGET) OPT=$(OPT) build_new

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
$(LIBLUA):
	cd $(LUA); make $(PLAT) local

$(LIBZLIB):
	cd $(ZLIB); make libz.a


bc_defs.py : $(SRC)/bc_defs.py
	cp $(SRC)/bc_defs.py .

e3_defs.py : $(SRC)/e3_defs.py
	cp $(SRC)/e3_defs.py .

e3_block.py : $(SRC)/e3_block.py
	cp $(SRC)/e3_block.py .

e3_grid.py : $(SRC)/e3_grid.py
	cp $(SRC)/e3_grid.py .

e3_flow.py : $(SRC)/e3_flow.py
	cp $(SRC)/e3_flow.py .

e3_render.py : $(SRC)/e3_render.py
	cp $(SRC)/e3_render.py .

e3_bc_util.py : $(SRC)/e3_bc_util.py
	cp $(SRC)/e3_bc_util.py .

flux_dict.py : $(SRC)/flux_dict.py
	cp $(SRC)/flux_dict.py .

#-----------------------------------------------------------------------
# Loadable library for Python.
# This is a combination of the libgas2 and libgeom2 modules.

# Dynamic link library for MINGW32 environment on MS-Windows
_libprep3.dll : libprep3_wrap.o $(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBPYTHON) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) -shared -o _libprep3.dll libprep3_wrap.o \
		-L$(PYTHON_DIR)/libs \
		$(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBPYTHON) $(LIBLUA) $(LIBZLIB) $(LLIB)

# Shared object for Unix/Linux
_libprep3.so : libprep3_wrap.o $(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) -shared -o _libprep3.so libprep3_wrap.o \
		$(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBLUA) $(LIBZLIB) $(LLIB)

# Dynamic library for Mac OS X
_libprep3.dylib : libprep3_wrap.o $(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(LFLAG) -shared -o _libprep3.dylib libprep3_wrap.o \
		$(LIBE3) $(LIBGEOM2) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBLUA) $(LIBZLIB) $(LLIB) -framework Python

# The SWIG interface file and corresponding documentation
# Normally, we would not try to rebuild this file because 
# some computers will not have swig.
# There are many more dependencies for libprep3 but 
# it is probably easier to make clean before a make.

$(SRC)/libprep3_wrap.cxx $(SRC)/libprep3.py : $(SRC)/libprep3.i \
		$(GAS)/models/libgas.i $(GEOM2_SRC)/libgeom2.i $(RAD_SRC)/radpy.i
	swig -python -c++ -w512 $(SRC)/libprep3.i

libprep3.py : $(SRC)/libprep3.py
	cp $(SRC)/libprep3.py .

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
libprep3_wrap.o : $(SRC)/libprep3_wrap.cxx
	$(CXXCOMPILE) $(CXXFLAG3) -I$(PYTHON_INCLUDE_DIR) -I$(LUA_INCLUDE_DIR) \
		-I$(GAS)/models -I$(GAS)/kinetics -I$(GEOM2_SRC) \
		$(SRC)/libprep3_wrap.cxx -I$(RAD_SRC)

#----------------------------------------------------------------------

# ----------------
# Executable files
# ----------------

e3prep.py : $(SRC)/e3prep.py
	cp $(SRC)/e3prep.py .
	chmod +x e3prep.py

e3post.py : $(SRC)/e3post.py
	cp $(SRC)/e3post.py .
	chmod +x e3post.py

turbo_post.py : $(SRC)/turbo_post.py
	cp $(SRC)/turbo_post.py .
	chmod +x turbo_post.py

cgns_grid.py : $(SRC)/cgns_grid.py
	cp $(SRC)/cgns_grid.py .
	chmod +x cgns_grid.py

e3cgns.py : $(SRC)/e3cgns.py
	cp $(SRC)/e3cgns.py .
	chmod +x e3cgns.py

e3history.py : $(SRC)/e3history.py
	cp $(SRC)/e3history.py .
	chmod +x e3history.py

e3console.tcl: $(SRC)/e3console.tcl
	cp $(SRC)/e3console.tcl .
	chmod +x e3console.tcl

e3shared.exe : e3shared.o $(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(PCA) $(LFLAG) e3shared.o \
		$(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB) $(LLIB) \
		-o e3shared.exe

e3mpi.exe : e3mpi.o $(E3_OBJECTS_2_MPI) \
		$(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(PCA) $(LFLAG) e3mpi.o $(E3_OBJECTS_2_MPI) \
		$(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LMPI) $(LIBLUA) $(LIBZLIB) $(LLIB) \
		-o e3mpi.exe

e3rad.exe : e3rad.o $(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB)
	$(CXXLINK) $(PCA) $(LFLAG) e3rad.o \
		$(LIBE3) $(LIBGAS) $(LIBRAD) $(LIBNM) $(LIBUTIL) $(LIBGEOM2) \
		$(LIBINIPARSER) $(LIBLUA) $(LIBZLIB) $(LLIB) \
		-o e3rad.exe

#---------------------------------------------------------------------
# Components that need to be compiled as C modules
#

time_to_go.o : $(UTIL_SRC)/time_to_go.c $(UTIL_SRC)/time_to_go.h
	$(COMPILE) $(CFLAG) $(UTIL_SRC)/time_to_go.c -o time_to_go.o

#---------------------------------------------------------------------
# Components that are to be compiled as C++ modules.
#

e3shared.o : $(SRC)/main.cxx $(SRC)/main.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh \
		$(SRC)/init.hh $(SRC)/exch2d.hh $(SRC)/visc.hh $(SRC)/piston.hh $(LIBLUA)
	$(CXXCOMPILE) $(PCA) $(CXXFLAG) -DSHARED -I/sw/include -I/opt/local/include -I$(LUA_INCLUDE_DIR) \
		$(SRC)/main.cxx -o e3shared.o

e3mpi.o : $(SRC)/main.cxx $(SRC)/main.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh \
		$(SRC)/init.hh $(SRC)/exch2d.hh $(SRC)/exch_mpi.hh \
		$(SRC)/visc.hh $(SRC)/piston.hh $(LIBLUA)
	$(CXXCOMPILE) $(PCA) $(CXXFLAG) $(CFLAG_MPI) -D_MPI -I/sw/include -I/opt/local/include -I$(LUA_INCLUDE_DIR) \
		$(SRC)/main.cxx -o e3mpi.o

e3rad.o : $(SRC)/main.cxx $(SRC)/main.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh \
		$(SRC)/init.hh $(SRC)/exch2d.hh $(SRC)/visc.hh $(SRC)/piston.hh \
		$(SRC)/radiation_transport.hh $(LIBLUA)
	$(CXXCOMPILE) $(PCA) $(CXXFLAG) -DSHARED -DE3RAD -I/sw/include -I/opt/local/include -I$(LUA_INCLUDE_DIR) \
		$(SRC)/main.cxx -o e3rad.o

exch2d.o : $(SRC)/exch2d.cxx $(SRC)/exch3d.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/exch2d.cxx -o exch2d.o

exch3d.o : $(SRC)/exch3d.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/exch3d.cxx -o exch3d.o

exch_mpi.o : $(SRC)/exch_mpi.cxx $(SRC)/exch_mpi.hh $(SRC)/block.hh \
		$(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(CFLAG_MPI) -D_MPI $(SRC)/exch_mpi.cxx -o exch_mpi.o

init.o : $(SRC)/init.cxx $(SRC)/init.hh $(SRC)/kernel.hh $(SRC)/block.hh \
		$(SRC)/bc.hh $(SRC)/bc_defs.hh $(SRC)/diffusion.hh $(SRC)/visc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/init.cxx -o init.o

bc.o : $(SRC)/bc.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) -I$(ZLIB) $(SRC)/bc.cxx -o bc.o

bc_catalytic.o : $(SRC)/bc_catalytic.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_catalytic.cxx -o bc_catalytic.o

bc_adjacent.o : $(SRC)/bc_adjacent.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_adjacent.cxx -o bc_adjacent.o

bc_supersonic_in.o : $(SRC)/bc_supersonic_in.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh $(SRC)/bc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_supersonic_in.cxx -o bc_supersonic_in.o

bc_extrapolate_out.o : $(SRC)/bc_extrapolate_out.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh $(SRC)/bc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_extrapolate_out.cxx -o bc_extrapolate_out.o

bc_shock_fitting_in.o : $(SRC)/bc_shock_fitting_in.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh $(SRC)/bc.hh $(SRC)/one_d_interp.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_shock_fitting_in.cxx -o bc_shock_fitting_in.o

bc_slip_wall.o : $(SRC)/bc_slip_wall.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_slip_wall.cxx -o bc_slip_wall.o

bc_adiabatic.o : $(SRC)/bc_adiabatic.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_adiabatic.cxx -o bc_adiabatic.o

bc_fixed_t.o : $(SRC)/bc_fixed_t.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_fixed_t.cxx -o bc_fixed_t.o

bc_fstc.o : $(SRC)/bc_fstc.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_fstc.cxx -o bc_fstc.o

bc_sliding_t.o : $(SRC)/bc_sliding_t.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_sliding_t.cxx -o bc_sliding_t.o

bc_subsonic_in.o : $(SRC)/bc_subsonic_in.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_subsonic_in.cxx -o bc_subsonic_in.o

bc_transient_uniform.o : $(SRC)/bc_transient_uniform.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_transient_uniform.cxx -o bc_transient_uniform.o

bc_static_profile.o : $(SRC)/bc_static_profile.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_static_profile.cxx -o bc_static_profile.o

bc_fixed_p_out.o : $(SRC)/bc_fixed_p_out.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_fixed_p_out.cxx -o bc_fixed_p_out.o

bc_ablating.o : $(SRC)/bc_ablating.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_ablating.cxx -o bc_ablating.o

bc_surface_energy_balance.o : $(SRC)/bc_surface_energy_balance.cxx $(SRC)/block.hh \
    $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_surface_energy_balance.cxx \
	    -o bc_surface_energy_balance.o

bc_user_defined.o : $(SRC)/bc_user_defined.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) \
    $(NM_SRC)/zero_finders.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_user_defined.cxx -o bc_user_defined.o

bc_menter_correction.o : $(SRC)/bc_menter_correction.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/bc_menter_correction.cxx -o bc_menter_correction.o

block.o : $(SRC)/block.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) $(LIBZLIB)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) -I$(ZLIB) $(SRC)/block.cxx -o block.o

block_filter.o : $(SRC)/block_filter.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/block_filter.cxx -o block_filter.o

block_io.o : $(SRC)/block_io.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) $(LIBZLIB)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) -I$(ZLIB) $(SRC)/block_io.cxx -o block_io.o

block_geometry.o : $(SRC)/block_geometry.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/block_geometry.cxx -o block_geometry.o

block_moving_grid.o : $(SRC)/block_moving_grid.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) $(LIBZLIB)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) -I$(ZLIB) $(SRC)/block_moving_grid.cxx -o block_moving_grid.o

block_bgk.o : $(SRC)/block_bgk.cxx $(SRC)/cell.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA) $(LIBZLIB)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) -I$(ZLIB) $(SRC)/block_bgk.cxx -o block_bgk.o

kernel.o : $(SRC)/kernel.cxx $(SRC)/cell.hh $(SRC)/kernel.hh $(SRC)/block.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/kernel.cxx -o kernel.o

block_invs.o : $(SRC)/block_invs.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/block_invs.cxx -o block_invs.o

visc.o : $(SRC)/visc.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(SRC)/bc_defs.hh \
		$(SRC)/diffusion.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/visc.cxx -o visc.o

visc3D.o : $(SRC)/visc3D.cxx $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(SRC)/bc_defs.hh \
		$(SRC)/diffusion.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/visc3D.cxx -o visc3D.o

$(SRC)/visc3D.cxx : $(SRC)/visc3D.cxx.m4
	m4 $(SRC)/visc3D.cxx.m4 > $(SRC)/visc3D.cxx

baldwin_lomax.o : $(SRC)/baldwin_lomax.cxx $(SRC)/block.hh $(SRC)/kernel.hh \
		$(SRC)/bc.hh $(SRC)/bc_defs.hh \
		$(SRC)/diffusion.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/baldwin_lomax.cxx -o baldwin_lomax.o

diffusion.o : $(SRC)/diffusion.cxx $(SRC)/diffusion.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/diffusion.cxx -o diffusion.o

piston.o : $(SRC)/piston.cxx $(SRC)/piston.hh $(SRC)/block.hh $(SRC)/kernel.hh $(SRC)/bc.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/piston.cxx -o piston.o

rivp.o : $(SRC)/rivp.cxx $(SRC)/flux_calc.hh $(NM_SRC)/qd_power.h
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/rivp.cxx -o rivp.o

efm.o : $(SRC)/efm.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/efm.cxx -o efm.o

ausmdv.o : $(SRC)/ausmdv.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/ausmdv.cxx -o ausmdv.o

adaptive_flux.o : $(SRC)/adaptive_flux.cxx $(SRC)/flux_calc.hh $(SRC)/kernel.hh $(SRC)/cell.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/adaptive_flux.cxx -o adaptive_flux.o

ausm_plus_up.o : $(SRC)/ausm_plus_up.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/ausm_plus_up.cxx -o ausm_plus_up.o

hlle.o : $(SRC)/hlle.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/hlle.cxx -o hlle.o

flux_calc.o : $(SRC)/flux_calc.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/flux_calc.cxx -o flux_calc.o

ausm.o : $(SRC)/ausm.cxx $(SRC)/flux_calc.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/ausm.cxx -o ausm.o

cell.o : $(SRC)/cell.cxx $(SRC)/flux_calc.hh $(SRC)/cell.hh $(SRC)/kernel.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/cell.cxx -o cell.o

one_d_interp.o : $(SRC)/one_d_interp.cxx $(SRC)/cell.hh $(SRC)/one_d_interp.hh $(SRC)/kernel.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/one_d_interp.cxx -o one_d_interp.o

flex_cell.o: $(SRC)/flex_cell.cxx $(SRC)/flex_cell.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/flex_cell.cxx -o flex_cell.o

c-flow-condition.o : $(SRC)/c-flow-condition.cxx $(SRC)/c-flow-condition.hh $(GAS)/models/gas_data.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/c-flow-condition.cxx

lua_service.o : $(UTIL_SRC)/lua_service.hh $(UTIL_SRC)/lua_service.cxx
	$(CXXCOMPILE) $(CXXFLAG) $(UTIL_SRC)/lua_service.cxx -I$(LUA_INCLUDE_DIR) -I$(ZLIB)

lua_service_for_e3.o : $(SRC)/lua_service_for_e3.cxx $(SRC)/lua_service_for_e3.hh
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/lua_service_for_e3.cxx

dbc_assert.o : $(UTIL_SRC)/dbc_assert.cxx $(UTIL_SRC)/dbc_assert.hh
	$(CXXCOMPILE) $(CXXFLAG) $(UTIL_SRC)/dbc_assert.cxx

radiation_transport.o : $(SRC)/radiation_transport.hh $(SRC)/radiation_transport.cxx \
		$(LUA_INCLUDE_DIR)
	$(CXXCOMPILE) $(PCA) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/radiation_transport.cxx \
		-I$(LUA_INCLUDE_DIR)

implicit.o : $(SRC)/implicit.cxx $(SRC)/implicit.hh $(LIBLUA)
	$(CXXCOMPILE) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/implicit.cxx -o implicit.o

cell_finder.o : $(SRC)/cell_finder.cxx $(SRC)/cell_finder.hh
	$(CXXCOMPILE) $(PCA) $(CXXFLAG) -I$(LUA_INCLUDE_DIR) $(SRC)/cell_finder.cxx -o cell_finder.o

mersenne.o : $(UTIL_SRC)/mersenne.cpp $(UTIL_SRC)/randomc.h
	$(CXXCOMPILE) $(CXXFLAG) $(UTIL_SRC)/mersenne.cpp

ray_tracing_pieces.o : $(SRC)/ray_tracing_pieces.cxx $(SRC)/ray_tracing_pieces.hh 
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/ray_tracing_pieces.cxx -I$(LUA_INCLUDE_DIR)

bgk.o : $(SRC)/bgk.cxx $(SRC)/bgk.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/bgk.cxx -o bgk.o

#
#---------------------------- end of makefile -----------------------------

