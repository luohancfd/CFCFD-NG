# \file    systems.mk
# \ingroup util
# \brief   Makefile component to customise compiler and linker options.
# \author  PA Jacobs
# \version 07-Oct-04: ripped out of the old makefile

SYSTEM := $(shell uname -s)
ARCH   := $(shell uname -m)



ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    TCL_DIR := /c/Tcl
    TCL_INCLUDE_DIR := $(TCL_DIR)/include
else
    # Assume that we are building in a Unix/Linux environment.
    # Debian
    TCL_INCLUDE_DIR := $(dir $(firstword $(wildcard /usr/include/tcl8*/tcl.h)))
    # Unix/Linux with ActiveState Tcl installed.
    ifeq ($(strip $(TCL_INCLUDE_DIR)),)
        TCL_INCLUDE_DIR := $(dir $(wildcard /usr/local/ActiveTcl/include/tcl.h))
    endif
    # Hopefully, there is a Tcl distribution in the usual place.
    ifeq ($(strip $(TCL_INCLUDE_DIR)),)
        TCL_INCLUDE_DIR := $(dir $(wildcard /usr/include/tcl.h))
    endif
endif

# Take a guess at where to find the Python header files.
PYTHON_MAJOR := $(shell python -c 'import sys; print sys.version_info[0]')
PYTHON_MINOR := $(shell python -c 'import sys; print sys.version_info[1]')
PYTHON_PREFIX := $(shell python -c 'import sys; print sys.prefix')
PYTHON_NN := python$(PYTHON_MAJOR)$(PYTHON_MINOR)
PYTHON_NdN := python$(PYTHON_MAJOR).$(PYTHON_MINOR)
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    PYTHON_DIR := /c/Python$(PYTHON_MAJOR)$(PYTHON_MINOR)
    LIBPYTHON := $(PYTHON_DIR)/libs/lib$(PYTHON_NN).a
    PYTHON_INCLUDE_DIR := $(PYTHON_DIR)/include
else
    # Assume that we are building in a Unix/Linux environment.
    # Fedora, Debian, MacOS-X, Cygwin
    PYTHON_INCLUDE_DIR := $(PYTHON_PREFIX)/include/python$(PYTHON_MAJOR).$(PYTHON_MINOR)
    PYTHON_BIN_DIR := $(PYTHON_PREFIX)/bin
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        LIBPYTHON := $(PYTHON_PREFIX)/lib/$(PYTHON_NdN)/config/lib$(PYTHON_NdN).dll.a
    endif
endif

# Take a guess at where to find the CGNS header and library files.
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    CGNS_LIB := /c/cgns/libcgns.a
    CGNS_INCLUDE_DIR := /c/cgns/include
else
    # Assume that we are building in a Unix/Linux environment.
    # Andrew has installed CGNS in /opt on gemini1 and blackhole.
    CGNS_LIB := $(wildcard /opt/cgns/lib/libcgns.a)
    CGNS_INCLUDE_DIR := $(wildcard /opt/cgns/include)
    # A standard install would be into /usr/local
    ifeq ($(strip $(CGNS_LIB)),)
        CGNS_LIB := $(wildcard /usr/local/lib/libcgns.a)
        CGNS_INCLUDE_DIR := $(wildcard /usr/local/include)
    endif
endif

#----------------------------------------------------------------------

TARGET ?= for_gnu

LMPI :=
PCA  :=
OPT ?= -O
AR := ar

# ------------------------------------
#            GNU Compilers
# ------------------------------------

# When using GNU compilers, user may pass in an MARCH flag
# NOTE: from the gcc man page --
#
#   specifying -march=cpu-type implies -mtune=cpu-type.
#
# Some examples:
#   Intel Xeon (64-bit)  -->  MARCH=nocona
#   AMD Opteron          -->  MARCH=opteron
#   Intel Core 2 Duo     -->  MARCH=core2
#

ifneq ($(MARCH),)
	MARCH_FLAG := -march=$(MARCH)
endif

ifeq ($(TARGET), for_clang)
    # UNIX/Linux workstation with the Clang C, C++ compiler
    COMPILE := clang 
    LINK    := clang
    CXX     := clang
    CXXLINK := clang 
    # Unix/Linux is default
    CFLAG   := -c $(OPT)
    LFLAG   := $(OPT)
    CXXFLAG := -c $(OPT) -std=c++11
    LLIB := -lstdc++ -lm
endif

ifeq ($(TARGET), for_gnu)
    # UNIX/Linux workstation with the default GNU C compiler
    # Don't specify the processor architecture.
    COMPILE := gcc 
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    # Unix/Linux is default
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic -finline-limit=2400 $(MARCH_FLAG)
    LFLAG   := $(OPT) -fPIC -finline-limit=2400 $(MARCH_FLAG)
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -pedantic $(MARCH_FLAG)
    ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
        # MINGW32 environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic $(MARCH_FLAG)
    endif
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic $(MARCH_FLAG)
    endif
    LLIB := -lm
    ifeq ($(WITH_SPRADIAN), 1)
        # Define the Fortran 90 compiler.
        F90 := gfortran
        # Compile without vector optimization for now.
        F90FLAG := -m64 -c -O0 -fPIC
        F90LFLAG := -m64 -lstdc++ -gnofor_main
        FLINK := -lgfortran
    endif
endif

ifeq ($(TARGET), for_gnu_gcc4)
    # UNIX/Linux workstation with the GNU C compiler version 4.x
    # Don't specify the processor architecture.
    COMPILE := gcc-4 
    LINK    := gcc-4
    CXX     := g++-4
    CXXLINK := g++-4
    # Unix/Linux is default
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic -finline-limit=2400 $(MARCH_FLAG)
    LFLAG   := $(OPT) -fPIC -finline-limit=2400 $(MARCH_FLAG)
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -pedantic $(MARCH_FLAG)
    ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
        # MINGW32 environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic $(MARCH_FLAG)
    endif
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic $(MARCH_FLAG)
    endif
    LLIB := -lm
endif

ifeq ($(TARGET), for_gnu_debug)
    # UNIX/Linux workstation with the GNU C compiler and debug+profiling
    COMPILE := gcc 
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    # Unix/Linux is default
    CFLAG   := -c -fPIC -W -Wall -pedantic -ggdb $(MARCH_FLAG)
    LFLAG   := -fPIC -pedantic -ggdb $(MARCH_FLAG)
    CXXFLAG := -c -fPIC -std=c++0x -Wall -pedantic -ggdb $(MARCH_FLAG)
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG   := -c -W -Wall -pedantic -g $(MARCH_FLAG)
        LFLAG   := -g -Wl,stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c -std=c++0x -Wall -pedantic -g $(MARCH_FLAG)
    endif
    ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
        # MINGW32 environment on MS-Windows
        CFLAG   := -c -W -Wall -pedantic -g $(MARCH_FLAG)
        LFLAG   := -g -Wl,stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c -std=c++0x -Wall -pedantic -g $(MARCH_FLAG)
    endif
    CFLAG   += -ftrapping-math -fsignaling-nans -DDEBUG
    CXXFLAG += -ftrapping-math -fsignaling-nans -DDEBUG
    LFLAG   += -ftrapping-math -fsignaling-nans -DDEBUG
    LLIB    := -lm
    # Have removed -lefence in favour of valgrind. 
endif

ifeq ($(TARGET), for_gprof)
    # UNIX/Linux/Cygwin workstation with the GNU C compiler and with profiling.
    COMPILE := gcc 
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    # Unix/Linux is default
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic -finline-limit=2400 -pg $(MARCH_FLAG)
    LFLAG   := $(OPT) -fPIC -finline-limit=2400 -pg $(MARCH_FLAG)
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -pedantic -pg $(MARCH_FLAG)
    ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
        # MINGW32 environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic -pg $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 -pg $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic -pg $(MARCH_FLAG)
    endif
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG   := -c $(OPT) -W -Wall -pedantic -pg $(MARCH_FLAG)
        LFLAG   := $(OPT) -Wl,-stack=0x8000000 -pg $(MARCH_FLAG)
        CXXFLAG := -c $(OPT) -std=c++0x -Wall -pedantic -pg $(MARCH_FLAG)
    endif
    LLIB    := -lm 
endif

ifeq ($(TARGET), for_gnu_opteron)
    # Linux Opteron server using GNU compiler
    COMPILE := gcc
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic -mtune=opteron -finline-limit=2400 
    LFLAG   := $(OPT) -fPIC -finline-limit=2400 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -mtune=opteron -pedantic 
    LLIB    := -lm
endif

ifeq ($(TARGET), for_macports_gnu)
    # OpenMPI on Mac OS-X with MacPorts
    # PJ edits for Ingo's Mac 06-June-2012.
    # Note that this works on LinuxMint13 as well as the default for_gnu.
    COMPILE := gcc
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    CFLAG   := -c $(OPT) -fPIC 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC 
    LFLAG   :=  $(OPT) -fPIC -lstdc++
    LLIB    := -lm -lstdc++
endif

ifeq ($(TARGET), for_gnu_openmp)
    # Linux Opteron server using GNU compiler
    COMPILE := gcc
    LINK    := gcc
    CXX     := g++
    CXXLINK := g++
    # Unix/Linux
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic -mtune=opteron -finline-limit=2400 $(MARCH_FLAG)
    LFLAG   := $(OPT) -fPIC -finline-limit=2400 $(MARCH_FLAG)
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -mtune=opteron -pedantic $(MARCH_FLAG)
    LLIB    := -lm
    PCA     := -fopenmp
endif

# ------------------------------------
#          Portland Compilers
# ------------------------------------

ifeq ($(TARGET), for_pgi)
    # UNIX/Linux workstation with the Portland-Group C compiler
    COMPILE := pgcc 
    LINK    := pgcc
    CXX     := pgCC
    CXXLINK := pgCC
    # Unix/Linux
    CFLAG   := -c -fast -fPIC 
    LFLAG   := -fast -fPIC 
    LLIB    := -lm
    CXXFLAG := -c -fast -fPIC 
#    PCA     := -mp
endif

ifeq ($(TARGET), for_pgi_debug)
    # UNIX/Linux workstation with the Portland-Group C compiler
    COMPILE := pgcc
    LINK    := pgcc
    CXX     := pgCC
    CXXLINK := pgCC
    # Unix/Linux
    CFLAG   := -c -g -p -fPIC 
    LFLAG   := -g -fPIC 
    LLIB    := -lm
    CXXFLAG := -c -g -p -fPIC 
#    PCA     := -mp
endif

ifeq ($(TARGET), for_pgi_opteron)
    # UNIX/Linux workstation with the Portland-Group C compiler
    COMPILE := pgcc 
    LINK    := pgcc
    CXX     := pgCC
    CXXLINK := pgCC
    # Unix/Linux
    CFLAG   := -c -fast -fPIC -tp amd64 
    LFLAG   := -fast -fPIC -tp amd64 
    LLIB    := -lm 
    CXXFLAG := -c -fast -fPIC -tp amd64 
#    PCA     := -mp
endif

ifeq ($(TARGET), for_pgi_opteron_O2)
    # UNIX/Linux workstation with the Portland-Group C compiler
    COMPILE := pgcc 
    LINK    := pgcc
    CXX     := pgCC
    CXXLINK := pgCC
    # Unix/Linux
    CFLAG   := -c -O2 -fPIC -tp amd64 
    LFLAG   := -O2 -fPIC -tp amd64 
    LLIB    := -lm 
    CXXFLAG := -c -O2 -fPIC -tp amd64 
#    PCA     := -mp
endif

# -------------------------------------
#  Intel compilers
# -------------------------------------

ifeq ($(TARGET), for_intel)
    # Intel Compiler 11.0 
    COMPILE := icc
    LINK    := icc
    CXX     := icpc
    CXXLINK := icpc
    CFLAG   := -c -wd161 -fPIC -xhost -ipo -std=c++11
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -fPIC -ipo -std=c++11
    LLIB    := -lm
    AR      := xiar
endif

ifeq ($(TARGET), for_intel_xeon)
    # Intel Compiler 9.0 
    COMPILE := icc
    LINK    := icc
    CXX     := icpc
    CXXLINK := icpc
    CFLAG   := -c -wd161 -fPIC -mtune=pentium4 
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl,
    LLIB    := -lm
#    PCA     := -openmp
endif

ifeq ($(TARGET), for_intel_itanium_2)
    # Intel compiler on APAC with linking to mpi library
    COMPILE := icc
    LINK    := icc
    CXX     := icpc
    CXXLINK := icpc
    CFLAG   := -c -wd161 -mcpu=itanium2 -O3
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -O3
    LLIB    := -lm
    LMPI   := -lmpi
    PCA     := 
endif

ifeq ($(TARGET), for_intel_mpi)
    # Intel Compiler, optimized for use on the barrine.hpcu.uq.edu.au cluster.
    # PJ, 08-Sep-2010 (at SGI course on Intel X86_64 Application Development).  
    COMPILE := mpiicc
    LINK    := mpiicc
    CXX     := mpiicpc
    CXXLINK := mpiicpc
    CFLAG   := -c -wd161 -fPIC -xhost -ipo -std=c++11
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -fPIC -ipo -std=c++11
    LLIB    := -lm
    LMPI    := -lmpi
    AR      := xiar
endif

ifeq ($(TARGET), for_intel_mpi_profile)
    # Intel Compiler, instrumented for profiling (gprof).  
    # PJ, 08-Sep-2010 (at SGI course on Intel X86_64 Application Development).  
    COMPILE := mpiicc
    LINK    := mpiicc
    CXX     := mpiicpc
    CXXLINK := mpiicpc
    CFLAG   := -c -wd161 -fPIC -ipo -opt-report -opt-report-phase=ipo -p -std=c++11
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -fPIC -ipo -opt-report -opt-report-phase=ipo -p -std=c++11
    LLIB    := -lm
    LMPI    := -lmpi
    AR      := xiar
endif

ifeq ($(TARGET), for_intel_mpi_trace)
    # Intel Compiler, optimized plus instrumented for use with Intel traceanalyzer.
    # Not that the build process will break when it gets to linking the shared libraries
    # for loading with Python, however, that doesn't matter since we really only want
    # the e3mpi.exe executable file.
    # PJ, 10-Sep-2010 (at SGI course on Intel X86_64 Application Development).  
    COMPILE := mpiicc
    LINK    := mpiicc
    CXX     := mpiicpc
    CXXLINK := mpiicpc
    CFLAG   := -c -wd161 -fPIC -ipo -trace -std=c++11
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -fPIC -ipo -trace -std=c++11
    LLIB    := -lm -lVT
    LMPI    := -lmpi
    AR      := xiar
endif

ifeq ($(TARGET), for_intel_openmp)
    # Intel Compiler 11.0 
    COMPILE := icc
    LINK    := icc
    CXX     := icpc
    CXXLINK := icpc
    CFLAG   := -c -wd161 -fPIC -xhost -ipo -std=c++11
    CXXFLAG := $(CFLAG)
    LFLAG   := -Wl, -fPIC -ipo -std=c++11
    LLIB    := -lm
    AR      := xiar
    PCA     := -openmp
endif

# -------------------------------------
#    MPI wrappers
# -------------------------------------

ifeq ($(TARGET), for_lam)
    # LAM on Linux. 
    # This was the usual target for the MPI code, until mid-2007.
    # Setting the Lam compiler-specific environment variables
    # will allow setting of options specific to each compiler.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) 
    CXXFLAG := -c $(OPT) -std=c++0x  
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
    LMPI    := -lmpi
endif

ifeq ($(TARGET), for_lam_debug)
    # LAM on Linux with debug on.
    # Setting the Lam compiler-specific environment variables
    # will allow setting of options specific to each compiler.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c -g
    CXXFLAG := -c -g -std=c++0x 
    LFLAG   :=  -fPIC -g
    LLIB    := -lm
    LMPI    := -lmpi
endif

ifeq ($(TARGET), for_openmpi)
    # OpenMPI on Linux.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) -fPIC 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
endif

ifeq ($(TARGET), for_openmpi_debug)
    # OpenMPI on Linux.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) -fPIC -g
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -g
    LFLAG   :=  $(OPT) -fPIC -g
    LLIB    := -lm
endif

ifeq ($(TARGET), for_macports_openmpi)
    # OpenMPI on Mac OS-X with MacPorts
    # PJ edits for Ingo's Mac 06-June-2012.
    # Note (1) This works on LinuxMint13 as well as the usual for_openmpi.
    # Note (2) The case sensitivity seems to work since Dan's note below,
    #          at least within the xcode unix-like environment that Ingo uses.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) -fPIC 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC 
    LFLAG   :=  $(OPT) -fPIC -lstdc++
    LLIB    := -lm -lstdc++ -lmpi_cxx
endif

ifeq ($(TARGET), for_mac_openmpi)
    # OpenMPI on Mac OS X as used by Dan Potter, PJ believes...
    # HFS+ is not case-sensitive, so need to use mpic++ instead of mpiCC.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpic++
    CXXLINK := mpic++
    CFLAG   := -c $(OPT) -fPIC 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
endif

ifeq ($(TARGET), for_mpich)
    # mpich on Linux.
    # Setting the Lam compiler-specific environment variables
    # will allow setting of options specific to each compiler.
    COMPILE := mpicc
    LINK    := mpicc
    CXX     := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) 
    CXXFLAG := -c $(OPT) -std=c++0x 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
    LMPI    := -lmpi
endif





