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
    # Debian default tcl
    TCL_INCLUDE_DIR := $(dir $(wildcard /usr/include/tcl/tcl.h))
    ifeq ($(strip $(TCL_INCLUDE_DIR)),)
        # If we were unsuccessful, there may be a specific version directory.
        TCL_INCLUDE_DIR := $(dir $(firstword $(wildcard /usr/include/tcl8*/tcl.h)))
    endif
    ifeq ($(strip $(TCL_INCLUDE_DIR)),)
        # Unix/Linux with ActiveState Tcl installed.
        TCL_INCLUDE_DIR := $(dir $(wildcard /usr/local/ActiveTcl/include/tcl.h))
    endif
    ifeq ($(strip $(TCL_INCLUDE_DIR)),)
        # Still haven't found the tcl header.
        # Hopefully, there is a Tcl distribution in the last place we look...
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
GNU_SUFFIX ?= 

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
    CXXCOMPILE := clang
    CXXLINK := clang 
    # Unix/Linux is default
    CFLAG   := -c $(OPT)
    LFLAG   := $(OPT)
    CXXFLAG := -c $(OPT) -std=c++11
    LLIB := -lstdc++ -lm
endif

ifeq ($(TARGET), for_gnu)
    # UNIX/Linux workstation with the GNU C compiler
    COMPILE := gcc$(GNU_SUFFIX)
    LINK    := gcc$(GNU_SUFFIX)
    CXXCOMPILE := g++$(GNU_SUFFIX)
    CXXLINK := g++$(GNU_SUFFIX)
    # Unix/Linux is default
    CFLAG   := -c $(OPT) -fPIC -W -Wall -pedantic $(MARCH_FLAG)
    LFLAG   := $(OPT) -fPIC $(MARCH_FLAG)
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
        F90 := gfortran$(GNU_SUFFIX)
        # Compile without vector optimization for now.
        F90FLAG := -m64 -c -O0 -fPIC
        F90LFLAG := -m64 -lstdc++ -gnofor_main
        FLINK := -lgfortran
    endif
endif

ifeq ($(TARGET), for_gnu_debug)
    # UNIX/Linux workstation with the GNU C compiler and debug+profiling
    # 2014-05-20 Let's make use of some gcc 4.8.x features, including the address sanitizer.
    # Note, however, that only the statically linked executables are useful from this build.
    # You should build without installing and then manually copy just the executables
    # into $(E3BIN).
    COMPILE := gcc$(GNU_SUFFIX) 
    LINK    := gcc$(GNU_SUFFIX)
    CXXCOMPILE := g++$(GNU_SUFFIX)
    CXXLINK := g++$(GNU_SUFFIX)
    # Unix/Linux is default
    CFLAG   := -c -fPIC -W -Wall -pedantic -ggdb $(MARCH_FLAG)
    LFLAG   := -fPIC -pedantic -ggdb $(MARCH_FLAG)
    CXXFLAG := -c -fPIC -std=c++11 -Wall -pedantic -ggdb $(MARCH_FLAG)
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        CFLAG   := -c -W -Wall -pedantic -g $(MARCH_FLAG)
        LFLAG   := -g -Wl,stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c -std=c++11 -Wall -pedantic -g $(MARCH_FLAG)
    endif
    ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
        # MINGW32 environment on MS-Windows
        CFLAG   := -c -W -Wall -pedantic -g $(MARCH_FLAG)
        LFLAG   := -g -Wl,stack=0x8000000 $(MARCH_FLAG)
        CXXFLAG := -c -std=c++11 -Wall -pedantic -g $(MARCH_FLAG)
    endif
    CFLAG   += -ftrapping-math -fsignaling-nans -Og -fsanitize=address -fno-omit-frame-pointer -DDEBUG
    CXXFLAG += -ftrapping-math -fsignaling-nans -Og -fsanitize=address -fno-omit-frame-pointer -DDEBUG
    LFLAG   += -ftrapping-math -fsignaling-nans -Og -fsanitize=address -DDEBUG
    LLIB    := -lm
endif

ifeq ($(TARGET), for_gprof)
    # UNIX/Linux/Cygwin workstation with the GNU C compiler and with profiling.
    COMPILE := gcc 
    LINK    := gcc
    CXXCOMPILE := g++
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
    CXXCOMPILE := g++
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
    CXXCOMPILE := g++
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
    CXXCOMPILE := g++
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
    CXXCOMPILE := pgCC
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
    CXXCOMPILE := pgCC
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
    CXXCOMPILE := pgCC
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
    CXXCOMPILE := pgCC
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
    CXXCOMPILE := icpc
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
    CXXCOMPILE := icpc
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
    CXXCOMPILE := icpc
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
    CXXCOMPILE := mpiicpc
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
    CXXCOMPILE := mpiicpc
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
    CXXCOMPILE := mpiicpc
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
    CXXCOMPILE := icpc
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

ifeq ($(TARGET), for_openmpi)
    # OpenMPI on Linux.
    COMPILE := mpicc
    LINK    := mpicc
    CXXCOMPILE := mpicxx
    CXXLINK := mpicxx
    CFLAG   := -c $(OPT) -fPIC -Wall -pedantic 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -pedantic 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
endif

ifeq ($(TARGET), for_openmpi_debug)
    # OpenMPI on Linux.
    # 2014-05-20 Let's make use of some gcc 4.8.x features, including the address sanitizer.
    # Note, however, that only the statically linked executables are useful from this build.
    # You should build without installing and then manually copy just the executables
    # into $(E3BIN).
    COMPILE := mpicc
    LINK    := mpicc
    CXXCOMPILE := mpicxx
    CXXLINK := mpicxx
    CFLAG   := -c -fPIC -ggdb -Wall -pedantic -Og -fsanitize=address -fno-omit-frame-pointer 
    CXXFLAG := -c -std=c++11 -fPIC -ggdb -Wall -pedantic -Og -fsanitize=address -fno-omit-frame-pointer  
    LFLAG   := -fPIC -ggdb -Og -fsanitize=address
    LLIB    := -lm
endif

ifeq ($(TARGET), for_openmpi_intel)
    # OpenMPI on Linux.
    COMPILE := mpicc
    LINK    := mpicc
    CXXCOMPILE := mpicxx
    CXXLINK := mpicxx
    CFLAG   := -c $(OPT) -fPIC -Wall -pedantic 
    CXXFLAG := -c $(OPT) -std=c++11 -fPIC -Wall -pedantic 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm -limf -lirc
endif

ifeq ($(TARGET), for_macports_openmpi)
    # OpenMPI on Mac OS-X with MacPorts
    # PJ edits for Ingo's Mac 06-June-2012.
    # Note (1) This works on LinuxMint13 as well as the usual for_openmpi.
    # Note (2) The case sensitivity seems to work since Dan's note below,
    #          at least within the xcode unix-like environment that Ingo uses.
    COMPILE := mpicc
    LINK    := mpicc
    CXXCOMPILE := mpiCC
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
    CXXCOMPILE := mpic++
    CXXLINK := mpic++
    CFLAG   := -c $(OPT) -fPIC -Wall -pedantic 
    CXXFLAG := -c $(OPT) -std=c++0x -fPIC -Wall -pedantic 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
endif

ifeq ($(TARGET), for_mpich)
    # mpich on Linux.
    # Setting the Lam compiler-specific environment variables
    # will allow setting of options specific to each compiler.
    COMPILE := mpicc
    LINK    := mpicc
    CXXCOMPILE := mpiCC
    CXXLINK := mpiCC
    CFLAG   := -c $(OPT) 
    CXXFLAG := -c $(OPT) -std=c++0x 
    LFLAG   :=  $(OPT) -fPIC 
    LLIB    := -lm
    LMPI    := -lmpi
endif





