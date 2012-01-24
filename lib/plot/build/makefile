# \file    makefile
# \ingroup plot
# \brief   Makefile for the plotting programs.
# \author  PA Jacobs
# \version 2.0
#
# \usage   make
#          make TARGET=for_gnu
#          make install
#          make clean
#          make tarfile
#----------------------------------------------------------------------

TARGET ?= for_gnu
INSTALL_DIR ?= $(HOME)/cfd_bin
include ../../util/source/systems.mk

SRC         := ../source
GEOM        := ../../geometry
GEOM_SRC    := $(GEOM)/source
LIBGEOM     := $(GEOM)/build/libgeom.a
USRC        := ../../util/source

#------------------------------------------------------------

ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    GIFMERGE_EXE := gifmerge.exe
else
    ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
        # CYGWIN environment on MS-Windows
        GIFMERGE_EXE := gifmerge.exe
    else
        # Unix/Linux
        GIFMERGE_EXE := gifmerge
    endif
endif
EXE_FILES := mb_cont.exe $(GIFMERGE_EXE) mb_view

SSP_OBJECTS := logfile.o \
	ssp.o \
	sspkern.o \
	sspmisc.o \
	ssp_colours.o \
	sspcont.o \
	ssp3d.o \
	gd.o

#------------------------------------------------------------

build_new: $(EXE_FILES)
	@echo "------------------------------------------"
	@echo "Plotting programs should be now up-to-date"
	@echo "------------------------------------------"

install: $(EXE_FILES) $(INSTALL_DIR)
	@echo Copy executable files to $(INSTALL_DIR)
	cp -f $(EXE_FILES) $(INSTALL_DIR)
	cp -f $(SRC)/mb_view.tcl $(INSTALL_DIR)
	chmod +x $(INSTALL_DIR)/mb_view.tcl

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)

clean:
	@echo Remove exe files
	- rm $(EXE_FILES)
	@echo Remove object files
	- rm $(SSP_OBJECTS) contour.o
	@echo Remove library file
	- rm ssp.a
	# geometry library
	cd $(GEOM)/build; make clean

tarfile:
	@echo Produce and compress a tar-file.
	- cd $(HOME); mv cfd_plot.tar.gz cfd_plot.tar.gz.old
	cd $(HOME); tar --exclude ".svn" --exclude "*~" -cvf cfd_plot.tar cfcfd/code/plot
	cd $(HOME); gzip cfd_plot.tar


#------------------------------------------------------------

# -----------------------
# Library of code modules
# -----------------------

ssp.a : $(SSP_OBJECTS)
	- rm ssp.a
	$(AR) rv ssp.a $(SSP_OBJECTS)
	- ranlib ssp.a

$(LIBGEOM) :
	cd $(GEOM)/build; make TARGET=$(TARGET) libgeom.a

# ----------------
# Executable files
# ----------------

mb_cont.exe : contour.o ssp.a $(LIBGEOM)
	$(LINK) $(LFLAG) contour.o ssp.a $(LIBGEOM) $(LLIB) \
		-o mb_cont.exe

$(GIFMERGE_EXE) : $(SRC)/gifmerge.c $(SRC)/gifmerge.h
	$(COMPILE) $(SRC)/gifmerge.c -o $(GIFMERGE_EXE)

mb_view :
	echo "wish $(INSTALL_DIR)/mb_view.tcl" > ./mb_view
	chmod +x ./mb_view

# Loadable module for Tcl/Tk -- Needs fixing sometime.
ssp_tk.so : sspd.o ssp.a $(LIBGEOM)
	$(LINK) $(LFLAG) sspd.so ssp.a $(LIBGEOM) $(LLIB) -o ssp_tk.so 

#-------------------------------------------------------------

# ------------
# Object files
# ------------

contour.o : $(SRC)/contour.c $(USRC)/compiler.h $(SRC)/ssp.h $(USRC)/logfile.h $(GEOM_SRC)/geom.h
	$(COMPILE) $(CFLAG) $(SRC)/contour.c -o contour.o

ssp.o     : $(SRC)/ssp.c $(SRC)/ssp.h $(GEOM_SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp.c -o ssp.o

sspmisc.o  : $(SRC)/sspmisc.c $(SRC)/ssp.h $(GEOM_SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/sspmisc.c -o sspmisc.o

sspkern.o  : $(SRC)/sspkern.c $(SRC)/ssp.h $(SRC)/ssp_colours.h $(GEOM_SRC)/geom.h \
	$(USRC)/compiler.h $(SRC)/gd.h
	$(COMPILE) $(CFLAG) $(SRC)/sspkern.c -o sspkern.o

ssp_colours.o  : $(SRC)/ssp_colours.c $(SRC)/ssp.h $(SRC)/ssp_colours.h $(GEOM_SRC)/geom.h \
	$(USRC)/compiler.h $(SRC)/gd.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp_colours.c -o ssp_colours.o

sspcont.o  : $(SRC)/sspcont.c $(SRC)/ssp.h $(GEOM_SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/sspcont.c -o sspcont.o

ssp3d.o    : $(SRC)/ssp3d.c $(SRC)/ssp.h $(GEOM_SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp3d.c -o ssp3d.o

ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
gd.o       : $(SRC)/gd.c $(SRC)/gd.h $(SRC)/mtables.c
	$(COMPILE) $(CFLAG) $(SRC)/gd.c -o gd.o
else
gd.o       : $(SRC)/gd.c $(SRC)/gd.h $(SRC)/mtables.c
	$(COMPILE) $(CFLAG) -I/usr/include/sys $(SRC)/gd.c -o gd.o
endif

logfile.o  : $(USRC)/logfile.c $(USRC)/logfile.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(USRC)/logfile.c -o logfile.o
#
#----------------------- end of makefile ----------------------------

