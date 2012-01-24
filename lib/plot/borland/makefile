# makefile
#
#    ************************************************
#    * Top-level Makefile for the plotting programs *
#    * using the Borland C/C++ compiler in Win32.   *
#    ************************************************
#
# P. Jacobs
# Department of Mechanical Engineering
# The University of Queensland
#
#------------------------------------------------------------

BCB=c:\progra~1\Borland\CBuilder5
COMPILE=bcc32
CFLAG=-c -DBCC32 -v -y -w-8004 -w-8057 -q
LINK=bcc32
LFLAG=-v -q

SRC=..\source
USRC=..\..\util\source
NMSRC=..\..\nm\source
INSTALL_DIR=..\..\..\cfd_bin

#------------------------------------------------------------

EXE_FILES = mb_cont.exe \
	gifmerge.exe \
	mb_view.tcl \
	mb_view.bat

SSP_OBJECTS = logfile.obj \
	geom.obj \
	ssp.obj \
	sspkern.obj \
	sspmisc.obj \
	ssp_colours.obj \
	sspcont.obj \
	ssp3d.obj \
	gd.obj

#------------------------------------------------------------

build_new: $(EXE_FILES)
	@echo "------------------------------------------"
	@echo "Plotting programs should be now up-to-date"
	@echo "------------------------------------------"

install: mb_cont.exe gifmerge.exe mb_view.bat mb_view.tcl
	@echo Copy executable files to $(INSTALL_DIR)
	- mkdir $(INSTALL_DIR)
	copy mb_cont.exe $(INSTALL_DIR)
	copy gifmerge.exe $(INSTALL_DIR)
	copy mb_view.bat $(INSTALL_DIR)
	copy mb_view.tcl $(INSTALL_DIR)

clean:
	@echo Remove exe files
	- del mb_cont.exe
	- del mb_view.bat
	- del gifmerge.exe
	@echo Remove object files
	- del contour.obj
	- del gifmerge.obj
	- del logfile.obj
	- del geom.obj
	- del ssp.obj
        - del sspkern.obj
	- del sspmisc.obj
	- del ssp_colours.obj
	- del sspcont.obj
	- del ssp3d.obj
	- del gd.obj
	@echo Remove library file
	- del ssp.lib
	@echo Remove Tcl scripts
	- del mb_view.tcl
	@echo Remove other compiler generated files
	- del *.tds
	- del ssp.bak

#------------------------------------------------------------

# -----------------------
# Library of code modules
# -----------------------

ssp.lib : $(SSP_OBJECTS)
	- del ssp.lib
	tlib ssp.lib +logfile.obj
	tlib ssp.lib +geom.obj
	tlib ssp.lib +ssp.obj
        tlib ssp.lib +sspkern.obj
	tlib ssp.lib +sspmisc.obj
	tlib ssp.lib +ssp_colours.obj
	tlib ssp.lib +sspcont.obj
	tlib ssp.lib +ssp3d.obj
	tlib ssp.lib +gd.obj

# ----------------
# Executable files
# ----------------

mb_cont.exe : contour.obj ssp.lib
	$(LINK) $(LFLAG) -emb_cont.exe contour.obj ssp.lib

gifmerge.exe : $(SRC)/gifmerge.c $(SRC)/gifmerge.h
	$(COMPILE) -egifmerge.exe $(SRC)/gifmerge.c

mb_view.bat : $(SRC)/mb_view.bat
	copy $(SRC)\mb_view.bat mb_view.bat

mb_view.tcl : $(SRC)/mb_view.tcl
	copy $(SRC)\mb_view.tcl mb_view.tcl

#-------------------------------------------------------------

# ------------
# Object files
# ------------

contour.obj : $(SRC)/contour.c $(USRC)/compiler.h $(SRC)/ssp.h \
		$(USRC)/logfile.h $(SRC)/geom.h
	$(COMPILE) $(CFLAG) $(SRC)/contour.c

geom.obj : $(SRC)/geom.c $(USRC)/compiler.h $(SRC)/geom.h
	$(COMPILE) $(CFLAG) $(SRC)/geom.c

ssp.obj     : $(SRC)/ssp.c $(SRC)/ssp.h $(SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp.c

sspmisc.obj  : $(SRC)/sspmisc.c $(SRC)/ssp.h $(SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/sspmisc.c

sspkern.obj  : $(SRC)/sspkern.c $(SRC)/ssp.h $(SRC)/ssp_colours.h $(SRC)/geom.h \
		$(USRC)/compiler.h $(SRC)/gd.h
	$(COMPILE) $(CFLAG) $(SRC)/sspkern.c

ssp_colours.obj  : $(SRC)/ssp_colours.c $(SRC)/ssp.h $(SRC)/ssp_colours.h \
		$(SRC)/geom.h $(USRC)/compiler.h $(SRC)/gd.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp_colours.c

sspcont.obj  : $(SRC)/sspcont.c $(SRC)/ssp.h $(SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/sspcont.c

ssp3d.obj    : $(SRC)/ssp3d.c $(SRC)/ssp.h $(SRC)/geom.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(SRC)/ssp3d.c

gd.obj       : $(SRC)/gd.c $(SRC)/gd.h $(SRC)/mtables.c
	$(COMPILE) $(CFLAG) $(SRC)/gd.c

logfile.obj  : $(USRC)/logfile.c $(USRC)/logfile.h $(USRC)/compiler.h
	$(COMPILE) $(CFLAG) $(USRC)/logfile.c


#
#=============== end of makefile ===================

