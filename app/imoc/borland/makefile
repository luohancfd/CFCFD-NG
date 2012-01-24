# makefile for Borland CBuilder 5.5
#
# This makefile assumes that SWIG has done its
# job of generating the wrapper-interface files.
# It also assumes that Tcl has been installed on
# the same drive in directory "\Program Files\Tcl"
#
# PJ, 07-May-2000

# Set the following macros to point to the right places.
SRC=..\source
# TCL=c:\progra~1\Tcl
TCL=c:\Tcl
TCLVER=84
BCB=d:\Borland\bcc55
# BCB=c:\progra~1\Borland\CBuilder5

IPATH=$(SRC);$(TCL)/include
COPTIONS=-c -3 -pc -tWD -tWM- -r- -X -Od -k -y -v -w-8004 -w-8057 -I$(IPATH)
LOPTIONS=-aa -Tpd -x -Gn -v -L$(BCB)\Lib

OBJECTS=imoc_wrap.obj \
	moc_kernel.obj \
	moc_gasdynamic.obj \
	moc_unitproc.obj \
	moc_wall.obj \
	moc_bezier.obj

all: imoc$(TCLVER).dll
        echo "Should have rebuilt imoc$(TCLVER).DLL."

imoc$(TCLVER).dll : $(OBJECTS) tcl$(TCLVER)_import.lib
        ilink32 $(LOPTIONS) c0d32.obj $(OBJECTS), imoc$(TCLVER)\
		,, tcl$(TCLVER)_import IMPORT32 CW32

# One of the extra tasks on Win32 is to generate an import library
# so that the object files can be linked against it.
tcl$(TCLVER)_import.lib : $(TCL)\bin\tcl$(TCLVER).dll
	impdef tcl.def0 $(TCL)\bin\tcl$(TCLVER).dll
	$(TCL)\bin\tclsh alias.tcl tcl.def tcl.def0
        implib tcl$(TCLVER)_import.lib tcl.def
	- del tcl.def
	- del tcl.def0

imoc_wrap.obj : $(SRC)/imoc_wrap.c $(SRC)/moc_kernel.h $(SRC)/moc_wall.h \
		$(SRC)/moc_gasdynamic.h $(SRC)/moc_unitproc.h
        bcc32 $(COPTIONS) -D__WIN32__ $(SRC)/imoc_wrap.c

moc_wall.obj : $(SRC)/moc_wall.c $(SRC)/moc_wall.h $(SRC)/moc_bezier.h
        bcc32 $(COPTIONS) $(SRC)/moc_wall.c

moc_bezier.obj : $(SRC)/moc_bezier.c $(SRC)/moc_bezier.h
        bcc32 $(COPTIONS) $(SRC)/moc_bezier.c

moc_gasdynamic.obj : $(SRC)/moc_gasdynamic.c $(SRC)/moc_gasdynamic.h
        bcc32 $(COPTIONS) $(SRC)/moc_gasdynamic.c

moc_unitproc.obj : $(SRC)/moc_unitproc.c $(SRC)/moc_unitproc.h $(SRC)/moc_gasdynamic.h $(SRC)/moc_kernel.h
        bcc32 $(COPTIONS) $(SRC)/moc_unitproc.c

moc_kernel.obj : $(SRC)/moc_kernel.c $(SRC)/moc_kernel.h $(SRC)/moc_gasdynamic.h
        bcc32 $(COPTIONS) $(SRC)/moc_kernel.c

# pkgIndex.tcl : moc_kernel.tcl moc_gui.tcl moc_plot.tcl \
#         moc_menu.tcl moc_placard.tcl moc_nodelist.tcl \
#         moc_scales.tcl moc_unitproc.tcl make_index.tcl
#         tclsh make_index.tcl

clean :: 
        - del *.obj
	- del *.tds

allclean ::
        - del *.dll
        - del *.obj
        - del *.lib
	- del *.tds
	- del imoc??.*




