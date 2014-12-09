# \brief   Makefile for the Method-of-Characteristics  program "IMOC".
# \author  PA Jacobs
# \version 2.0
#
# \usage   make
#          make TARGET=for_gnu
#          make install
#          make clean
#          make allclean
#          make tarfile
#          make zipfile
#----------------------------------------------------------------------

TARGET ?= for_gnu

SRC=../source
DOC=../doc
UTIL=../../../lib/util/source
INSTALL_DIR = $(HOME)/imoc_bin

include $(UTIL)/systems.mk

# Assume Unix/GNU-Linux
LOADABLE_MODULE := imoc.so
COMPILE := gcc -c -O -fPIC -DUSE_TCL_STUBS -I$(TCL_INCLUDE_DIR)
LINK := gcc -shared
LLIBS :=
# For Solaris at DLR, Goettingen
# COMPILE=cc -c -I/usr/local/tktcl-8.3.2/include/
# LINK=ld -G

ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    LINK += -L$(TCL_DIR)/lib
    LLIBS += -lm -ltclstub84
    LOADABLE_MODULE := imoc84.dll
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    LINK += -L$(TCL_DIR)/lib
    LLIBS += -lm -ltclstub84
    LOADABLE_MODULE := imoc84.dll
endif
ifeq ($(findstring Darwin, $(SYSTEM)), Darwin)
    # Darwin environment on a Mac OS X
    LINK := gcc -dynamiclib -DUSE_TCL_STUBS -framework Tcl
    LOADABLE_MODULE := imoc.dylib
endif

#----------------------------------------------------------------------

all: $(LOADABLE_MODULE) $(SRC)/pkgIndex.tcl documentation $(SRC)/moc_startup.tcl
	@echo "Tcl headers were found in $(TCL_INCLUDE_DIR)"
	@echo "Should have built shared-object and package-Index file."

$(LOADABLE_MODULE) : imoc_wrap.o moc_kernel.o moc_gasdynamic.o \
                moc_unitproc.o moc_wall.o moc_bezier.o
	$(LINK) -o $(LOADABLE_MODULE) \
		imoc_wrap.o moc_kernel.o moc_gasdynamic.o moc_wall.o \
		moc_bezier.o moc_unitproc.o $(LLIBS)

imoc_wrap.o : $(SRC)/imoc_wrap.c $(SRC)/moc_kernel.h $(SRC)/moc_wall.h \
		$(SRC)/moc_gasdynamic.h $(SRC)/moc_unitproc.h
	$(COMPILE) $(SRC)/imoc_wrap.c

$(SRC)/imoc_wrap.c : $(SRC)/imoc.i
	swig -tcl -pkgversion 0.1 $(SRC)/imoc.i

moc_wall.o : $(SRC)/moc_wall.c $(SRC)/moc_wall.h $(SRC)/moc_bezier.h
	$(COMPILE) $(SRC)/moc_wall.c

moc_bezier.o : $(SRC)/moc_bezier.c $(SRC)/moc_bezier.h
	$(COMPILE) $(SRC)/moc_bezier.c

moc_gasdynamic.o : $(SRC)/moc_gasdynamic.c $(SRC)/moc_gasdynamic.h
	$(COMPILE) $(SRC)/moc_gasdynamic.c

moc_unitproc.o : $(SRC)/moc_unitproc.c $(SRC)/moc_unitproc.h \
		$(SRC)/moc_gasdynamic.h $(SRC)/moc_kernel.h
	$(COMPILE) $(SRC)/moc_unitproc.c

moc_kernel.o : $(SRC)/moc_kernel.c $(SRC)/moc_kernel.h $(SRC)/moc_gasdynamic.h
	$(COMPILE) $(SRC)/moc_kernel.c

#----------------------------------------------------------------------
$(SRC)/pkgIndex.tcl : $(LOADABLE_MODULE) \
		$(SRC)/moc_kernel.tcl $(SRC)/moc_gui.tcl $(SRC)/moc_plot.tcl \
		$(SRC)/moc_menu.tcl $(SRC)/moc_placard.tcl $(SRC)/moc_nodelist.tcl \
		$(SRC)/moc_scales.tcl $(SRC)/moc_unitproc.tcl $(SRC)/make_index.tcl \
		$(SRC)/moc_syn_cmds.tcl
	@echo Building the Tcl autoloading package file
	cp $(LOADABLE_MODULE) $(SRC)
	cd $(SRC); tclsh make_index.tcl
	rm $(SRC)/$(LOADABLE_MODULE)

#----------------------------------------------------------------------
$(SRC)/moc_syn_cmds.tcl : $(SRC)/synonyms.table
	@echo Generating the synonymous command procedures
	cd $(SRC); awk -f make_syn_cmds.awk synonyms.table > moc_syn_cmds.tcl

#----------------------------------------------------------------------
documentation : $(UTIL)/doctcl.awk $(UTIL)/docgen.awk \
		$(DOC)/moc_kernel_c.html $(DOC)/moc_gasdynamic_c.html \
		$(DOC)/moc_wall_c.html $(DOC)/moc_bezier_c.html \
		$(DOC)/moc_unitproc_c.html \
		$(DOC)/moc_gui_tcl.html $(DOC)/moc_kernel_tcl.html \
		$(DOC)/moc_menu_tcl.html $(DOC)/moc_nodelist_tcl.html \
		$(DOC)/moc_placard_tcl.html $(DOC)/moc_scales_tcl.html \
		$(DOC)/moc_unitproc_tcl.html $(DOC)/moc_plot_tcl.html \
		$(DOC)/moc_syn_cmds_tcl.html
	@echo Should have built all HTML files

$(DOC)/moc_kernel_c.html : $(SRC)/moc_kernel.c 
	awk -f $(UTIL)/docgen.awk $(SRC)/moc_kernel.c > $(DOC)/moc_kernel_c.html

$(DOC)/moc_gasdynamic_c.html : $(SRC)/moc_gasdynamic.c 
	awk -f $(UTIL)/docgen.awk $(SRC)/moc_gasdynamic.c > $(DOC)/moc_gasdynamic_c.html

$(DOC)/moc_wall_c.html : $(SRC)/moc_wall.c 
	awk -f $(UTIL)/docgen.awk $(SRC)/moc_wall.c > $(DOC)/moc_wall_c.html

$(DOC)/moc_bezier_c.html : $(SRC)/moc_bezier.c 
	awk -f $(UTIL)/docgen.awk $(SRC)/moc_bezier.c > $(DOC)/moc_bezier_c.html

$(DOC)/moc_unitproc_c.html : $(SRC)/moc_unitproc.c 
	awk -f $(UTIL)/docgen.awk $(SRC)/moc_unitproc.c > $(DOC)/moc_unitproc_c.html

$(DOC)/moc_gui_tcl.html : $(SRC)/moc_gui.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_gui.tcl > $(DOC)/moc_gui_tcl.html

$(DOC)/moc_kernel_tcl.html : $(SRC)/moc_kernel.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_kernel.tcl > $(DOC)/moc_kernel_tcl.html

$(DOC)/moc_menu_tcl.html : $(SRC)/moc_menu.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_menu.tcl > $(DOC)/moc_menu_tcl.html

$(DOC)/moc_nodelist_tcl.html : $(SRC)/moc_nodelist.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_nodelist.tcl > $(DOC)/moc_nodelist_tcl.html

$(DOC)/moc_placard_tcl.html : $(SRC)/moc_placard.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_placard.tcl > $(DOC)/moc_placard_tcl.html

$(DOC)/moc_plot_tcl.html : $(SRC)/moc_plot.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_plot.tcl > $(DOC)/moc_plot_tcl.html

$(DOC)/moc_scales_tcl.html : $(SRC)/moc_scales.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_scales.tcl > $(DOC)/moc_scales_tcl.html

$(DOC)/moc_unitproc_tcl.html : $(SRC)/moc_unitproc.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_unitproc.tcl > $(DOC)/moc_unitproc_tcl.html

$(DOC)/moc_syn_cmds_tcl.html : $(SRC)/moc_syn_cmds.tcl 
	awk -f $(UTIL)/doctcl.awk $(SRC)/moc_syn_cmds.tcl > $(DOC)/moc_syn_cmds_tcl.html

#----------------------------------------------------------------------

install:
	- mkdir $(INSTALL_DIR); mkdir $(INSTALL_DIR)/source
	- mkdir $(INSTALL_DIR)/doc; mkdir $(INSTALL_DIR)/demo
	- mkdir $(INSTALL_DIR)/lib
	cp -f $(LOADABLE_MODULE) $(INSTALL_DIR)/source
	cp -f $(SRC)/*.tcl $(INSTALL_DIR)/source
	cp -f ../demo/*.* $(INSTALL_DIR)/demo
	cp -f ../lib/*.* $(INSTALL_DIR)/lib
	cp -f $(DOC)/*.html $(INSTALL_DIR)/doc

clean :
	@echo Remove the shared-object and object files
	- rm $(LOADABLE_MODULE) *.o

allclean :
	@echo Remove all files from this directory that were machine generated
	- rm $(LOADABLE_MODULE) *.o
	- rm $(SRC)/pkgIndex.tcl $(SRC)/imoc_wrap.c $(SRC)/imoc_wrap.doc

tarfile :
	@echo Build a compressed tar-file in the HOME directory
	- cd $(HOME); mv imoc.tar.gz imoc.tar.gz.old
	cd $(HOME); tar --exclude ".svn" --exclude "*~" -cvf imoc.tar cfcfd/code/imoc
	cd $(HOME); gzip imoc.tar

zipfile:
	@echo Build a ZIP file in the home directory
	-cd $(HOME); mv imoc.zip imoc.zip.old
	cd $(HOME); zip -r imoc cfcfd/code/imoc -x \*.svn\* \*\~
