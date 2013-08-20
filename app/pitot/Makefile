# makefile for pitot
# lines to install the gas library and cfpylib taken from the nenzfr / estcj make file
#
# Chris James (c.james4@uq.edu.au) - 07-May-2013

TARGET ?= for_gnu
OPT ?= -O2

INSTALL_DIR ?= $(HOME)/e3bin
GAS     := ../../lib/gas
PROGRAMS := pitot.py pitot_input_utils.py pitot_flow_functions.py pitot_output_utils.py pitot_area_ratio_check.py

install : $(INSTALL_DIR) $(PROGRAMS)
	@echo Copying cfpylib files to $(INSTALL_DIR)
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd ../../lib/; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)
	@echo Copying pitot files to $(INSTALL_DIR)
	- cp -f $(PROGRAMS) $(INSTALL_DIR)
	- chmod +x $(INSTALL_DIR)/pitot.py
	@echo -----------------------------------------
	@echo pitot should be now installed.
	@echo -----------------------------------------
	
$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)
