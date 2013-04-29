# makefile for pitot
# lines to install the gas library and cfpylib taken from the nenzfr / estcj make file
#
# Chris James (c.james4@uq.edu.au) - 22-Apr-2013

TARGET ?= for_gnu
OPT ?= -O2

INSTALL_DIR ?= $(HOME)/e3bin
GAS     := ../../lib/gas

install : $(INSTALL_DIR) pitot.py
	@echo Copying cfpylib files to $(INSTALL_DIR)
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd ../../lib/; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)
	@echo Copying pitot to $(INSTALL_DIR)
	mkdir $(INSTALL_DIR)/pitot
	- cp -f pitot.py $(INSTALL_DIR)/pitot/pitot.py
	- chmod +x $(INSTALL_DIR)/pitot/pitot.py
	@echo -----------------------------------------
	@echo pitot should be now installed.
	@echo -----------------------------------------
	
$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)
