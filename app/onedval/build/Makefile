# Makefile for onedval
#
# RJG, 12-Feb-2013

include ../../../lib/util/source/systems.mk
INSTALL_DIR ?= $(HOME)/e3bin
SOURCE := ../source
EILMER3 := ../../eilmer3
GAS := ../../../lib/gas
TARGET := for_gnu
OPT := -O2
PROGRAM := onedval.py
MODULES := cell.py \
	prop_avg.py

install : $(INSTALL_DIR) $(PROGRAM) 
	@echo Copying executable files to $(INSTALL_DIR)
	- cp -f $(PROGRAM) $(INSTALL_DIR)/onedval
	- chmod +x $(INSTALL_DIR)/onedval
	@echo Copying module files
	cp $(MODULES) $(INSTALL_DIR)
	@echo Building lib/gas and installing to $(INSTALL_DIR)
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	@echo Copying cfpylib files to $(INSTALL_DIR)
	cd ../../../lib/; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)
	@echo Building eilmer3 to get its support routines
	cd $(EILMER3)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	@echo --------------------------------
	@echo onedval should now be installed.
	@echo --------------------------------

$(INSTALL_DIR) :
	mkdir $(INSTALL_DIR)

$(PROGRAM) : $(SOURCE)/onedval.py $(MODULES)
	cp $(SOURCE)/onedval.py .

cell.py : $(SOURCE)/cell.py
	cp $(SOURCE)/cell.py .

prop_avg.py : $(SOURCE)/prop_avg.py
	cp $(SOURCE)/prop_avg.py .

clean :
	@echo Deleting object files
	- rm *.pyc *.py
	@echo Cleaning up other areas
	- cd $(EILMER3)/build; make clean
	- cd $(GAS)/build; make clean
