# makefile for nenzfr and estcj
# These Python programs depend upon Eilmer3 and the NASA CEA2 code.
#
# PJ, 30-Sep-2011

include ../../lib/util/source/systems.mk
INSTALL_DIR ?= $(HOME)/e3bin
EILMER3 := ../eilmer3
GAS     := ../../lib/gas
TARGET := for_gnu
OPT := -O2
PROGRAMS := estcj.py nenzfr.py nenzfr_compute_viscous_data.py nenzfr_quasi_transient.py nenzfr_sensitivity.py nenzfr_perturbed.py nenzfr_RSA.py

install : $(INSTALL_DIR) $(PROGRAMS)
	@echo Copying executable files to $(INSTALL_DIR)
	cd $(EILMER3)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd $(GAS)/build; make TARGET=$(TARGET) OPT=$(OPT) INSTALL_DIR=$(INSTALL_DIR) install
	cd ../../lib/; rsync -av --exclude=".svn" --exclude="*~" cfpylib $(INSTALL_DIR)
	rsync -av --exclude=".svn" --exclude="*~" nenzfr_data_files $(INSTALL_DIR)
	- cp -f $(PROGRAMS) $(INSTALL_DIR)
	- chmod +x $(INSTALL_DIR)/nenzfr.py
	- chmod +x $(INSTALL_DIR)/estcj.py
	@echo -----------------------------------------
	@echo nenzfr and estcj should be now installed.
	@echo -----------------------------------------

clean :
	@echo Deleting object files
	- rm *.pyc
	@echo Cleaning up other areas
	- cd $(EILMER3)/build; make clean
	- cd $(GAS)/build; make clean

$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)
