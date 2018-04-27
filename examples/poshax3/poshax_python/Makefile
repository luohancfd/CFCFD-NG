# makefile for poshax python, just copies the poshax_python file, basically.
# this is based on the pitot one, which was based on the nenzfr one
# (so props to whoever did that)
#
# Chris James (c.james4@uq.edu.au) - 21-Apr-2014

INSTALL_DIR ?= $(HOME)/e3bin
 
install : $(INSTALL_DIR)
	@echo Copying poshax_python.py to the install directory $(INSTALL_DIR)
	- cp -f poshax_python.py $(INSTALL_DIR)
	
$(INSTALL_DIR):
	mkdir $(INSTALL_DIR)
