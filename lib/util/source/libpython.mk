# \file    libpython.mk
# \ingroup util
# \brief   Build a Python library that can linked to the CFCFD codes.
# \author  PA Jacobs
#
# On MS-Windows, Python typically doesn't come with an import library
# that is ready to use with MinGW.
#
# The following instructions have been gleaned from the web at
# http://sebsauvage.net/python/mingw.html
#
# \usage   make -f libpython.mk
#          make -f libpython.mk clean
#----------------------------------------------------------------------

TARGET ?= for_gnu
include ../../util/source/systems.mk

$(LIBPYTHON) : lib$(PYTHON_NN).a
	cp -p lib$(PYTHON_NN).a $(PYTHON_DIR)/libs

lib$(PYTHON_NN).a : $(PYTHON_NN).dll $(PYTHON_NN).def
	dlltool --dllname $(PYTHON_NN).dll --def $(PYTHON_NN).def \
		--output-lib lib$(PYTHON_NN).a

$(PYTHON_NN).def : $(PYTHON_NN).dll
	# We may need to run the following command from the Win32 cmd shell
	# because the MinGW bash shell doesn't seem to pass the command-line
	# arguments to pexports.
	# Also, don't strip off the leading underscore as described on some
	# web pages.
	pexports $(PYTHON_NN).dll > $(PYTHON_NN).def

$(PYTHON_NN).dll : /c/WINDOWS/system32/$(PYTHON_NN).dll
	cp -p /c/WINDOWS/system32/$(PYTHON_NN).dll $(PYTHON_NN).dll
