# makefile for td_browser filter extension.
#
# Updated 06-Jun-2007 to use the newer cfcfd build files.

TARGET ?= for_gnu
include ../../../lib/util/source/systems.mk


# For GNU/Linux
# COMPILE=gcc -c -Wall
# LINK=ld -shared

# For Solaris at DLR, Goettingen
# TCLTK=$(HOME)/local/ActiveTcl/include
# BLT=$(HOME)/local/blt/include
# COMPILE=cc -c -I$(TCLTK) -I$(BLT)
# LINK=ld -G

td_br_filter.so : td_br_filter_wrap.o td_br_filter.o
	$(LINK) $(LFLAG) -shared -o td_br_filter.so td_br_filter.o td_br_filter_wrap.o

td_br_filter.o : td_br_filter.c td_br_filter.h
	$(COMPILE) $(CFLAG) -I$(TCL_INCLUDE_DIR) td_br_filter.c

td_br_filter_wrap.o : td_br_filter_wrap.c
	$(COMPILE) $(CFLAG) -I$(TCL_INCLUDE_DIR) td_br_filter_wrap.c

# We want to update wrapper code.
td_br_filter_wrap.c : td_br_filter.i td_br_filter.h
	swig -tcl td_br_filter.i

tags:
	etags ./*.tcl ./*.c

clean:
	- rm td_br_filter.o td_br_filter_wrap.o
