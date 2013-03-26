## \file makefile
## \ingroup util
## \brief Makefile for the testing of util source
## \author RJG
## \version 12-Feb-06
##
## \usage   make
##          make TARGET=for_gnu
##          make install
##          make clean
#----------------------------------------------------------------------

TARGET ?= for_gnu
OPT ?= -O2
include ../source/systems.mk

SRC             := ../source
INSTALL_DIR     := $(HOME)/cfd_bin

#----------------------------------------------------------------------

EXE_FILES := config_parser_test.x dbc_assert_test.x

UTIL_OBJECTS :=	config_parser.o 

#------------------------------------------------------------

build_new: $(EXE_FILES)
	@echo "Util testing should be up to date."

clean:
	-rm $(EXE_FILES) 
	-rm *.o

#--------------- EXE Files ----------------------------------

config_parser_test.x : config_parser_test.o config_parser.o
	$(CXXLINK) -o config_parser_test.x config_parser_test.o config_parser.o

dbc_assert_test.x : 
	$(CXXLINK) -o dbc_assert_test.x $(SRC)/dbc_assert_test.cxx

#--------------- Object files -------------------------------

config_parser.o : $(SRC)/config_parser.cxx $(SRC)/config_parser.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/config_parser.cxx


config_parser_test.o : $(SRC)/config_parser_test.cxx $(SRC)/config_parser.hh
	$(CXXCOMPILE) $(CXXFLAG) $(SRC)/config_parser_test.cxx
