## \file makefile
## \ingroup util
## \brief Makefile for the libutil module.
## \author RJG
## \version 04-May-2007
##
#----------------------------------------------------

TARGET ?= for_gnu
OPT ?= -O2
include ../source/systems.mk

SRC     := ../source
LIBUTIL := libutil.a

UTIL_OBJECTS := dbc_assert.o \
	config_parser.o \
	string_util.o

#----------------------------------------------------

build_new: $(LIBUTIL)
	@echo "Util library should be up to date."

clean:
	- rm $(LIBUTIL)
	- rm *.o

#----------------------------------------------------

$(LIBUTIL) : $(UTIL_OBJECTS)
ifneq ($(strip $(wildcard $(LIBUTIL))),)
	- rm $(LIBUTIL)
endif
	$(AR) rv $(LIBUTIL) $(UTIL_OBJECTS)
	- ranlib $(LIBUTIL)

#----------------------------------------------------

config_parser.o : $(SRC)/config_parser.cxx $(SRC)/config_parser.hh
	$(CXX) $(CXXFLAG) $(SRC)/config_parser.cxx

string_util.o : $(SRC)/string_util.cxx $(SRC)/string_util.hh
	$(CXX) $(CXXFLAG) $(SRC)/string_util.cxx

dbc_assert.o : $(SRC)/dbc_assert.cxx $(SRC)/dbc_assert.hh
	$(CXX) $(CXXFLAG) $(SRC)/dbc_assert.cxx



