#
# iniparser Makefile
#
SYSTEM := $(shell uname -s)

CC      = gcc
CFLAGS  = -O3 -fPIC
ifeq ($(findstring MINGW32, $(SYSTEM)), MINGW32)
    # MINGW32 environment on MS-Windows
    CFLAGS  = -O3
endif
ifeq ($(findstring CYGWIN, $(SYSTEM)), CYGWIN)
    # CYGWIN environment on MS-Windows
    CFLAGS  = -O3
endif
AR	    = ar
ARFLAGS = rcv
RM      = rm -f


# Implicit rules

SUFFIXES = .o .c .h .a .so .sl

COMPILE.c=$(CC) $(CFLAGS) -c
.c.o:
	@(echo "compiling $< ...")
	@($(COMPILE.c) -o $@ $<)


SRCS = src/iniparser.c \
	   src/dictionary.c \
	   src/strlib.c

OBJS = $(SRCS:.c=.o)


default:	libiniparser.a

libiniparser.a:	$(OBJS)
	@($(AR) $(ARFLAGS) libiniparser.a $(OBJS))
	- ranlib libiniparser.a

clean:
	$(RM) $(OBJS) libiniparser.a

veryclean:
	$(RM) $(OBJS) libiniparser.a
	rm -rf ./html ; mkdir html
	cd test ; $(MAKE) veryclean

docs:
	@(cd doc ; $(MAKE))

check:
	@(cd test ; $(MAKE))
