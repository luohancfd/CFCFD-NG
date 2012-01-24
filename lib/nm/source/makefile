# makefile
#
# This makefile compiles some of the testing routines within
# the numerical methods directory.
#
# R.J.Gollan
# 18-Aug-04

COMPILE=gcc
LINK=gcc
CFLAG=-fPIC -c -O2 -Wall -pedantic
LFLAG=-O -Wall -pedantic
LLIB=-lm

all: ode_test newt_raph_test vecmat_test

ode_test: ode_solvers_test.exe

newt_raph_test: newt_raph_test.exe

vecmat_test: vecmat_test.exe

ode_solvers_test.exe: Tode_solvers.o vecmat.o newt_raph.o
	$(LINK)  Tode_solvers.o vecmat.o newt_raph.o -o ode_solvers_test.exe $(LLIB)

Tode_solvers.o: ode_solvers.c ode_solvers.h vecmat.c vecmat.h
	$(COMPILE) -DTEST_ODE $(CFLAG) ode_solvers.c -o Tode_solvers.o

newt_raph_test.exe: Tnewt_raph.o vecmat.o
	$(LINK) $(LFLAG) Tnewt_raph.o vecmat.o -o newt_raph_test.exe $(LLIB)

Tnewt_raph.o: newt_raph.c newt_raph.h vecmat.c vecmat.h
	$(COMPILE) -DTEST_NEWT_RAPH $(CFLAG) newt_raph.c -o Tnewt_raph.o

newt_raph.o : newt_raph.c newt_raph.h vecmat.c vecmat.h
	$(COMPILE) $(CFLAG) newt_raph.c -o newt_raph.o

Tvecmat.o : vecmat.c vecmat.h
	$(COMPILE) -DTEST_MAIN $(CFLAG) vecmat.c -o Tvecmat.o

vecmat.o: vecmat.c vecmat.h
	$(COMPILE) $(CFLAG) vecmat.c -o vecmat.o

vecmat_test.exe: Tvecmat.o
	$(LINK) $(LFLAG) Tvecmat.o -o vecmat_test.exe $(LLIB)

clean:
	rm *.o *.exe

