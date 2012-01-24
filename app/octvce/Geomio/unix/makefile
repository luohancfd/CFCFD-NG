# \file    makefile
# \ingroup Geomio
# \brief   Makefile for geometry engine (point inclusion test for polygons/polyhedrons).

CC = gcc
CFLAGS = 
SRC_DIR = ../source/
OUT_DIR = ../tests/

#Main makes
#------------------------------------------------------

polyh_test : polyh_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o
	$(CC) $(CFLAGS) polyh_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)polyh_test.exe -lm

vtk_test: vtk_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o
	$(CC) $(CFLAGS) vtk_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)vtk_test.exe -lm

poly_test : poly_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o
	$(CC) $(CFLAGS) poly_test.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)poly_test.exe -lm

test_collinear : test_collinear.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o
	$(CC) $(CFLAGS) test_collinear.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)test_collinear.exe

test_sideloc : test_sideloc.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o
	$(CC) $(CFLAGS) test_sideloc.o gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)test_sideloc.exe -lm

clean :
	rm *.o
	rm *~*

#------------------------------------------------------

vtk_test.o : $(SRC_DIR)vtk_test.c 
	$(CC) -c $(CFLAGS) $(SRC_DIR)vtk_test.c -o vtk_test.o

polyh_test.o : $(SRC_DIR)polyh_test.c $(SRC_DIR)gio_io.c
	$(CC) -c $(CFLAGS) $(SRC_DIR)polyh_test.c -o polyh_test.o

poly_test.o : $(SRC_DIR)poly_test.c
	$(CC) -c $(CFLAGS) $(SRC_DIR)poly_test.c -o poly_test.o

test_collinear.o : $(SRC_DIR)test_collinear.c
	$(CC) -c $(CFLAGS) $(SRC_DIR)test_collinear.c -o test_collinear.o

test_sideloc.o : $(SRC_DIR)test_sideloc.c
	$(CC) -c $(CFLAGS) $(SRC_DIR)test_sideloc.c -o test_sideloc.o

gio_kernel.o : $(SRC_DIR)gio_kernel.c $(SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(SRC_DIR)gio_kernel.c -o gio_kernel.o

gio_adts.o : $(SRC_DIR)gio_adts.c $(SRC_DIR)gio_adts.h $(SRC_DIR)gio_kernel.h $(SRC_DIR)gio_lists.h
	$(CC) -c $(CFLAGS) $(SRC_DIR)gio_adts.c -o gio_adts.o

gio_io.o : $(SRC_DIR)gio_io.c $(SRC_DIR)gio_io.h $(SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(SRC_DIR)gio_io.c -o gio_io.o

gio_lists.o : $(SRC_DIR)gio_lists.c $(SRC_DIR)gio_lists.h $(SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(SRC_DIR)gio_lists.c -o gio_lists.o

#------------------------------------------------------
