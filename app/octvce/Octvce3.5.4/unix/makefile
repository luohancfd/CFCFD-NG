# \file    makefile
# \ingroup Octvce
# \brief   Makefile for Octree Virtual Cell Embedding blast wave simulator.

CC = gcc 
CFLAGS = -O2 
OV_SRC_DIR = ../source/
GM_SRC_DIR = ../../Geomio/source/
OUT_DIR = ../tests/

#------------------------------------------------------

octvce.exe : octvce.o ov_adts.o ov_adapt.o ov_connect.o ov_setgeom.o ov_icbc.o ov_io.o ov_lists.o \
	     ov_flow.o ov_merge.o time_keeper.o ov_recon.o ov_thermo.o ov_kernel.o \
	     gio_kernel.o gio_adts.o gio_io.o gio_lists.o		  
	$(CC) $(CFLAGS) octvce.o ov_adts.o ov_adapt.o ov_connect.o ov_setgeom.o \
	ov_icbc.o ov_io.o ov_lists.o ov_flow.o ov_merge.o ov_recon.o ov_thermo.o time_keeper.o ov_kernel.o \
	gio_kernel.o gio_adts.o gio_io.o gio_lists.o -o $(OUT_DIR)octvce.exe -lm

clean :
	rm -f *.o
	rm -f *.oo
	rm -f *~*

octvce.o : $(OV_SRC_DIR)octvce.c $(OV_SRC_DIR)ov_adts.h $(OV_SRC_DIR)ov_adapt.h $(OV_SRC_DIR)ov_io.h \
	   $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_flow.h $(OV_SRC_DIR)time_keeper.h \
	   $(OV_SRC_DIR)ov_lists.h		
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)octvce.c -o octvce.o

time_keeper.o : $(OV_SRC_DIR)time_keeper.c $(OV_SRC_DIR)time_keeper.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)time_keeper.c -o time_keeper.o

ov_flow.o : $(OV_SRC_DIR)ov_flow.c $(OV_SRC_DIR)ov_flow.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_adapt.h \
	    $(OV_SRC_DIR)ov_connect.h $(OV_SRC_DIR)ov_icbc.h $(OV_SRC_DIR)ov_io.h \
	    $(OV_SRC_DIR)ov_merge.h $(OV_SRC_DIR)ov_setgeom.h $(OV_SRC_DIR)ov_recon.h $(OV_SRC_DIR)ov_thermo.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_flow.c -o ov_flow.o

ov_adts.o : $(OV_SRC_DIR)ov_adts.c $(OV_SRC_DIR)ov_adts.h $(OV_SRC_DIR)ov_kernel.h 
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_adts.c -o ov_adts.o 

ov_adapt.o : $(OV_SRC_DIR)ov_adapt.c $(OV_SRC_DIR)ov_adapt.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_connect.h \
	     $(OV_SRC_DIR)ov_lists.h $(OV_SRC_DIR)ov_setgeom.h $(OV_SRC_DIR)ov_icbc.h $(OV_SRC_DIR)ov_merge.h \
	     $(OV_SRC_DIR)ov_recon.h $(OV_SRC_DIR)ov_thermo.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_adapt.c -o ov_adapt.o 

ov_connect.o : $(OV_SRC_DIR)ov_connect.c $(OV_SRC_DIR)ov_connect.h $(OV_SRC_DIR)ov_kernel.h \
	       $(OV_SRC_DIR)ov_lists.h $(OV_SRC_DIR)ov_setgeom.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_connect.c -o ov_connect.o 

ov_setgeom.o : $(OV_SRC_DIR)ov_setgeom.c $(OV_SRC_DIR)ov_setgeom.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_adts.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_setgeom.c -o ov_setgeom.o 

ov_icbc.o : $(OV_SRC_DIR)ov_icbc.c $(OV_SRC_DIR)ov_icbc.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_thermo.h  
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_icbc.c -o ov_icbc.o 

ov_io.o : $(OV_SRC_DIR)ov_io.c $(OV_SRC_DIR)ov_io.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_recon.h \
	  $(OV_SRC_DIR)ov_thermo.h $(OV_SRC_DIR)ov_adapt.h $(OV_SRC_DIR)ov_lists.h 
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_io.c -o ov_io.o 

ov_lists.o : $(OV_SRC_DIR)ov_lists.c $(OV_SRC_DIR)ov_lists.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_setgeom.h 
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_lists.c -o ov_lists.o 

ov_merge.o : $(OV_SRC_DIR)ov_merge.c $(OV_SRC_DIR)ov_merge.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_lists.h 
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_merge.c -o ov_merge.o

ov_recon.o : $(OV_SRC_DIR)ov_recon.c $(OV_SRC_DIR)ov_recon.h $(OV_SRC_DIR)ov_kernel.h $(OV_SRC_DIR)ov_thermo.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_recon.c -o ov_recon.o

ov_thermo.o : $(OV_SRC_DIR)ov_thermo.c $(OV_SRC_DIR)ov_thermo.h $(OV_SRC_DIR)ov_kernel.h
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_thermo.c -o ov_thermo.o

ov_kernel.o : $(OV_SRC_DIR)ov_kernel.c $(OV_SRC_DIR)ov_kernel.h \
	      $(GM_SRC_DIR)gio_kernel.h $(GM_SRC_DIR)gio_adts.h $(GM_SRC_DIR)gio_io.h $(GM_SRC_DIR)gio_lists.h   
	$(CC) -c $(CFLAGS) $(OV_SRC_DIR)ov_kernel.c -o ov_kernel.o

#Compile geometry engine ------------------------------

gio_kernel.o : $(GM_SRC_DIR)gio_kernel.c $(GM_SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(GM_SRC_DIR)gio_kernel.c -o gio_kernel.o

gio_adts.o : $(GM_SRC_DIR)gio_adts.c $(GM_SRC_DIR)gio_adts.h $(GM_SRC_DIR)gio_kernel.h $(GM_SRC_DIR)gio_lists.h
	$(CC) -c $(CFLAGS) $(GM_SRC_DIR)gio_adts.c -o gio_adts.o

gio_io.o : $(GM_SRC_DIR)gio_io.c $(GM_SRC_DIR)gio_io.h $(GM_SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(GM_SRC_DIR)gio_io.c -o gio_io.o

gio_lists.o : $(GM_SRC_DIR)gio_lists.c $(GM_SRC_DIR)gio_lists.h $(GM_SRC_DIR)gio_kernel.h
	$(CC) -c $(CFLAGS) $(GM_SRC_DIR)gio_lists.c -o gio_lists.o

#------------------------------------------------------
