# \file Makefile 
# \brief Builds the executable file for the databox device-driver.
# \author ZS and PJ 
# \version October 2004
#          Generalise Zane's specific kernel 2.4 build commands.
#          Will build on kernel 2.6 (and later, hopefully).
# \version December 2004
#          Assumes Emstar is available.
#          Commands for setuid are included.

KERNEL_RELEASE := $(shell uname -r)
KERNEL_INCLUDE_DIR := /lib/modules/$(KERNEL_RELEASE)/build/include
EMSTAR_HOME := $(HOME)/emstar-09-oct-2004

dbox_driver : dbox_driver.c
	gcc -O2 -Wall -isystem $(KERNEL_INCLUDE_DIR) \
		-I$(EMSTAR_HOME)/fusd/include \
		-o dbox_driver dbox_driver.c \
		$(EMSTAR_HOME)/obj.i686-linux/lib/libfusd.a

permission : dbox_driver
	sudo chown root dbox_driver
	sudo chmod 06775 dbox_driver

clean :
	- rm dbox_driver
