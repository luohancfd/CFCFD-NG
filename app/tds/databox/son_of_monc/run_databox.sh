#!/bin/bash
#
# run_databox.sh
#
# This shell script was written to change the permissions on the serial and 
# usb ports and then to open the databox program as xlabuser
#
# Carolyn Jacobs
# 16-Feb-07

# Step 1: Change the permissions on the ports (as sudo)
sudo chmod 777 /dev/ttyS*
sudo chmod 777 /dev/ttyUSB0

# Step 2: Run the databox program as xlabuser
cd /home/xlabuser/databox/son_of_monc
./dbox_view_usb.py