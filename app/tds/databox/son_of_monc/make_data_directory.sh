#!/bin/bash
#
# make_data_directory.sh
#
# This shell script was written to make the necessary directory in 
# /data/tunnels for a new shot
#
# Carolyn Jacobs
# 12-Apr-07

# Step 1: User defines the shot number to be created
read -p "Shot number to be created (eg. x2s100):" -r SHOT

# Step 2: Make the directory
mkdir /data/X2/$SHOT

printf "Created directory for shot number %s \n" $SHOT


##################
#  NOTE TO SELF
#    - need to allow for directories already existing since it can't write over
##################