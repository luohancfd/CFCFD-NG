#! /bin/sh
# post.sh

# USED TO GET Y WALL DATA FOR Y+ DETERMINATION THROUGH compute_viscous_data_simple.py
# CAN REMOVE ALL OTHER POST PROCESSING

JOB=inject
TINDX=31

#--tindx=$TINDX
# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=$JOB-y-slice-UPSTREAM-$TINDX.dat  \
    --slice-list="0,:,0,0" --tindx=$TINDX

e3post.py --job=$JOB --output-file=$JOB-y-slice-DOWNSTREAM-$TINDX.dat  \
    --slice-list="48,:,0,0;56,:,0,0" --tindx=$TINDX

# Runs python script to compute viscous parameters.

echo "Postprocessing completed."
