#! /bin/sh
# post.sh

JOB=flat_runup

e3post.py --job=$JOB --add-mach --tindx=all --vtk-xml

#--tindx=$TINDX
# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=$JOB-y-wall.dat  \
    --slice-list="0,:,0,:;8,:,0,:;16,:,0,:;24,:,0,:;32,:,0,:;40,:,0,:;48,:,0,:;56,:,0,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-206mm-slice.dat  \
    --slice-list="56:63,-19,:,0"

e3post.py --job=$JOB --output-file=$JOB-204mm-slice.dat  \
    --slice-list="56:63,-20,:,0"

# Runs python script to compute viscous parameters.
echo "Computing viscous parameters .."
python compute_viscous_data_simple.py

echo "Postprocessing completed."
