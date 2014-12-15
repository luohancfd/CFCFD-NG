#! /bin/sh
# post.sh

JOB=tc2update

e3post.py --job=$JOB --add-mach --tindx=all --vtk-xml

# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=$JOB-y-wall.dat  \
    --slice-list="0,:,0,0;16,:,0,0;32,:,0,0;48,:,0,0" \

# Runs python script to compute viscous parameters.
echo "Computing viscous parameters .."
python compute_viscous_data_simple.py

echo "Postprocessing completed."
