#! /bin/sh
# post.sh

JOB=flat_runup

e3post.py --job=$JOB --add-mach --tindx=all --vtk-xml

# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=$JOB-y-wall.dat  \
    --slice-list="0,:,0,:;4,:,0,:;8,:,0,:;12,:,0,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-endghost1.dat  \
    --slice-list="12,-1,:,:;13,-1,:,:;14,-1,:,:;15,-1,:,:" 

e3post.py --job=$JOB --output-file=$JOB-endghost2.dat  \
    --slice-list="12,-2,:,:;13,-2,:,:;14,-2,:,:;15,-2,:,:" \

# Runs python script to compute viscous parameters.
echo "Computing viscous parameters .."
python compute_viscous_data_simple.py

echo "Postprocessing completed."
