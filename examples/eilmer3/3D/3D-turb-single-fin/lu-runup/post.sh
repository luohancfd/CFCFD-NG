#! /bin/sh
# post.sh

# USED TO GET Y WALL DATA FOR Y+ DETERMINATION THROUGH compute_viscous_data_simple.py
# CAN REMOVE ALL OTHER POST PROCESSING

JOB=flat_runup

e3post.py --job=$JOB --add-mach --tindx=all --vtk-xml

#--tindx=$TINDX
# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=$JOB-y-wall.dat  \
    --slice-list="0,:,0,:;4,:,0,:;8,:,0,:;12,:,0,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-178mm-slice.dat  \
    --slice-list="12,34,:,:;13,34,:,:;14,34,:,:;15,34,:,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-178mm-slice2.dat  \
    --slice-list="12,38,:,:;13,38,:,:;14,38,:,:;15,38,:,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-178mm-slice3.dat  \
    --slice-list="12,44,:,:;13,44,:,:;14,44,:,:;15,44,:,:" \
    --add-pitot-p --add-mach

e3post.py --job=$JOB --output-file=$JOB-178mm-slice4.dat  \
    --slice-list="12,50,:,:;13,50,:,:;14,50,:,:;15,50,:,:" \
    --add-pitot-p --add-mach

#--tindx=$TINDX
# Runs python script to compute viscous parameters.
echo "Computing viscous parameters .."
python compute_viscous_data_simple.py

echo "Postprocessing completed."
