#! /bin/sh
# post_march.sh

JOB=turb_flat_plate_march
TINDX=0001

# Extracts Paraview or Visit files for flowfield visualisation
e3post.py --job=$JOB --vtk-xml --tindx=$TINDX --add-pitot-p --add-mach

# Extracts slices at x=0.368m. This will allow us to
# examine the boundary layer profiles at these locations.
# Have looked in viscous_data.dat file to determine which 
# i-index to select.
e3post.py --job=$JOB --output-file=turb_flat_plate-x-368mm.dat --tindx=$TINDX \
    --slice-list="22:23,-6,:,0" --add-pitot-p --add-mach

# Extracts a slice of the nearest cells (to the wall) along 
# the plate. This will allow us to examine viscous properties,
# like skin friction and y+ values.
e3post.py --job=$JOB --output-file=turb_flat_plate-y-wall.dat --tindx=$TINDX \
    --slice-list="1,:,-1,:;3,:,-1,:;5,:,-1,:;7,:,-1,:;9,:,-1,:;11,:,-1,:;13,:,-1,:;15,:,-1,:;17,:,-1,:;19,:,-1,:;21,:,-1,:;23,:,-1,:" \
    --add-pitot-p --add-mach

# Runs python script to compute viscous parameters.
echo "Computing viscous parameters .."
python compute_viscous_data_simple.py

# Compute boundary layer thickness based on 0.99 of freestream velocity
echo "Estimating boundary layer thickness based on 0.99 of freestream velocity .."
python estimate-bl-thickness.py

# Compute displacement and momentum thicknesses
echo "Estimating displacement and momentum thicknesses .."
awk -f integral-thicknesses.awk turb_flat_plate-x-368mm.dat 

# Plots numerical results against experimental data.
echo "Plotting results .."
./plot.sh

echo "Postprocessing completed."
