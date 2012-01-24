#! /bin/sh
# kiock_post.sh

# Around the blade, suction surface
e3post.py --job=kiock --output-file=simulation_results.dat --tindx=9999 \
    --slice-list="0,:,0,0; 1,:,0,0; 2,:,0,0; 3,-1,:,0; 4,:,-1,0; 5,:,-1,0"
#                 0south   1south   2south   3east     4north    5north

# Extract the solution data over whole flow domain and reformat.
e3post.py --job=kiock --vtk-xml --add-mach

# Calculate average flow properties at inlet and outlet
turbo_post.py kiock

# Python script to convert x-coordinate surface Mach to chord position Mach number
python plot_surface_Mach.py

echo At this point, we should have a plotted data.

