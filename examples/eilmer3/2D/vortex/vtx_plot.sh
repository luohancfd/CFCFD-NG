#! /bin/sh
# vtx_plot.sh

# Generate the ideal profile.
awk -f make_profile.awk

# Extract the flow data 45 degrees around.
e3post.py --job=vtx --tindx=9999 --output-file=vtx_profile_45.dat \
    --slice-list="0,39,:,0"

# Extract the flow data 90 degrees around.
e3post.py --job=vtx --tindx=9999 --output-file=vtx_profile_90.dat \
    --slice-list="0,-1,:,0"

awk -f extract_radial.awk vtx_profile_45.dat > radial_profile_45.dat
awk -f extract_radial.awk vtx_profile_90.dat > radial_profile_90.dat

# Generate postscript plots of the radial profiles.
gnuplot radial_profile.gnu

echo At this point, we should have a plotted the solution

