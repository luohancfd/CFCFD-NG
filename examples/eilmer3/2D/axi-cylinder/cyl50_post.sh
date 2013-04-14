#! /bin/sh
# cyl50_post.sh

e3post.py --job=cyl50 --tindx=2 --vtk-xml

# Extract the profile near the downstream end of the cylinder.
e3post.py --job=cyl50 --tindx=2 --slice-list="2,22,:,0;3,22,:,0"

echo "At this point, we should have a new solution"
echo "Run cyl50_plot.sh next"

