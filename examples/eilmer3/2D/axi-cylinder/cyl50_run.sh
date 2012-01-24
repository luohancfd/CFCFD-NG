#! /bin/sh
# cyl50_run.sh
e3prep.py --job=cyl50 --do-svg
time e3shared.exe --job=cyl50 --run
e3post.py --job=cyl50 --tindx=2 --vtk-xml

# Extract the profile near the downstream end of the cylinder.
e3post.py --job=cyl50 --tindx=2 --slice-list="0,47,:,0"

echo "At this point, we should have a new solution"
echo "Run cyl50_plot.sh next"

