#! /bin/sh
# ram2_prep.sh

e3prep.py --job=ram2 --do-svg
e3post.py --job=ram2 --tindx=0 --vtk-xml

echo "At this point, we should have a starting grid"
echo "Run ram2_run.sh next"

