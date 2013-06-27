#! /bin/bash
# prep.sh
e3prep.py --job=dbl-cone --do-svg
e3loadbalance.py --job=dbl-cone -n 14
e3post.py --job=dbl-cone --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

