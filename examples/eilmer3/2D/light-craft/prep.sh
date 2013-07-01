#! /bin/bash
# prep.sh
e3prep.py --job=lc --do-svg
e3loadbalance.py --job=lc -n 4
e3post.py --job=lc --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

