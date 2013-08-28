#! /bin/bash
# prep.sh
e3prep.py --job=convex-ramp --do-svg
e3loadbalance.py --job=convex-ramp -n 4
e3post.py --job=convex-ramp --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

