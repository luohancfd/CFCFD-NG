#! /bin/bash
# prep.sh
e3prep.py --job=cubic-ramp --do-svg
e3loadbalance.py --job=cubic-ramp -n 4
e3post.py --job=cubic-ramp --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

