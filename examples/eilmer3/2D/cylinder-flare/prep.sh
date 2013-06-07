#! /bin/bash
# prep.sh
e3prep.py --job=cyl-flare --do-svg
e3loadbalance.py --job=cyl-flare -n 12
e3post.py --job=cyl-flare --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

