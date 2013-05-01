#! /bin/bash
# prep.sh
e3prep.py --job=swlbli --do-svg
e3loadbalance.py --job=swlbli -n 4
e3post.py --job=swlbli --tindx=0 --vtk-xml

echo "At this point, we should have a grid."
echo "Use run.sh next"

