#!/bin/sh
# post_simulation.sh

e3post.py --job=titan_x2_shell --vtk-xml

e3post.py --job=titan_x2_shell --output-file=titan-surface \
    --surface-list="0,BOTTOM;1,TOP;2,BOTTOM;3,BOTTOM;4,BOTTOM;5,BOTTOM;6,TOP;7,TOP;8,BOTTOM;9,BOTTOM;10,TOP;11,BOTTOM;12,BOTTOM"

echo "Done."