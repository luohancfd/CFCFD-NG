#!/bin/bash
# post.sh
# Extract the flowfield data along the stagnation line of cells and one beside.
e3post.py --job=hayabusa --tindx=9999 --output-file=line0.data \
    --slice-list="0,:,0,0;3,:,0,0;6,:,0,0" --add-pitot-p --add-mach
e3post.py --job=hayabusa --tindx=9999 --output-file=line1.data \
    --slice-list="0,:,1,0;3,:,1,0;6,:,1,0" --add-pitot-p --add-mach