#!/bin/bash

# 1. Output solution in vtk-xml format
e3post.py --job=hayabusa --tindx=9999 --vtk-xml

# 2. Slice along the stagnation line
e3post.py --job=hayabusa --tindx=9999 --slice-list="0,:,0,:;2,:,0,:"

# 3. Heat flux profile along vehicle surface
e3post.py --job=hayabusa --tindx=9999 --heat-flux-list="2:3,1,:,:,:"
