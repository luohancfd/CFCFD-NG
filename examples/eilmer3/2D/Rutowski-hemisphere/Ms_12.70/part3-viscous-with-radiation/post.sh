#!/bin/bash 

# 0. Set the tindx to work on
TINDX=10024

# 1. Extract vtk solution
e3post.py --job=hemisphere --tindx=$TINDX --vtk-xml

# 2. Extract stagnation streamline profile
e3post.py --job=hemisphere --tindx=$TINDX --slice-list="0,:,0,:;3,:,0,:;6,:,0,:;9,:,0,:"

# 3. Perform tangent-slab calculation along the stagnation stream line
#radmodel.py -i argon-radiators.py -L rad-model.lua
#e3post.py --job=hemisphere --tindx=$TINDX --tangent-slab-list="0,:,0,:;3,:,0,:;6,:,0,:;9,:,0,:"

# 3. Extract surface heating profile
e3post.py --job=hemisphere --tindx=$TINDX --heat-flux-list="9:11,1,:,:,:"
