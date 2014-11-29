#!/bin/bash
# post.sh

e3post.py --job=rfp --vtk-xml --tindx=all --add-mach

e3post.py --job=rfp --vtk-xml --tindx=last --add-mach \
          --slice-list="0,:,0,0;2,:,0,0" --output-file="surface.data"

