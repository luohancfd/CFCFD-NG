#!/bin/bash
# PJ, 26-Feb-2011
# Just enough to demonstrate that we can load the grids
# and generate the flow field.
e3prep.py --job=jones-stator
e3post.py --job=jones-stator --tindx=0 --vtk-xml
