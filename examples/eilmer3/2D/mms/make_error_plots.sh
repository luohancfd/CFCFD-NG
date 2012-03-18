#!/bin/bash
echo "Make VTK plot showing simulation error:"
e3post.py --job=mms --tindx=20 \
    --gmodel-file="very-viscous-air.lua" \
    --ref-function=analytic_solution_wrapper.py \
    --vtk-xml

