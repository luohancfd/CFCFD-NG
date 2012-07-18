#!/bin/bash
echo "Check that simulation has converged by comparing solution instances:"
e3post.py --job=mms --tindx=6 --gmodel-file="very-viscous-air.lua" \
    --compare-job=mms --compare-tindx=20
e3post.py --job=mms --tindx=7 --gmodel-file="very-viscous-air.lua" \
    --compare-job=mms --compare-tindx=20

echo "----------------------------------------------------------------------"
echo "Check simulation against analytical data:"
e3post.py --job=mms --tindx=20 --gmodel-file="very-viscous-air.lua" \
    --ref-function=analytic_solution_wrapper.py \
    --global-norm-list="rho,L2"

echo "----------------------------------------------------------------------"
echo "Generate VTK files for plotting:"
e3post.py --job=mms --tindx=20 --gmodel-file="very-viscous-air.lua" \
    --vtk-xml

