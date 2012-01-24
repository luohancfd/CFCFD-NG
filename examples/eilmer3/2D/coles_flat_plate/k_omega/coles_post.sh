#! /bin/sh
# coles_post.sh

e3post.py --job=coles --vtk-xml --tindx=5

e3post.py --job=coles --output-file=coles-ix21.dat --tindx=5 \
    --slice-list="3,21,:,:" --add-pitot-p
e3post.py --job=coles --output-file=coles-ix3.dat --tindx=5 \
    --slice-list="3,3,:,:" --add-pitot-p

awk -f integral-thicknesses.awk coles-ix3.dat > thicknesses-ix3.txt
awk -f integral-thicknesses.awk coles-ix21.dat > thicknesses-ix21.txt
