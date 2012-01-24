#! /bin/sh
# colesBL_post.sh

JOB=colesBL
TINDX=5

e3post.py --job=$JOB --vtk-xml --tindx=$TINDX

e3post.py --job=$JOB --output-file=$JOB-ix21.dat --tindx=$TINDX \
    --slice-list="3,21,:,:" --add-pitot-p
e3post.py --job=$JOB --output-file=$JOB-ix3.dat --tindx=$TINDX \
    --slice-list="3,3,:,:" --add-pitot-p

awk -f integral-thicknesses.awk $JOB-ix3.dat > thicknesses-ix3.txt
awk -f integral-thicknesses.awk $JOB-ix21.dat > thicknesses-ix21.txt
