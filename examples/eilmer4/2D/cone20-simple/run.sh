#!/bin/bash
# run.sh
e4shared --prep --job=cone20
e4shared --run --job=cone20 --verbosity=1
cp cone20.list block_labels.list
e3post.py --job=cone20 --tindx=all --vtk-xml
rm block_labels.list
