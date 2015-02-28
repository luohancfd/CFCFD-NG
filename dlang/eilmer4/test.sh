#!/bin/bash
# test.sh
./e4shared --prep --job=test
./e4shared --run --job=test --verbosity=1
cp test.list block_labels.list
e3post.py --job=test --tindx=all --vtk-xml
rm block_labels.list
