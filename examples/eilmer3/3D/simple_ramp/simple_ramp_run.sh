#! /bin/sh
# simple_ramp_run.sh
e3prep.py --job=simple_ramp
time e3shared.exe --job=simple_ramp --run --verbose
e3post.py --job=simple_ramp --vtk-xml --tindx=all
