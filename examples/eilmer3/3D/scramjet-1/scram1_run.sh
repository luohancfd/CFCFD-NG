#! /bin/sh
# scram1_run.sh
e3prep.py --job=scram1
time e3shared.exe --job=scram1 --run
e3post.py --job=scram1 --vtk-xml
