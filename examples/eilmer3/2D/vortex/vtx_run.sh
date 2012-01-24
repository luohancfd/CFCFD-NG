#! /bin/sh
# vtx_run.sh
e3prep.py --job=vtx --do-svg
time e3shared.exe --job=vtx --run
e3post.py --job=vtx --vtk-xml
