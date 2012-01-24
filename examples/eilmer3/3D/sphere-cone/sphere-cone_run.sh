#! /bin/bash
# sphere-cone_run.sh
e3prep.py --job=sphere-cone
time e3shared.exe --job=sphere-cone --run

e3post.py --job=sphere-cone --vtk-xml
e3post.py --job=sphere-cone --surface-list="0,EAST" --output-file="sphere-cone-surface"
