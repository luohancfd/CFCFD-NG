#!/bin/bash -l
e3prep.py --job=gcl-3d

e3shared.exe --job=gcl-3d --run

e3post.py --job=gcl-3d --vtk-xml --moving-grid --tindx=all
