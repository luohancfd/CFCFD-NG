#!/bin/bash -l

e3prep.py --job=gcl-2d

e3shared.exe --job=gcl-2d --run

e3post.py --job=gcl-2d --vtk-xml --moving-grid --tindx=all
