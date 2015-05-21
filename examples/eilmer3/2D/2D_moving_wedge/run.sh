#!/bin/bash -l

e3prep.py --job=plate_2d --do-svg

e3shared.exe --job=plate_2d --run

e3post.py --job=plate_2d --vtk-xml --moving-grid --tindx=all
