#!/bin/bash
# prepare-grid.sh
e3prep.py --job=the_minimal_grid --do-svg
e3post.py --job=the_minimal_grid --tindx=0 --vtk-xml
e3prep.py --job=the_plain_bottle --do-svg
e3post.py --job=the_plain_bottle --tindx=0 --vtk-xml
e3prep.py --job=the_clustered_bottle --do-svg
e3post.py --job=the_clustered_bottle --tindx=0 --vtk-xml

