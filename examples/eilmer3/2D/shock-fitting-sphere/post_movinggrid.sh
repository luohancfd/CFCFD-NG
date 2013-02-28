#!/bin/bash
# post_movinggrid.sh

stages="0 1 2 3"
for STAGE in ${stages} 
do
    echo "Stage $STAGE:"

    e3post.py --job=sphere${STAGE} --tindx=all --vtk-xml --add-mach --moving-grid
        
done

