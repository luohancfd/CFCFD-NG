#! /bin/sh
# sc10_3D_post.sh

e3post.py --job=sc10_3D --zip-files --tindx=all --add-mach --add-total-p \
    --surface-list="0,south;1,east;2,north;3,north;0,bottom;1,bottom;2,bottom;3,bottom;\
4,bottom;5,bottom;6,bottom;7,bottom;8,bottom;9,bottom;10,bottom;11,bottom;12,bottom;13,bottom"

e3post.py --job=sc10_3D --zip-files --tindx=all --add-mach --add-total-p --vtk-xml

