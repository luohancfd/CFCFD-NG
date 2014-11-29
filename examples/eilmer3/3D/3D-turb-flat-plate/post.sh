#! /bin/sh
# post.sh

JOB=turb_flat_plate

e3post.py --job=$JOB --vtk-xml --add-mach --add-pitot-p --tindx=all

#Plate-Surface Data (Cf)

e3post.py --job=$JOB --output-file=$JOB-y-wall0.dat --slice-list="1,:,0,-1;5,:,0,-1;9,:,0,-1;13,:,0,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall1.dat --slice-list="1,:,1,-1;5,:,1,-1;9,:,1,-1;13,:,1,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall2.dat --slice-list="1,:,2,-1;5,:,2,-1;9,:,2,-1;13,:,2,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall3.dat  --slice-list="1,:,3,-1;5,:,3,-1;9,:,3,-1;13,:,3,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall4.dat --slice-list="1,:,4,-1;5,:,4,-1;9,:,4,-1;13,:,4,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall5.dat --slice-list="3,:,0,-1;7,:,0,-1;11,:,0,-1;15,:,0,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall6.dat --slice-list="3,:,1,-1;7,:,1,-1;11,:,1,-1;15,:,1,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall7.dat --slice-list="3,:,2,-1;7,:,2,-1;11,:,2,-1;15,:,2,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall8.dat --slice-list="3,:,3,-1;7,:,3,-1;11,:,3,-1;15,:,3,-1" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-y-wall9.dat --slice-list="3,:,4,-1;7,:,4,-1;11,:,4,-1;15,:,4,-1" --add-pitot-p --add-mach

#End-Slice Data (Mach Number, Temp, Pressure)

e3post.py --job=$JOB --output-file=$JOB-x-368mm0.dat --slice-list="12:13,26,0,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm1.dat --slice-list="12:13,26,1,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm2.dat --slice-list="12:13,26,2,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm3.dat --slice-list="12:13,26,3,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm4.dat --slice-list="12:13,26,4,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm5.dat --slice-list="14:15,26,0,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm6.dat --slice-list="14:15,26,1,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm7.dat --slice-list="14:15,26,2,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm8.dat --slice-list="14:15,26,3,:" --add-pitot-p --add-mach
e3post.py --job=$JOB --output-file=$JOB-x-368mm9.dat --slice-list="14:15,26,4,:" --add-pitot-p --add-mach

python compute_viscous_data_master.py

./plot.sh
./plot2Dv3D.sh
