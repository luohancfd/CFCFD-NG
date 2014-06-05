job=backstep
tindx=0005   #9999

# xh_-10 = -3.175e-3 m
# xh_-01 = -0.3175e-3 m
# xh_17 = 5.3975e-3 m
# xh_30 = 9.525e-3 m
# xh_39 = 12.3825e-3 m
# xh_67 = 21.2725e-3 m
# xh_108 = 34.29e-3 m


# vtk or visit files for visualisation
e3post.py --job=$job --tindx=$tindx --vtk-xml

# Boundary layer slices to compare with experimental data
e3post.py --job=$job --tindx=$tindx --output-file=xh_-10.dat --slice-list="15,0,:,:;16,0,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_-01.dat --slice-list="15,-5,:,:;16,-5,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_17.dat --slice-list="1,21,:,:;7,21,:,:;8,21,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_30.dat --slice-list="2,14,:,:;9,14,:,:;10,14,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_39.dat --slice-list="2,25,:,:;9,25,:,:;10,25,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_67.dat --slice-list="3,24,:,:;11,24,:,:;12,24,:,:"
e3post.py --job=$job --tindx=$tindx --output-file=xh_108.dat --slice-list="4,-1,:,:;13,-1,:,:;14,-1,:,:"

