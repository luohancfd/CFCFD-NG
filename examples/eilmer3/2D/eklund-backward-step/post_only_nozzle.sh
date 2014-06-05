job=backstep
tindx=9999

#e3post.py --job=$job --tindx=$tindx --vtk-xml

# To get profiles for checking inflow into backward facing step
e3post.py --job=$job --tindx=$tindx --zip-files --output-file=xh_-10.dat --slice-list="0,0,:,:;1,0,:,:"
