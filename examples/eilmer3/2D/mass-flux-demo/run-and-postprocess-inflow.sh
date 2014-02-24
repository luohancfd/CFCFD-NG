e3prep.py --job=mfd-inflow
e3shared.exe --job=mfd-inflow --run
e3post.py --job=mfd-inflow --tindx=all --add-mach --vtk-xml