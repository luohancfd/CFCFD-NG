# testcase_full_run.sh
# 
# Before runnig setup the testcase_setup.py file

# Prepare the simulation input files (parameter, grid and initial flow data).
e3prep.py --job=testcase --do-vrml # do-vrml : render 3D blocks to a virtual

# Integrate the solution in time, to produce flow data at subsequent time
e3shared.exe --job=testcase --run 

# Extract the flow solution data (to produce one unstructured grid) and reformat.
e3post.py --job=testcase --vtk-xml --tindx=all --omegaz="[0.0, 2000.0, 0.0]"
#e3post.py --job=testcase --vtk-xml --tindx=all --omegaz="[2000.0]"
#e3post.py --job=testcase --vtk-xml --tindx=all --omegaz="[0.0, 0.0, 0.0, 0.0, 2000.0, 2000.0, 2000.0, 2000.0]"
#e3post.py --job=testcase --vtk-xml --tindx=all --omegaz="[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]"
# Analyse and save data (pressure and velocity at MP) with testcase_output.py
#./testcase_output.py
