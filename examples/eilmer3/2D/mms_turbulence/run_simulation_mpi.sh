#!/bin/bash                                                                                                                          
# 8X8_grids_2x2_blocks_run_mpi.sh                                                                                                                 
# do the order of accuracy testing for Eilmer turbulence model under multiple blocks.                                                                                                                                                            

# Prepare the simulation input files (parameter, grid and initial flow data).                                                       
# The SVG file provides us with a graphical check on the geometry   
python make_source_terms.py
cp mms-regular.py mms.py
e3prep.py --job=mms

# Integrate the solution in time using the parallel processing code. 
time mpirun -np 4 e3mpi.exe --job=mms --run

