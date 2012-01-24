#!/bin/sh
#  back_run.sh
#  Exercise the Navier-Stokes solver for the conical nozzle 
#  as used by Back, Massier and Gier (1965) AIAA J. 3(9):1606-1614.
e3prep.py --job=back --do-svg
time e3shared.exe --job=back --run
e3post.py --job=back --tindx=9999 --vtk-xml
