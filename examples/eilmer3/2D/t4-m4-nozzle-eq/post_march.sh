#! /bin/sh
# post_march.sh

JOB=t4-m4-nozzle
TINDX=0001
GMODEL=cea-lut-air.lua.gz

# Exit plane slice.
e3post.py --job=$JOB --tindx=$TINDX --gmodel-file=$GMODEL \
          --output-file=nozzle-exit.data --slice-list='-2:-1,-2,:,0' \
          --add-mach --add-pitot --add-total-enthalpy --add-total-p

# Generate files for plotting with Paraview.
e3post.py --job=$JOB --vtk-xml --tindx=$TINDX --gmodel-file=$GMODEL \
          --add-mach --add-pitot --add-total-enthalpy --add-total-p

# Compute viscous data at the nozzle wall.
nenzfr_compute_viscous_data.py --job=$JOB --nbj=2
