#!/bin/sh
# prepare_grid.sh

gzip -d icem_grid_plot3d.fmt
import_grid.py --input=icem_grid_plot3d.fmt \
               --output=icem_grid \
               --plot3dplanes
gzip icem_grid_plot3d.fmt

echo "Done."