#!/bin/sh
# run.sh
e3prep.py --job=t4noz --do-svg > LOGFILE_PREP
e3shared.exe --job=t4noz --run > LOGFILE_RUN
e3post.py --job=t4noz --tindx=9999 --vtk-xml > LOGFILE_POST
e3post.py --job=t4noz --tindx=9999 --output-file=t4noz-exit.data \
    --slice-list="19,-2,:,0" --add-pitot-p --add-mach --add-total-enthalpy \
    --add-total-p > LOGFILE_SLICE
./stats.py t4noz-exit.data > LOGFILE_STATS
