#! /bin/bash -I
# run_march.sh
echo "Start time: "; date
e3march.py --job=lam_flat_plate_march --run --nbj=4 > LOGFILE_march
echo "Finish time: "; date
