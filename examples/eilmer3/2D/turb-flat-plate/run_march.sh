#! /bin/bash -I
# run_march.sh
echo "Start time: "; date
e3march.py --job=turb_flat_plate_march --run --nbj=2 > LOGFILE_march
echo "Finish time: "; date
