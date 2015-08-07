#! /bin/bash -I
# run_march.sh

echo "Start time: "; date
e3march.py --job=t4-m4-nozzle --run --nbj=2 > LOGFILE
echo "Finish time: "; date
