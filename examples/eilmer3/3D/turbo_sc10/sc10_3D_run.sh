#! /bin/sh
# sc10_3D_run.sh

nohup e3shared.exe --job=sc10_3D --run --zip-files 1> LOGFILE 2>&1 &

echo "Simulation should now be running in background"
echo "with STDOUT and STDERR redirected to LOGFILE."
