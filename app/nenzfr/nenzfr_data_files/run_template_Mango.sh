# This file should be modified as required however the --pe and --job options must be left as --pe=$pe
# and --job=$jobName respectively.
# Luke Doherty
# 24-02-2012


nohup nenzfr.py --gas="air5species" --T1=300.0 --p1=250000.0 --Vs=2195.0 --pe=$pe --chem="neq" --area=1581.165 --job=$jobName --cfile="contour-t4-m10.data" --gfile="None" --exitfile="nozzle-exit.data" --nni=2700 --nnj=150 --nbi=270 --bx=1.05 --by=1.002 --max-time=0.006 --max-step=800000 > LOGFILE &




