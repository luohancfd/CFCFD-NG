#!/bin/bash
# run_sumulation.sh
# catch both stdout and stderr
nohup time e3shared.exe --job=bar --run &> LOGFILE &
