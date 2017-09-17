#! /bin/sh
# run.sh

l_script.py -f valve
l1d.exe -f valve -prep


echo 
echo Beginning simulation...
echo Console output is being caught in file valve.log.
time l1d.exe -f valve > valve.log
echo Finished simulation.
echo 

