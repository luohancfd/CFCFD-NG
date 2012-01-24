#! /bin/sh
for D in mdot_0.000 mdot_0.002 mdot_0.009 
do
cd $D
./clean.sh
cd ..
done
