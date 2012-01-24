#! /bin/sh
for D in mdot_0.000 mdot_0.004 mdot_0.010 mdot_0.013
do
cd $D
./clean.sh
cd ..
done
