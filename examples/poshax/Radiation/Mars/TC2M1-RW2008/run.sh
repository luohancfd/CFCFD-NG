#!/bin/bash
# Comparison of TC2M1 models with EAST data

# 0. Clear working directories
rm -R -f work-dir/*
# cp -R results results-saved
# cp -R plots plots-saved
# rm -R -f results/*
rm -R -f plots/*

# 1.   Gas-model loop
for GM in 2T-basic-Lee 2T-full-Lee 2T-full-Park 2T-full-PL 2T-full-PLF 2T-full-PLFM 3T-full-PLF
do
cp $GM/ParkX.inp work-dir/ParkX.inp
cp $GM/ParkX.py work-dir/ParkX.py

#   b. Radiation feature loop
for RF in VUV C-line CN-Violet-A CN-Violet-B UV C2-Swan-A C2-Swan-B O-line
do
cp radiation-files/CNO-radiators-$RF.inp work-dir/CNO-radiators.inp
cp cfg_files/EAST_$RF.cfg work-dir/EAST.cfg
cd work-dir/
script_noneq.py -i ParkX.py -p mTg.pgm -r ParkX.rsi
poshax.x EAST.cfg
mv Irad_smeared_V_x.dat ../results/$GM-$RF\_IvX.dat
cd ..
# end rad loop
done
mv work-dir/TC2M1.data results/TC2M1-$GM.data
# end gas loop
done

# 2. Special cases: coupled radiation
for GM in 2T-full-PLF
do
cp $GM/ParkX.inp work-dir/ParkX.inp
cp $GM/ParkX.py work-dir/ParkX.py

#   b. NEQ Radiation feature loop
for RF in VUV C-line CN-Violet-A CN-Violet-B UV C2-Swan-A C2-Swan-B O-line
do
cp radiation-files/CNO-radiators-$RF.inp work-dir/CNO-radiators.inp
cp cfg_files/EAST_$RF-RC.cfg work-dir/EAST.cfg
cd work-dir/
script_noneq.py -i ParkX.py -p mTg.pgm -r ParkX.rsi
poshax.x EAST.cfg
mv Irad_smeared_V_x.dat ../results/$GM-RC-NEQ-$RF\_IvX.dat
cd ..
# end rad loop
done

#   b. EQ Radiation feature loop
for RF in VUV C-line CN-Violet-A CN-Violet-B UV C2-Swan-A C2-Swan-B O-line
do
cp radiation-files/CNO-radiators-$RF-EQ.inp work-dir/CNO-radiators.inp
cp cfg_files/EAST_$RF-RC.cfg work-dir/EAST.cfg
cd work-dir/
script_noneq.py -i ParkX.py -p mTg.pgm -r ParkX.rsi
poshax.x EAST.cfg
mv Irad_smeared_V_x.dat ../results/$GM-RC-EQ-$RF\_IvX.dat
cd ..
# end rad loop
done

# end gas loop
done

# 3. plot results
rm work-dir/*
cp results/* plot_files/* EAST_data/* work-dir/
cd work-dir/

for RF in VUV C-line CN-Violet-A CN-Violet-B UV C2-Swan-A C2-Swan-B O-line
do
gnuplot plot_$RF-intensities.gnu
done

gnuplot plot_basics.gnu

mv *.eps ../plots/

echo "Finished running through all model variations."
