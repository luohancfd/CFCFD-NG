#!/bin/bash
#
# Radiation argon shock layer test case.
#
# DFP, 2-June-2014

echo "Run the inviscid stage"
cd part1-inviscid/
e3prep.py --job=hemisphere > LOGFILE_PREP
mpirun -np 4 e3mpi.exe -f hemisphere -r > LOGFILE_RUN
cd ..

echo "Run the viscous stage"
cd part2-viscous/
e3prep.py --job=hemisphere > LOGFILE_PREP
echo "Adding viscous effects"
mpirun -np 4 e3mpi.exe -f hemisphere -q -r > LOGFILE_RUN
echo "Increasing CFL number to 0.5"
set_control_parameter.py hemisphere.control cfl 0.5
set_control_parameter.py hemisphere.control max_time 1.88964e-05
mpirun -np 4 e3mpi.exe -f hemisphere -t 1 -q -r >> LOGFILE_RUN

echo "Run the viscous with radiation stage"
cd part3-viscous-with-radiation/
radmodel.py -i Ar-nonequilibrium-radiation.py -L rad-model.lua > LOGFILE_PREP
e3prep.py --job=hemisphere >> LOGFILE_PREP

echo "Run e3mpi for one body length on new grid"
mpirun -np 4 e3mpi.exe -f hemisphere -q -r > LOGFILE_RUN
get_residuals.py 0 residuals-0.txt
mv e3mpi*.log e3mpi.log.part1/
echo "First radiation transport calculation"
e3rad.exe -f hemisphere -q -t 1 -r > LOGFILE_RUN

echo "Run e3mpi for another body length with radiation coupling"
set_control_parameter.py hemisphere.control max_time 7.55855e-06
mpirun -np 4 e3mpi.exe -f hemisphere -q -t 2 -r >> LOGFILE_RUN
get_residuals.py 0 residuals-1.txt
echo "Final radiation transport calculation"
radmodel.py -i Ar-nonequilibrium-radiation-180to6000nm.py -L rad-model.lua >> LOGFILE_PREP
e3rad.exe -f hemisphere -q -t 3 -r >> LOGFILE_RUN

echo "Extract final surface heat flux profile"
e3post.py --job=hemisphere --tindx=4 --heat-flux-list="2:3,1,:,:,:" > LOGFILE_POST
cd ..

# Compute the radiative heat flux error with respect to the experiment
# measurement using average of the two experimental datapoints with
# blackened gauges at Ms=12.7 from Figure 9
./compute_qrad_error.py part3-viscous-with-radiation/hf_profile.data 5.5188e7 > LOGFILE_COMPARE
