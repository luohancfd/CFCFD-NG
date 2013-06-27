#PBS -S /bin/bash
#PBS -l walltime=09:00:00
#PBS -l select=1:ncpus=14
#PBS -W group_list=seses
#PBS -N CIC_CO_CH

## -j oe : standard error and standard output are merged into standard output.
#PBS -j oe
## Send an email when a job is *a*borted by the batch system, 
## when a job *b*egins execution and when a job t*e*rminates.
#PBS -m ae
## Send email information to the following address
#PBS -M nikhil_banerji@epfl.ch
#

echo "----------------------------------"
echo "Begin MPI job..."
cd $PBS_O_WORKDIR

mkdir outfiles_CO-CH

cluster_CIC.py --species='CO CH' --T='300' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_300 &
cluster_CIC.py --species='CO CH' --T='500' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_500 &
cluster_CIC.py --species='CO CH' --T='1000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_1000 &
cluster_CIC.py --species='CO CH' --T='2000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_2000 &
cluster_CIC.py --species='CO CH' --T='4000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_4000 &
cluster_CIC.py --species='CO CH' --T='5000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_5000 &
cluster_CIC.py --species='CO CH' --T='8000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_8000 &
cluster_CIC.py --species='CO CH' --T='11000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_11000 &
cluster_CIC.py --species='CO CH' --T='12000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_12000 &
cluster_CIC.py --species='CO CH' --T='15000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_15000 &
cluster_CIC.py --species='CO CH' --T='16000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_16000 &
cluster_CIC.py --species='CO CH' --T='20000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_20000 &
cluster_CIC.py --species='CO CH' --T='26000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_26000 &
cluster_CIC.py --species='CO CH' --T='32000' --usr='nbanerji' --reactype='mm' > outfiles_CO-CH/outfile_CO-CH_32000 &
wait 

echo "End MPI job."
date
