#! /bin/bash 
# run-block-marching.sh
#
e3prep.py --job=bump --clean-start > LOGFILE_PREP
echo "Start time: "; date
e3march.py --job=bump --nbj=2 --run  > LOGFILE_RUN_MARCHING
echo "Finish time: "; date
e3post.py --job=bump --tindx=1 --add-mach --vtk-xml
