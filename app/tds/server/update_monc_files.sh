#!/bin/sh
# update_monc_files.sh
#
# Looks in the (remote) UQSPACE directories and updates
# the local copy of the MONC data files for each facility.
# Once the files have been updated and tidied up, 
# it runs the td_scan.tcl script.
#
# PJ, 28-Apr-03
#

SRC=/mnt/uqspace
DEST=/data/tunnels
FACILITY_LIST="T4 X2 X3 X1"

# --------------------------------------------------------
echo "Stage 1: copy new files from UQSPACE backup area."

if ! [ -d $SRC/T4/DATA ]
then
    echo "$SRC directories cannot be seen."
    echo "Use ncpmount to make them accessible."
    exit -1
else
    echo "$SRC directories are visible, let us do some work."
fi

if ! [ -d $DEST ]
then
    echo "Cannot see the destination directory $DEST."
    exit -1
else
    echo "Can also see the destination directory $DEST."
fi

echo "T4:"
rsync -av $SRC/T4/DATA/ $DEST/T4/
echo "X3:"
rsync -av $SRC/XTUBE/X3/ $DEST/X3/
echo "X2:"
rsync -av $SRC/XTUBE/X2/ $DEST/X2/
echo "X1 General:"
rsync -av $SRC/XTUBE/X1/ $DEST/X1/
echo "X1 Sam Chiu Series 13:"
rsync -av $SRC/XTUBE/X1/CHIU/SER13/ $DEST/X1/
echo "X1 Dillon Hunt:"
rsync -av $SRC/XTUBE/X1/DILLON/ $DEST/X1/

# -------------------------------------------------------
echo "Stage 2: Clean-up."

CRAP_LIST="T4/5848/5849"
for CRAP in $CRAP_LIST
do
    echo "    Checking for $DEST/$CRAP"
    if [ -d $DEST/$CRAP ]
    then
	echo "    Deleting $DEST/$CRAP"
	rm -rf $DEST/$CRAP
    fi
done

echo "Copy some of the description files to lower-case directories."
rsync -av $DEST/T4/RUNDESC/ $DEST/T4/rundesc/
rsync -av $DEST/T4/DESCRIPT/ $DEST/T4/descript/
rsync -av $DEST/X1/RUNDESC/ $DEST/X1/rundesc/
rsync -av $DEST/X1/DESCRIPT/ $DEST/X1/descript/
rsync -av $DEST/X2/RUNDESC/ $DEST/X2/rundesc/
rsync -av $DEST/X2/DESCRIPT/ $DEST/X2/descript/
rsync -av $DEST/X3/RUNDESC/ $DEST/X3/rundesc/
rsync -av $DEST/X3/DESCRIPT/ $DEST/X3/descript/

# -------------------------------------------------------
echo "Stage 3: scanning for new MONC files in the archive."

if [ -f td_scan.tcl ] 
then
    TODAYS_DATE="`date +%d%b%y`"
    echo "Todays Date = $TODAYS_DATE"
    for FACILITY in $FACILITY_LIST
    do
        echo "    Scanning for Facility $FACILITY"
        ./td_scan.tcl $DEST/$FACILITY/ new
        mv td_scan.log $DEST/$FACILITY/td_scan_$TODAYS_DATE.log
        ./td_scan.tcl $DEST/$FACILITY/ index
    done
else
    echo "Cannot see the td_scan.tcl script in the current directory."
    exit -1
fi
