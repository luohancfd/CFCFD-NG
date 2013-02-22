#!/bin/sh
# scan_monc_files.sh
#
# Runs the td_scan.tcl script over the collection of monc files.
#
# PJ, 28-Apr-03, copied out of update_monc_files.sh
#

SRC=/mnt/uqspace
DEST=/data/tunnels
FACILITY_LIST="T4 X2 X3 X1"

# --------------------------------------------------------

if ! [ -d $DEST ]
then
    echo "Cannot see the destination directory $DEST."
    exit -1
else
    echo "Can also see the destination directory $DEST."
fi

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
